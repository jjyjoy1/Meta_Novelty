### 1  |  Key factors to bake into any **paired metagenome ↔ metaproteome** workflow

| Factor to plan for                      | Why it matters                                                                                           | Practical guard-rails                                                                                                                                    |
| --------------------------------------- | -------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Sample pairing metadata**             | All statistics must respect the “mates” (e.g. tumour vs adjacent normal from same patient, or T0 vs T1). | Keep a `subject_id` field in every table and use *GroupKFold* or mixed-effects models so that both halves of a pair always live in the **same** CV fold. |
| **Compositionality differences**        | Metagenome = read counts; metaproteome = spectral intensities → different scale & zero structure.        | CLR-transform taxa/gene tables; variance-stabilising normalisation (VSN) or MaxLFQ for proteins, then optionally re-CLR if you merge them.               |
| **Database completeness & consistency** | MAG-derived proteins + external refs must share the **same taxonomy backbone** or mappings will clash.   | Build one custom FASTA from your MAG ORFs; label entries with GTDB R207 IDs or SGBs (MetaPhlAn 4) for both gene and protein layers. ([Nature][1])        |
| **Small / lineage-specific genes**      | Many peptides come from short ORFs missed by generic predictors.                                         | Run a lineage-aware gene caller (e.g. LS-Gene) before ORF→protein export. ([PMC][2])                                                                     |
| **Mapping key between layers**          | You need a reversible link MAG → contig → ORF → peptide.                                                 | Store ORF IDs in the FASTA header and keep the mapping table in SQLite/Parquet; most M2P pipelines will propagate it automatically. ([PMC][3])           |
| **Zero inflation / missingness**        | Proteins may be undetected even when the gene copy exists.                                               | Treat “missing” as *censored* (use hurdle or Tobit models) or carry a binary “detected / not detected” feature alongside intensity.                      |
| **Batch effects**                       | Library prep and LC–MS runs introduce layer-specific noise.                                              | Include batch covariates in every linear model **and** feed the same covariate matrix to multi-omics tools so correction is coherent.                    |
| **Phylogeny vs taxonomy**               | GTDB/NCBI labels sometimes disagree with tree-based placement.                                           | For diversity, compute distances on a concatenated marker-gene tree (PhyloPhlAn / pplacer) and treat taxonomy as a loose label. ([PMC][4])               |

---

### 2  |  When **gene, protein, taxonomy or phylogenetic calls disagree**, walk through this decision tree

| Inconsistency observed                                                                                      | Typical root causes                                                                                                                       | Triage & reconciliation                                                                                                                                                                                                                 |
| ----------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Gene present, protein absent**                                                                            | (a) Low transcription/translation<br>(b) Short ORF missed by proteomics engine<br>(c) Split MAG scaffold → peptide DB lacks full sequence | • Check RNA (if available) to see if gene is transcribed.<br>• Re-search spectra with **open search + semi-tryptic** to catch truncated peptides.<br>• Verify gene integrity: re-assemble with long-reads if coverage is low.           |
| **Protein detected, gene copy not found**                                                                   | (a) Strain variant not captured by MAG bin<br>(b) Peptide shared across paralogs and mis-assigned<br>(c) Contamination / carry-over       | • Re-map peptide to the entire contig set (not just MAGs).<br>• Lower FDR to 0.5 % and require ≥ 2 unique peptides.<br>• BLAST peptide against public db; if best hit < 95 % identity call it an *orphan* and flag for manual curation. |
| **Taxonomy mismatch between layers**<br>(e.g. gene says *Bacteroides*, protein search returns *Prevotella*) | (a) Out-of-date refs, different nomenclatures<br>(b) LCA algorithm depth thresholds differ<br>(c) Horizontal gene transfer                | • Harmonise both layers on **GTDB** species-level bins (or SGBs).<br>• Use a majority-vote or confidence-weighted consensus; if disagreement persists, collapse to the **lowest common ancestor**.                                      |
| **Phylogenetic tree places a MAG far from its taxonomic label**                                             | MAG mis-binned or chimeric; taxonomy DB lagging behind GTDB                                                                               | • Re-run bin QC (CheckM 2, GUNC).<br>• If ANI to close genomes < 95 %, register as *uSGB* (unknown species) and treat novelty as biological, not error.                                                                                 |
| **Copy-number ≠ protein abundance direction**                                                               | Regulation, post-translational turnover, operon structure                                                                                 | • Don’t force reconciliation—this is biology!<br>• Summarise to pathway level (e.g. KEGG module completeness × median peptide intensity) to integrate both layers. ([PMC][5])                                                           |

---

### 3  |  Statistical integration tricks that tolerate layer discordance

| Approach                                  | One-liner recipe                                                                                                                        |                                                                                                     |
| ----------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| **Hurdle models for proteins**            | \`glmmTMB(abundance \~ group + (1                                                                                                       | subject), ziformula = \~1, family = truncated\_gaussian)\` — treats non-detected as censored zeros. |
| **Multi-omics factor analysis (MOFA v2)** | Learns latent factors that capture covariation; discordant features load weakly and can be down-weighted without being dropped.         |                                                                                                     |
| **Graph-based reconciliation**            | Build a bipartite graph MAG↔protein, weight edges by peptide evidence; run *Label Propagation* to smooth inconsistent assignments.      |                                                                                                     |
| **Pathway-level roll-up**                 | Collapse gene and protein tables to KEGG modules / CAZy families; use module scores as features in paired tests (Wilcoxon signed-rank). |                                                                                                     |
| **SHAP discrepancy flag**                 | Train LightGBM on combined table; any feature with high SHAP *and* high layer-discordance gets a manual QC flag for wet-lab validation. |                                                                                                     |

---

### 4  |  Minimal QC checklist before you lock the analysis

1. **Re-validate MAG bin purity** (completeness ≥ 90 %, contamination ≤ 5 %).
2. **Cross-check ORF coordinates**: no premature stops or frame-shifts in top biomarkers.
3. **Positive-control peptides**: spike‐in or house-keeping proteins detected in ≥ 95 % of samples.
4. **Phylogeny sanity plot**: tree tip labels coloured by layer—look for stray placements.
5. **Permutation test of pairing**: shuffle subject IDs; AUROC should crash to \~0.5, proving no leakage.

---

### 5  |  Take-home

> **Perfect concordance is the exception, not the rule.**
> Treat layer mismatches first as a QC signal (possible artefact), then as *biology* (regulation, novelty) once artefacts are ruled out. Harmonise databases, use phylogeny-aware tools, and choose statistical models that don’t assume every gene must yield a peptide. That way your paired metagenome-metaproteome study remains both **robust** and **biologically insightful** in the face of inevitable inconsistencies.

[1]: https://www.nature.com/articles/s41587-023-01688-w?utm_source=chatgpt.com "Extending and improving metagenomic taxonomic profiling with ..."
[2]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11968815/?utm_source=chatgpt.com "Lineage-specific microbial protein prediction enables large-scale ..."
[3]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11131966/?utm_source=chatgpt.com "Pairing metagenomics and metaproteomics to characterize ..."
[4]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9040630/?utm_source=chatgpt.com "Phylogeny-Aware Analysis of Metagenome Community Ecology ..."
[5]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6372841/?utm_source=chatgpt.com "Evaluating Metagenomic Prediction of the Metaproteome in a 4.5 ..."



