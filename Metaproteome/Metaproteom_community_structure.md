Below is a **community-centric roadmap** that layers diversity analysis on top of the metagenome + metaproteome novelty pipeline we sketched earlier.  I split it into:

1. **Hypothesis sets you can actually test** (case-control *vs.* longitudinal)
2. **Step-by-step workflow** pointing to recent methods, with notes on where diversity, differential-abundance, and novelty intersect.

---

## 1  |  Hypotheses that put “community change” front-and-centre

| Study design     | Example, ready-to-test hypotheses                                                                                                                                                             | Key response variables                                                                                              |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **Case–control** | *Inflammatory-bowel-disease (IBD) patients have lower phylogenetic α-diversity but higher “novel protein” richness than matched healthy controls.*                                            | Shannon / Faith PD (16S or MAGs); richness of proteins with BLAST E-value > 1e-5; MaAsLin3-derived log-fold changes |
|                  | *The relative abundance of carbohydrate-active MAGs (and their expressed CAZyme peptides) is enriched in high-fiber diet responders.*                                                         | Differential MAG abundance (ANCOM-BC2); spectral-LFQ CAZyme intensities                                             |
| **Time-series**  | *After short-course ciprofloxacin, β-diversity (Aitchison distance) rebounds by week 6, but the fraction of outlier proteins (ESM-IsolationForest score > 0.9) keeps rising through week 12.* | Aitchison distances; novelty index per sample; MetaLonDA-identified time windows                                    |
|                  | *Novel fungal peptides appear **prior** to a bloom of fungal MAGs, suggesting transcriptional response precedes detectable genome-level shifts.*                                              | Cross-correlation of fungal peptide counts vs MAG relative abundance                                                |

Anchor on one primary hypothesis and phrase the rest as secondary so your statistical corrections aren’t brutal.

---

## 2  |  Analysis roadmap with diversity / differential-abundance hooks

### A. Study design & metadata

* ≥ 3–5 biological replicates per group or per time-point.
* Record diet, antibiotics, batch, instrument run-order (covariates for MaAsLin3 models).

---

### B. Pre-processing (same as before)

* **Metagenomes:** QC → assembly → *VAMB* binning → *GTDB-Tk* taxonomy / ANI.
* **Metaproteomes:** RAW → mzML → *MSFragger+Philosopher* search against the custom MAG-ORF FASTA.
* **Integration:** use *M2P* or the Metagenome-Metaproteome (M2P) pipeline to map proteins back to MAGs and produce feature tables (taxa, genes, KO functions, protein groups). ([Nature][1])

---

### C. Feature tables you’ll feed into diversity stats

| Table                                                   | Transformation                    | Diversity statistics               |
| ------------------------------------------------------- | --------------------------------- | ---------------------------------- |
| MAG relative abundance (reads-per-kilo)                 | Centered log-ratio (CLR)          | α: Shannon, Faith PD; β: Aitchison |
| Protein LFQ intensities (grouped by MAG or KO)          | CLR or Variance-stabilizing (VSN) | same α/β; functional Bray-Curtis   |
| Novel protein indicator (0 = reference-like, 1 = novel) | Raw counts per sample             | Richness, binomial mixed models    |

---

### D. Differential-abundance & association tests

| Scenario                           | Recommended method                                                                                                                                                        | Why                                                               |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| **Cross-sectional (case–control)** | **ANCOM-BC 2** for taxa and KO counts ([PMC][2]); **MaAsLin 3** for multivariable models with covariates ([GitHub][3], [huttenhower.sph.harvard.edu][4])                  | Handles compositionality, zero inflation; supports random effects |
| **Time-series**                    | **MetaLonDA** for per-feature time-interval DE ([GitHub][5], [ResearchGate][6]); **GLMM/GAMM** (e.g., lme4 or mgcv) or Bayesian compositional GLMMs ([BioMed Central][7]) | Captures subject-specific trajectories & non-linear trends        |
| **Proteomic intensities**          | **limma-voom** (case-control) or **dream mixed model** (time-series) after VSN                                                                                            | Well-validated for mass-spec counts                               |

> **Tip:** export the same sample-level covariate matrix for all tables so MaAsLin3 and limma share the exact fixed/random effects.

---

### E. Diversity & novelty trend tests

1. **Alpha diversity:** compare by Wilcoxon (case-control) or mixed-effects ANOVA (time-series).
2. **Beta diversity:** PERMANOVA + PERMDISP to confirm dispersion effects.
3. **Novelty index:**

   * Per-sample proportion of *novel proteins* = `sum(novel_PSMs) / total_PSMs`.
   * Fit a mixed-effects model: `novel_prop ~ group*time + (1|subject)`.
   * Overlay on β-diversity PCoA ordination to see if “drift” is driven by novelty.

---

### F. Multi-omics coupling (taxon ↔ protein)

| Goal                                            | Tool                                                                 |
| ----------------------------------------------- | -------------------------------------------------------------------- |
| Link MAG abundance shifts to expressed proteins | **DIABLO** block-PLS in **mixOmics**                                 |
| Co-trajectory ordination                        | **MOFA 2** latent factors (handles missing views)                    |
| Visualise time-lagged relationships             | Cross-correlation heat-map; **dynOmics** for time-warped multi-omics |

---

### G. Visual reporting staples

* **Stacked-area time-series** of top 20 taxa + overlaid line for novelty index.
* **PCoA / UMAP** coloured by time or case status with arrows showing subject trajectories.
* **Volcano / effect-size plot** for differential proteins (case-control) or **MetaLonDA “heat-ribbon”** showing significant time intervals.
* **Sankey**: MAG → KO → protein group, highlighting those with novelty-score > 0.9.

---

## Pitfalls & pro-tips

* **Batch in time-series**: randomise extraction & MS run-order, include iRT peptides.
* **Compositional zeros**: add a small-count pseudo before CLR or use Bayesian multiplicative replacement.
* **Spectral depth**: supplement DDA with DIA or BoxCar to avoid missing low-abundance proteins.
* **False discovery inflation**: MaAsLin3/ANCOM-BC2 already do per-feature multiple-testing; still run Benjamini–Hochberg on any ad-hoc tests.
* **Interpreting novelty**: filter truncated ORFs (< 60 aa) and coverage outliers so “novel” ≠ artefact.

---

### Recent exemplars to crib from

* **Population-based metaproteomics reveals functional associations with microbial diversity** (bioRxiv 2024) – combines metaproteome LFQ with metagenomic MAG profiles, finds metaproteome captures shifts earlier than metagenome. ([bioRxiv][8])
* **Longitudinal serum proteome mapping** (Nat Metab 2024) – nice template for mixed-model treatment of protein trajectories. ([Nature][1])

---

### Key schematic

> **Community change = Taxonomic shifts (metagenome) + Functional shifts (metaproteome) + Novelty trajectories.**
> Use compositional mixed-models (ANCOM-BC2, MetaLonDA, MaAsLin3) to nail *where/when* shifts occur, then overlay **novelty scores** to see *how much* of that change is driven by previously unseen genes or proteins.

Follow that recipe and you’ll report not just “the microbiome changed,” but **which taxa/functions drive the change, when they emerge, and whether they represent genuinely new biology**—the story reviewers love.

[1]: https://www.nature.com/articles/s42255-024-01185-7?utm_source=chatgpt.com "Longitudinal serum proteome mapping reveals biomarkers for ..."
[2]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10187376/?utm_source=chatgpt.com "Multi-group Analysis of Compositions of Microbiomes with Covariate ..."
[3]: https://github.com/biobakery/maaslin3?utm_source=chatgpt.com "MaAsLin3: Microbiome Multivariate Association with Linear Models"
[4]: https://huttenhower.sph.harvard.edu/maaslin3/?utm_source=chatgpt.com "MaAsLin 3 - The Huttenhower Lab"
[5]: https://github.com/aametwally/MetaLonDA?utm_source=chatgpt.com "aametwally/MetaLonDA: METAgenomic LONgitudinal Differential ..."
[6]: https://www.researchgate.net/publication/323157312_MetaLonDA_A_flexible_R_package_for_identifying_time_intervals_of_differentially_abundant_features_in_metagenomic_longitudinal_studies?utm_source=chatgpt.com "(PDF) MetaLonDA: A flexible R package for identifying time intervals ..."
[7]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-025-06114-3?utm_source=chatgpt.com "Bayesian compositional generalized linear mixed models for ..."
[8]: https://www.biorxiv.org/content/10.1101/2024.11.26.625569v2.full-text?utm_source=chatgpt.com "Population-based metaproteomics reveals functional associations ..."



