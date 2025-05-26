Phase 1: Single Sample Processing (Revised)

Key Innovations

Assembly-First Strategy: Better captures novel organisms than read-based approaches
DNABert + Isolation Forest: Modern ML approach to detect novel sequence patterns
Integrated Novelty Scoring: Combines multiple detection methods for robust results
Contig-Based Analysis: All diversity calculations use assembled sequences
Comprehensive Validation: Automated testing ensures pipeline reliability

 Expected Outputs

Assembly statistics and quality metrics
Novel sequence candidates with confidence scores and method agreement
MAG-level taxonomy and abundance profiles
Interactive visualizations for result exploration
Comprehensive reports in HTML and markdown formats

1.1 Quality Control & Preprocessing

FastQC → quality assessment
Trimmomatic/Fastp → adapter removal, quality trimming
Host decontamination → BWA/Bowtie2 against host genome

1.2 Assembly & Binning

MetaSPAdes/MEGAHIT → de novo assembly
Assembly quality assessment → MetaQUAST, assembly statistics
Read mapping back → BWA/Bowtie2 to calculate contig coverage
MetaBAT2/CONCOCT/MaxBin2 → binning for MAGs
Bin refinement → DAS Tool for consensus binning
CheckM2/GUNC → MAG quality and contamination assessment

1.3 Taxonomic Profiling (Assembly-Based)

GTDB-Tk → taxonomic classification of MAGs
CAT/BAT → contig/bin taxonomic classification
Contig-level classification → DIAMOND/MMseqs2 against GTDB/NCBI
Phylogenetic placement → for novel lineages

1.4 Abundance Estimation

Coverage calculation → from read mapping to contigs/MAGs
MAG abundance → coverage-based relative abundance
Read recruitment → percentage of reads assigned to each MAG

1.5 Enhanced Novelty Detection (Single Sample)
Assembly → Contig Quality Filter → Multi-layered Novelty Detection:

A. Homology-Based:
   - DIAMOND vs NCBI/GTDB
   - Low-identity/coverage contigs

B. ML-Based (Your Addition):
   - DNABert embedding generation
   - Isolation Forest outlier detection
   - Anomaly score ranking

C. Compositional:
   - GC content deviation
   - Codon usage bias
   - k-mer signature analysis

1.6 Contig Classification & Abundance

Novel contig prioritization → combine homology + ML scores
Taxonomic placement → phylogenetic analysis for high-scoring outliers
Abundance calculation → coverage-based, including novel contigs

Recommended Hybrid Approach:

Primary pipeline: Assembly-first for comprehensive novelty detection
Validation track: Read-based profiling for comparison and validation
Integration: Combine insights from both approaches
Quality control: Ensure assembly-based results are robust


