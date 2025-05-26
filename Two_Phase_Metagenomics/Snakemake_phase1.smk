#!/usr/bin/env python3

"""
Metagenomics Shotgun Pipeline - Single Sample Processing
Focus: Community Structure + Novelty Detection (Assembly-First Approach)
"""

import os
import pandas as pd
from snakemake.utils import min_version

min_version("7.0")

# Configuration
configfile: "config/config.yaml"

# Sample information
SAMPLES = config["samples"]
READS = ["R1", "R2"]

# Reference databases paths
GTDB_DB = config["databases"]["gtdb"]
HOST_GENOME = config["databases"]["host_genome"]
NCBI_DB = config["databases"]["ncbi"]
CHECKM2_DB = config["databases"]["checkm2"]

# Parameters
THREADS = config["threads"]
MEMORY = config["memory"]

# Output directories
RESULTS_DIR = config["output_dir"]
QC_DIR = f"{RESULTS_DIR}/01_quality_control"
ASSEMBLY_DIR = f"{RESULTS_DIR}/02_assembly"
BINNING_DIR = f"{RESULTS_DIR}/03_binning"
TAXONOMY_DIR = f"{RESULTS_DIR}/04_taxonomy"
NOVELTY_DIR = f"{RESULTS_DIR}/05_novelty"
ABUNDANCE_DIR = f"{RESULTS_DIR}/06_abundance"

# Target rule
rule all:
    input:
        # Quality control outputs
        expand(f"{QC_DIR}/{{sample}}/{{sample}}_{{read}}_fastqc.html", sample=SAMPLES, read=READS),
        expand(f"{QC_DIR}/{{sample}}/{{sample}}_trimmed_{{read}}.fastq.gz", sample=SAMPLES, read=READS),
        expand(f"{QC_DIR}/{{sample}}/{{sample}}_decontaminated_{{read}}.fastq.gz", sample=SAMPLES, read=READS),
        
        # Assembly outputs
        expand(f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta", sample=SAMPLES),
        expand(f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_assembly_stats.txt", sample=SAMPLES),
        
        # Binning outputs
        expand(f"{BINNING_DIR}/{{sample}}/refined_bins/", sample=SAMPLES),
        expand(f"{BINNING_DIR}/{{sample}}/checkm2_results/quality_report.tsv", sample=SAMPLES),
        
        # Taxonomy outputs
        expand(f"{TAXONOMY_DIR}/{{sample}}/gtdbtk/gtdbtk.bac120.summary.tsv", sample=SAMPLES),
        expand(f"{TAXONOMY_DIR}/{{sample}}/cat_bat/{{sample}}_contigs.classification", sample=SAMPLES),
        
        # Novelty detection outputs
        expand(f"{NOVELTY_DIR}/{{sample}}/homology_novelty.tsv", sample=SAMPLES),
        expand(f"{NOVELTY_DIR}/{{sample}}/ml_novelty_scores.tsv", sample=SAMPLES),
        expand(f"{NOVELTY_DIR}/{{sample}}/combined_novelty_report.tsv", sample=SAMPLES),
        
        # Abundance outputs
        expand(f"{ABUNDANCE_DIR}/{{sample}}/contig_abundance.tsv", sample=SAMPLES),
        expand(f"{ABUNDANCE_DIR}/{{sample}}/mag_abundance.tsv", sample=SAMPLES)

# ============================================================================
# QUALITY CONTROL & PREPROCESSING
# ============================================================================

rule fastqc_raw:
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["R1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["R2"]
    output:
        r1_html = f"{QC_DIR}/{{sample}}/{{sample}}_R1_fastqc.html",
        r2_html = f"{QC_DIR}/{{sample}}/{{sample}}_R2_fastqc.html",
        r1_zip = f"{QC_DIR}/{{sample}}/{{sample}}_R1_fastqc.zip",
        r2_zip = f"{QC_DIR}/{{sample}}/{{sample}}_R2_fastqc.zip"
    params:
        outdir = f"{QC_DIR}/{{sample}}"
    threads: 4
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {threads} {input.r1} {input.r2}
        """

rule trim_adapters:
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["R1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["R2"]
    output:
        r1 = f"{QC_DIR}/{{sample}}/{{sample}}_trimmed_R1.fastq.gz",
        r2 = f"{QC_DIR}/{{sample}}/{{sample}}_trimmed_R2.fastq.gz",
        report = f"{QC_DIR}/{{sample}}/{{sample}}_fastp_report.html"
    params:
        min_length = config["trimming"]["min_length"],
        quality_threshold = config["trimming"]["quality_threshold"]
    threads: 8
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -h {output.report} \
              -q {params.quality_threshold} \
              -l {params.min_length} \
              --thread {threads} \
              --detect_adapter_for_pe
        """

rule host_decontamination:
    input:
        r1 = f"{QC_DIR}/{{sample}}/{{sample}}_trimmed_R1.fastq.gz",
        r2 = f"{QC_DIR}/{{sample}}/{{sample}}_trimmed_R2.fastq.gz",
        host_index = f"{HOST_GENOME}.1.bt2"
    output:
        r1 = f"{QC_DIR}/{{sample}}/{{sample}}_decontaminated_R1.fastq.gz",
        r2 = f"{QC_DIR}/{{sample}}/{{sample}}_decontaminated_R2.fastq.gz",
        log = f"{QC_DIR}/{{sample}}/{{sample}}_decontamination.log"
    params:
        host_prefix = HOST_GENOME
    threads: 16
    shell:
        """
        bowtie2 -x {params.host_prefix} \
                -1 {input.r1} -2 {input.r2} \
                --un-conc-gz {QC_DIR}/{wildcards.sample}/{wildcards.sample}_decontaminated_R%.fastq.gz \
                -p {threads} \
                --very-sensitive \
                -S /dev/null 2> {output.log}
        """

# ============================================================================
# ASSEMBLY & BINNING
# ============================================================================

rule metaspades_assembly:
    input:
        r1 = f"{QC_DIR}/{{sample}}/{{sample}}_decontaminated_R1.fastq.gz",
        r2 = f"{QC_DIR}/{{sample}}/{{sample}}_decontaminated_R2.fastq.gz"
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta",
        scaffolds = f"{ASSEMBLY_DIR}/{{sample}}/scaffolds.fasta"
    params:
        outdir = f"{ASSEMBLY_DIR}/{{sample}}",
        memory = MEMORY
    threads: 32
    shell:
        """
        metaspades.py -1 {input.r1} -2 {input.r2} \
                      -o {params.outdir} \
                      -t {threads} \
                      -m {params.memory} \
                      --meta
        
        # Rename and filter contigs (>1kb)
        seqtk seq -L 1000 {params.outdir}/contigs.fasta > {output.contigs}
        cp {params.outdir}/scaffolds.fasta {output.scaffolds}
        """

rule assembly_stats:
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta"
    output:
        stats = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_assembly_stats.txt"
    script:
        "scripts/assembly_stats.py"

rule map_reads_to_contigs:
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta",
        r1 = f"{QC_DIR}/{{sample}}/{{sample}}_decontaminated_R1.fastq.gz",
        r2 = f"{QC_DIR}/{{sample}}/{{sample}}_decontaminated_R2.fastq.gz"
    output:
        bam = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_mapped.bam",
        bai = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_mapped.bam.bai"
    threads: 16
    shell:
        """
        # Index contigs
        bwa index {input.contigs}
        
        # Map reads
        bwa mem -t {threads} {input.contigs} {input.r1} {input.r2} | \
        samtools sort -@ {threads} -o {output.bam}
        
        # Index BAM
        samtools index {output.bam}
        """

rule calculate_contig_depth:
    input:
        bam = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_mapped.bam"
    output:
        depth = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contig_depth.txt"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        """

rule metabat2_binning:
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta",
        depth = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contig_depth.txt"
    output:
        directory(f"{BINNING_DIR}/{{sample}}/metabat2")
    params:
        prefix = f"{BINNING_DIR}/{{sample}}/metabat2/bin"
    threads: 16
    shell:
        """
        mkdir -p {output}
        metabat2 -i {input.contigs} \
                 -a {input.depth} \
                 -o {params.prefix} \
                 -t {threads} \
                 -m 1500
        """

rule concoct_binning:
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta",
        bam = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_mapped.bam"
    output:
        directory(f"{BINNING_DIR}/{{sample}}/concoct")
    threads: 16
    shell:
        """
        mkdir -p {output}
        
        # Cut contigs into chunks
        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b {output}/contigs_10K.bed > {output}/contigs_10K.fa
        
        # Generate coverage table
        concoct_coverage_table.py {output}/contigs_10K.bed {input.bam} > {output}/coverage_table.tsv
        
        # Run CONCOCT
        concoct --composition_file {output}/contigs_10K.fa \
                --coverage_file {output}/coverage_table.tsv \
                -b {output}/ \
                -t {threads}
        
        # Merge clustering
        merge_cutup_clustering.py {output}/clustering_gt1000.csv > {output}/clustering_merged.csv
        
        # Extract bins
        extract_fasta_bins.py {input.contigs} {output}/clustering_merged.csv --output_path {output}/bins/
        """

rule maxbin2_binning:
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta",
        r1 = f"{QC_DIR}/{{sample}}/{{sample}}_decontaminated_R1.fastq.gz",
        r2 = f"{QC_DIR}/{{sample}}/{{sample}}_decontaminated_R2.fastq.gz"
    output:
        directory(f"{BINNING_DIR}/{{sample}}/maxbin2")
    threads: 16
    shell:
        """
        mkdir -p {output}
        
        # Create abundance file
        cat {input.r1} {input.r2} > {output}/reads.fastq.gz
        
        run_MaxBin.pl -contig {input.contigs} \
                      -reads {output}/reads.fastq.gz \
                      -out {output}/bin \
                      -thread {threads}
        """

rule das_tool_refinement:
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta",
        metabat2 = f"{BINNING_DIR}/{{sample}}/metabat2",
        concoct = f"{BINNING_DIR}/{{sample}}/concoct",
        maxbin2 = f"{BINNING_DIR}/{{sample}}/maxbin2"
    output:
        directory(f"{BINNING_DIR}/{{sample}}/refined_bins")
    params:
        prefix = f"{BINNING_DIR}/{{sample}}/das_tool"
    threads: 16
    shell:
        """
        # Prepare binning results for DAS Tool
        Fasta_to_Scaffolds2Bin.sh -i {input.metabat2} -e fa > {params.prefix}_metabat2.tsv
        Fasta_to_Scaffolds2Bin.sh -i {input.concoct}/bins -e fa > {params.prefix}_concoct.tsv
        Fasta_to_Scaffolds2Bin.sh -i {input.maxbin2} -e fasta > {params.prefix}_maxbin2.tsv
        
        # Run DAS Tool
        DAS_Tool -i {params.prefix}_metabat2.tsv,{params.prefix}_concoct.tsv,{params.prefix}_maxbin2.tsv \
                 -l metabat2,concoct,maxbin2 \
                 -c {input.contigs} \
                 -o {params.prefix} \
                 --write_bins 1 \
                 -t {threads}
        
        # Move final bins
        mkdir -p {output}
        if [ -d {params.prefix}_DASTool_bins ]; then
            cp {params.prefix}_DASTool_bins/* {output}/
        fi
        """

rule checkm2_quality:
    input:
        bins = f"{BINNING_DIR}/{{sample}}/refined_bins"
    output:
        report = f"{BINNING_DIR}/{{sample}}/checkm2_results/quality_report.tsv"
    params:
        outdir = f"{BINNING_DIR}/{{sample}}/checkm2_results",
        db = CHECKM2_DB
    threads: 16
    shell:
        """
        checkm2 predict --threads {threads} \
                       --input {input.bins} \
                       --output-directory {params.outdir} \
                       --database_path {params.db}
        """

# ============================================================================
# TAXONOMIC PROFILING
# ============================================================================

rule gtdbtk_classify:
    input:
        bins = f"{BINNING_DIR}/{{sample}}/refined_bins"
    output:
        summary = f"{TAXONOMY_DIR}/{{sample}}/gtdbtk/gtdbtk.bac120.summary.tsv"
    params:
        outdir = f"{TAXONOMY_DIR}/{{sample}}/gtdbtk",
        db = GTDB_DB
    threads: 32
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}
        
        gtdbtk classify_wf --genome_dir {input.bins} \
                          --out_dir {params.outdir} \
                          --cpus {threads} \
                          --extension fa
        """

rule cat_bat_classification:
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta"
    output:
        classification = f"{TAXONOMY_DIR}/{{sample}}/cat_bat/{{sample}}_contigs.classification"
    params:
        outdir = f"{TAXONOMY_DIR}/{{sample}}/cat_bat",
        db = config["databases"]["cat_bat"]
    threads: 16
    shell:
        """
        mkdir -p {params.outdir}
        
        CAT contigs -c {input.contigs} \
                   -d {params.db}/CAT_database \
                   -t {params.db}/CAT_taxonomy \
                   -o {params.outdir}/{wildcards.sample}_contigs \
                   --force \
                   -n {threads}
        """

# ============================================================================
# NOVELTY DETECTION
# ============================================================================

rule homology_based_novelty:
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta"
    output:
        blast_results = f"{NOVELTY_DIR}/{{sample}}/blast_results.tsv",
        novelty_report = f"{NOVELTY_DIR}/{{sample}}/homology_novelty.tsv"
    params:
        db = NCBI_DB
    threads: 16
    shell:
        """
        mkdir -p {NOVELTY_DIR}/{wildcards.sample}
        
        # DIAMOND BLAST against NCBI
        diamond blastx -q {input.contigs} \
                      -d {params.db} \
                      -o {output.blast_results} \
                      -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
                      --threads {threads} \
                      --max-target-seqs 1 \
                      --evalue 1e-5
        
        # Analyze novelty based on homology
        python scripts/analyze_homology_novelty.py \
               --contigs {input.contigs} \
               --blast {output.blast_results} \
               --output {output.novelty_report}
        """

rule ml_novelty_detection:
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta"
    output:
        embeddings = f"{NOVELTY_DIR}/{{sample}}/dnabert_embeddings.npy",
        novelty_scores = f"{NOVELTY_DIR}/{{sample}}/ml_novelty_scores.tsv"
    threads: 8
    shell:
        """
        # Generate DNABert embeddings and run Isolation Forest
        python scripts/ml_novelty_detection.py \
               --contigs {input.contigs} \
               --embeddings {output.embeddings} \
               --scores {output.novelty_scores} \
               --threads {threads}
        """

rule combine_novelty_results:
    input:
        homology = f"{NOVELTY_DIR}/{{sample}}/homology_novelty.tsv",
        ml_scores = f"{NOVELTY_DIR}/{{sample}}/ml_novelty_scores.tsv"
    output:
        combined = f"{NOVELTY_DIR}/{{sample}}/combined_novelty_report.tsv"
    script:
        "scripts/combine_novelty_results.py"

# ============================================================================
# ABUNDANCE ESTIMATION
# ============================================================================

rule calculate_contig_abundance:
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_contigs.fasta",
        bam = f"{ASSEMBLY_DIR}/{{sample}}/{{sample}}_mapped.bam"
    output:
        abundance = f"{ABUNDANCE_DIR}/{{sample}}/contig_abundance.tsv"
    script:
        "scripts/calculate_contig_abundance.py"

rule calculate_mag_abundance:
    input:
        bins = f"{BINNING_DIR}/{{sample}}/refined_bins",
        contig_abundance = f"{ABUNDANCE_DIR}/{{sample}}/contig_abundance.tsv"
    output:
        mag_abundance = f"{ABUNDANCE_DIR}/{{sample}}/mag_abundance.tsv"
    script:
        "scripts/calculate_mag_abundance.py"

