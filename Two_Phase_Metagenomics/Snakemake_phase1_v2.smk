# =============================================================================
# 🧬 Metagenomics Pipeline - Phase 1: Single-Sample Processing
# Snakefile for comprehensive metagenomic analysis with ML-enhanced novelty detection
# =============================================================================

import os
from pathlib import Path
import pandas as pd

# Load configuration
configfile: "config/config.yaml"

# Get sample information
SAMPLES = list(config["samples"].keys())

# Define output directories
RESULTS_DIR = config["directories"]["results"]
TEMP_DIR = config["directories"]["temp"]
LOGS_DIR = config["directories"]["logs"]

# Create directories
for directory in [RESULTS_DIR, TEMP_DIR, LOGS_DIR]:
    Path(directory).mkdir(parents=True, exist_ok=True)

# =============================================================================
# Target Rules
# =============================================================================

rule all:
    """Main target rule - runs complete Phase 1 pipeline"""
    input:
        # Quality control
        expand("{results}/01_quality_control/{sample}_fastqc_report.html", 
               results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/01_quality_control/{sample}_trimmed_R1.fastq.gz", 
               results=RESULTS_DIR, sample=SAMPLES),
        
        # Assembly
        expand("{results}/02_assembly/{sample}_contigs.fasta", 
               results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/02_assembly/{sample}_assembly_stats.json", 
               results=RESULTS_DIR, sample=SAMPLES),
        
        # Binning
        expand("{results}/03_binning/{sample}/final_bins", 
               results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/03_binning/{sample}_bin_quality.tsv", 
               results=RESULTS_DIR, sample=SAMPLES),
        
        # Taxonomy
        expand("{results}/04_taxonomy/{sample}_gtdbtk_summary.tsv", 
               results=RESULTS_DIR, sample=SAMPLES),
        
        # Novelty detection
        expand("{results}/05_novelty/{sample}_combined_novelty.tsv", 
               results=RESULTS_DIR, sample=SAMPLES),
        
        # Abundance calculation
        expand("{results}/06_abundance/{sample}_contig_abundance.tsv", 
               results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/06_abundance/{sample}_mag_abundance.tsv", 
               results=RESULTS_DIR, sample=SAMPLES),
        
        # Final report
        expand("{results}/07_reports/{sample}_phase1_summary.html", 
               results=RESULTS_DIR, sample=SAMPLES)

rule quality_only:
    """Run only quality control and assembly"""
    input:
        expand("{results}/02_assembly/{sample}_assembly_stats.json", 
               results=RESULTS_DIR, sample=SAMPLES)

rule novelty_only:
    """Run only novelty detection (requires assembly)"""
    input:
        expand("{results}/05_novelty/{sample}_combined_novelty.tsv", 
               results=RESULTS_DIR, sample=SAMPLES)

rule abundance_only:
    """Run only abundance calculation (requires assembly)"""
    input:
        expand("{results}/06_abundance/{sample}_mag_abundance.tsv", 
               results=RESULTS_DIR, sample=SAMPLES)

# =============================================================================
# Quality Control Rules
# =============================================================================

rule fastqc_raw:
    """Run FastQC on raw reads"""
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample] + "_R1.fastq.gz",
        r2 = lambda wildcards: config["samples"][wildcards.sample] + "_R2.fastq.gz"
    output:
        html = "{results}/01_quality_control/{sample}_fastqc_report.html",
        zip = "{results}/01_quality_control/{sample}_fastqc_data.zip"
    params:
        outdir = "{results}/01_quality_control",
        threads = config["quality_control"]["fastqc"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        runtime = 60
    log:
        "{logs}/fastqc/{sample}.log"
    conda:
        "environment.yml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc {input.r1} {input.r2} \
            --outdir {params.outdir} \
            --threads {params.threads} \
            --extract \
            > {log} 2>&1
        
        # Combine reports
        cat {params.outdir}/*_fastqc.html > {output.html}
        """

rule trim_reads:
    """Trim and filter reads with fastp"""
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample] + "_R1.fastq.gz",
        r2 = lambda wildcards: config["samples"][wildcards.sample] + "_R2.fastq.gz"
    output:
        r1 = "{results}/01_quality_control/{sample}_trimmed_R1.fastq.gz",
        r2 = "{results}/01_quality_control/{sample}_trimmed_R2.fastq.gz",
        html = "{results}/01_quality_control/{sample}_fastp_report.html",
        json = "{results}/01_quality_control/{sample}_fastp_report.json"
    params:
        qualified_quality_phred = config["quality_control"]["fastp"]["qualified_quality_phred"],
        unqualified_percent_limit = config["quality_control"]["fastp"]["unqualified_percent_limit"],
        length_required = config["quality_control"]["fastp"]["length_required"],
        threads = config["quality_control"]["fastp"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        runtime = 120
    log:
        "{logs}/fastp/{sample}.log"
    conda:
        "environment.yml"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --unqualified_percent_limit {params.unqualified_percent_limit} \
            --length_required {params.length_required} \
            --detect_adapter_for_pe \
            --correction \
            --html {output.html} \
            --json {output.json} \
            --thread {params.threads} \
            > {log} 2>&1
        """

# =============================================================================
# Assembly Rules
# =============================================================================

rule metaspades_assembly:
    """Assemble reads with MetaSPAdes"""
    input:
        r1 = "{results}/01_quality_control/{sample}_trimmed_R1.fastq.gz",
        r2 = "{results}/01_quality_control/{sample}_trimmed_R2.fastq.gz"
    output:
        contigs = "{results}/02_assembly/{sample}_contigs.fasta",
        scaffolds = "{results}/02_assembly/{sample}_scaffolds.fasta",
        assembly_graph = "{results}/02_assembly/{sample}_assembly_graph.gfa"
    params:
        outdir = "{results}/02_assembly/{sample}_metaspades",
        k_list = config["assembly"]["metaspades"]["k_list"],
        threads = config["assembly"]["metaspades"]["threads"],
        memory = config["assembly"]["metaspades"]["memory"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 128000,
        runtime = 1440  # 24 hours
    log:
        "{logs}/assembly/{sample}_metaspades.log"
    conda:
        "environment.yml"
    shell:
        """
        metaspades.py \
            -1 {input.r1} -2 {input.r2} \
            -o {params.outdir} \
            -k {params.k_list} \
            --meta \
            --threads {params.threads} \
            --memory {params.memory} \
            > {log} 2>&1
        
        # Copy main outputs
        cp {params.outdir}/contigs.fasta {output.contigs}
        cp {params.outdir}/scaffolds.fasta {output.scaffolds}
        cp {params.outdir}/assembly_graph.gfa {output.assembly_graph}
        """

rule filter_contigs:
    """Filter contigs by length and quality"""
    input:
        contigs = "{results}/02_assembly/{sample}_contigs.fasta"
    output:
        filtered = "{results}/02_assembly/{sample}_contigs_filtered.fasta"
    params:
        min_length = config["assembly"]["quality_filters"]["min_contig_length"],
        min_coverage = config["assembly"]["quality_filters"]["min_coverage"]
    resources:
        mem_mb = 8000,
        runtime = 30
    log:
        "{logs}/assembly/{sample}_filter_contigs.log"
    conda:
        "environment.yml"
    shell:
        """
        python scripts/filter_contigs.py \
            --input {input.contigs} \
            --output {output.filtered} \
            --min-length {params.min_length} \
            --min-coverage {params.min_coverage} \
            > {log} 2>&1
        """

rule assembly_statistics:
    """Calculate assembly statistics"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta"
    output:
        stats = "{results}/02_assembly/{sample}_assembly_stats.json",
        summary = "{results}/02_assembly/{sample}_assembly_summary.tsv",
        contigs = "{results}/02_assembly/{sample}_contig_stats.tsv"
    params:
        outdir = "{results}/02_assembly",
        sample = "{sample}",
        min_length = config["assembly"]["quality_filters"]["min_contig_length"],
        plots = "--plots" if config.get("advanced", {}).get("detailed_logs", False) else ""
    resources:
        mem_mb = 16000,
        runtime = 60
    log:
        "{logs}/assembly/{sample}_stats.log"
    conda:
        "environment.yml"
    shell:
        """
        python scripts/assembly_stats.py \
            --assembly {input.assembly} \
            --output {params.outdir} \
            --sample-name {params.sample} \
            --min-length {params.min_length} \
            {params.plots} \
            > {log} 2>&1
        """

# =============================================================================
# Binning Rules
# =============================================================================

rule bwa_index:
    """Create BWA index for abundance calculation"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta"
    output:
        index = "{results}/02_assembly/{sample}_contigs_filtered.fasta.bwt"
    resources:
        mem_mb = 16000,
        runtime = 60
    log:
        "{logs}/binning/{sample}_bwa_index.log"
    conda:
        "environment.yml"
    shell:
        """
        bwa index {input.assembly} > {log} 2>&1
        """

rule map_reads_for_binning:
    """Map reads to assembly for binning depth calculation"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta",
        index = "{results}/02_assembly/{sample}_contigs_filtered.fasta.bwt",
        r1 = "{results}/01_quality_control/{sample}_trimmed_R1.fastq.gz",
        r2 = "{results}/01_quality_control/{sample}_trimmed_R2.fastq.gz"
    output:
        bam = "{results}/02_assembly/{sample}_mapped.bam",
        depth = "{results}/02_assembly/{sample}_depth.txt"
    params:
        threads = config["binning"]["metabat2"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        runtime = 240
    log:
        "{logs}/binning/{sample}_mapping.log"
    conda:
        "environment.yml"
    shell:
        """
        # Map reads
        bwa mem -t {params.threads} {input.assembly} {input.r1} {input.r2} | \
        samtools sort -@ {params.threads} -o {output.bam} -
        
        # Index BAM
        samtools index {output.bam}
        
        # Calculate depth
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {output.bam}
        
        > {log} 2>&1
        """

rule metabat2_binning:
    """Bin contigs with MetaBAT2"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta",
        depth = "{results}/02_assembly/{sample}_depth.txt"
    output:
        bins_dir = directory("{results}/03_binning/{sample}/metabat2"),
        bin_list = "{results}/03_binning/{sample}/metabat2_bins.txt"
    params:
        min_contig = config["binning"]["metabat2"]["min_contig_length"],
        threads = config["binning"]["metabat2"]["threads"]
    resources:
        mem_mb = 32000,
        runtime = 120
    log:
        "{logs}/binning/{sample}_metabat2.log"
    conda:
        "environment.yml"
    shell:
        """
        mkdir -p {output.bins_dir}
        
        metabat2 \
            -i {input.assembly} \
            -a {input.depth} \
            -o {output.bins_dir}/bin \
            -m {params.min_contig} \
            -t {params.threads} \
            > {log} 2>&1
        
        # List created bins
        ls {output.bins_dir}/*.fa > {output.bin_list} 2>/dev/null || touch {output.bin_list}
        """

rule maxbin2_binning:
    """Bin contigs with MaxBin2"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta",
        r1 = "{results}/01_quality_control/{sample}_trimmed_R1.fastq.gz",
        r2 = "{results}/01_quality_control/{sample}_trimmed_R2.fastq.gz"
    output:
        bins_dir = directory("{results}/03_binning/{sample}/maxbin2"),
        bin_list = "{results}/03_binning/{sample}/maxbin2_bins.txt"
    params:
        min_contig = config["binning"]["maxbin2"]["min_contig_length"],
        threads = config["binning"]["maxbin2"]["threads"]
    resources:
        mem_mb = 32000,
        runtime = 180
    log:
        "{logs}/binning/{sample}_maxbin2.log"
    conda:
        "environment.yml"
    shell:
        """
        mkdir -p {output.bins_dir}
        
        # Create reads list file
        echo "{input.r1}" > {output.bins_dir}/reads_list.txt
        echo "{input.r2}" >> {output.bins_dir}/reads_list.txt
        
        run_MaxBin.pl \
            -contig {input.assembly} \
            -reads_list {output.bins_dir}/reads_list.txt \
            -out {output.bins_dir}/bin \
            -min_contig_length {params.min_contig} \
            -thread {params.threads} \
            > {log} 2>&1
        
        # List created bins
        ls {output.bins_dir}/*.fasta > {output.bin_list} 2>/dev/null || touch {output.bin_list}
        """

rule concoct_binning:
    """Bin contigs with CONCOCT"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta",
        bam = "{results}/02_assembly/{sample}_mapped.bam"
    output:
        bins_dir = directory("{results}/03_binning/{sample}/concoct"),
        bin_list = "{results}/03_binning/{sample}/concoct_bins.txt"
    params:
        chunk_size = config["binning"]["concoct"]["chunk_size"],
        threads = config["binning"]["concoct"]["threads"]
    resources:
        mem_mb = 32000,
        runtime = 180
    log:
        "{logs}/binning/{sample}_concoct.log"
    conda:
        "environment.yml"
    shell:
        """
        mkdir -p {output.bins_dir}
        
        # Cut contigs into chunks
        cut_up_fasta.py {input.assembly} \
            -c {params.chunk_size} \
            -o 0 \
            --merge_last \
            -b {output.bins_dir}/contigs_10K.bed \
            > {output.bins_dir}/contigs_10K.fa
        
        # Generate coverage table
        concoct_coverage_table.py \
            {output.bins_dir}/contigs_10K.bed \
            {input.bam} \
            > {output.bins_dir}/coverage_table.tsv
        
        # Run CONCOCT
        concoct \
            --composition_file {output.bins_dir}/contigs_10K.fa \
            --coverage_file {output.bins_dir}/coverage_table.tsv \
            -b {output.bins_dir}/ \
            -t {params.threads}
        
        # Extract bins
        merge_cutup_clustering.py \
            {output.bins_dir}/clustering_gt1000.csv \
            > {output.bins_dir}/clustering_merged.csv
        
        extract_fasta_bins.py \
            {input.assembly} \
            {output.bins_dir}/clustering_merged.csv \
            --output_path {output.bins_dir}/
        
        # List created bins
        ls {output.bins_dir}/*.fa > {output.bin_list} 2>/dev/null || touch {output.bin_list}
        
        > {log} 2>&1
        """

rule dastool_refinement:
    """Refine bins with DAS Tool"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta",
        metabat2_bins = "{results}/03_binning/{sample}/metabat2_bins.txt",
        maxbin2_bins = "{results}/03_binning/{sample}/maxbin2_bins.txt",
        concoct_bins = "{results}/03_binning/{sample}/concoct_bins.txt"
    output:
        refined_bins = directory("{results}/03_binning/{sample}/final_bins"),
        summary = "{results}/03_binning/{sample}/dastool_summary.tsv"
    params:
        outdir = "{results}/03_binning/{sample}/dastool",
        threads = config["binning"]["dastool"]["threads"],
        score_threshold = config["binning"]["dastool"]["score_threshold"]
    resources:
        mem_mb = 32000,
        runtime = 120
    log:
        "{logs}/binning/{sample}_dastool.log"
    conda:
        "environment.yml"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Prepare input files for DAS Tool
        python scripts/prepare_dastool_input.py \
            --metabat2 {input.metabat2_bins} \
            --maxbin2 {input.maxbin2_bins} \
            --concoct {input.concoct_bins} \
            --output {params.outdir}
        
        # Run DAS Tool
        DAS_Tool \
            -i {params.outdir}/metabat2.scaffolds2bin.tsv,{params.outdir}/maxbin2.scaffolds2bin.tsv,{params.outdir}/concoct.scaffolds2bin.tsv \
            -l metabat2,maxbin2,concoct \
            -c {input.assembly} \
            -o {params.outdir}/DAS_Tool \
            --score_threshold {params.score_threshold} \
            --threads {params.threads} \
            --write_bins
        
        # Copy final bins
        mkdir -p {output.refined_bins}
        cp {params.outdir}/DAS_Tool_DASTool_bins/*.fa {output.refined_bins}/ 2>/dev/null || echo "No bins to copy"
        
        # Copy summary
        cp {params.outdir}/DAS_Tool_DASTool_summary.txt {output.summary} 2>/dev/null || touch {output.summary}
        
        > {log} 2>&1
        """

rule checkm_quality:
    """Assess bin quality with CheckM"""
    input:
        bins_dir = "{results}/03_binning/{sample}/final_bins"
    output:
        quality = "{results}/03_binning/{sample}_bin_quality.tsv",
        checkm_dir = directory("{results}/03_binning/{sample}/checkm")
    params:
        threads = config["binning"]["checkm"]["threads"],
        reduced_tree = "--reduced_tree" if config["binning"]["checkm"]["reduced_tree"] else ""
    resources:
        mem_mb = 64000,
        runtime = 240
    log:
        "{logs}/binning/{sample}_checkm.log"
    conda:
        "environment.yml"
    shell:
        """
        mkdir -p {output.checkm_dir}
        
        checkm lineage_wf \
            {input.bins_dir} \
            {output.checkm_dir} \
            --threads {params.threads} \
            {params.reduced_tree} \
            --tab_table \
            -f {output.quality} \
            > {log} 2>&1
        """

# =============================================================================
# Taxonomy Rules
# =============================================================================

rule gtdbtk_classify:
    """Classify MAGs with GTDB-Tk"""
    input:
        bins_dir = "{results}/03_binning/{sample}/final_bins"
    output:
        summary = "{results}/04_taxonomy/{sample}_gtdbtk_summary.tsv",
        ar122_summary = "{results}/04_taxonomy/{sample}/classify/gtdbtk.ar122.summary.tsv",
        bac120_summary = "{results}/04_taxonomy/{sample}/classify/gtdbtk.bac120.summary.tsv"
    params:
        outdir = "{results}/04_taxonomy/{sample}",
        database = config["taxonomy"]["gtdbtk"]["database_path"],
        threads = config["taxonomy"]["gtdbtk"]["threads"]
    resources:
        mem_mb = 64000,
        runtime = 480
    log:
        "{logs}/taxonomy/{sample}_gtdbtk.log"
    conda:
        "environment.yml"
    shell:
        """
        export GTDBTK_DATA_PATH={params.database}
        
        gtdbtk classify_wf \
            --genome_dir {input.bins_dir} \
            --out_dir {params.outdir} \
            --cpus {params.threads} \
            --extension fa
        
        # Combine summaries
        cat {output.ar122_summary} {output.bac120_summary} > {output.summary} 2>/dev/null || touch {output.summary}
        
        > {log} 2>&1
        """

rule kraken2_classify:
    """Additional classification with Kraken2"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta"
    output:
        classification = "{results}/04_taxonomy/{sample}_kraken2_classification.txt",
        report = "{results}/04_taxonomy/{sample}_kraken2_report.txt"
    params:
        database = config["taxonomy"]["kraken2"]["database"],
        threads = config["taxonomy"]["kraken2"]["threads"],
        confidence = config["taxonomy"]["kraken2"]["confidence"]
    resources:
        mem_mb = 32000,
        runtime = 60
    log:
        "{logs}/taxonomy/{sample}_kraken2.log"
    conda:
        "environment.yml"
    shell:
        """
        kraken2 \
            --db {params.database} \
            --threads {params.threads} \
            --confidence {params.confidence} \
            --output {output.classification} \
            --report {output.report} \
            {input.assembly} \
            > {log} 2>&1
        """

# =============================================================================
# Novelty Detection Rules
# =============================================================================

rule homology_novelty:
    """Detect novelty using homology-based methods"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta"
    output:
        results = "{results}/05_novelty/{sample}_homology_novelty.tsv",
        stats = "{results}/05_novelty/{sample}_homology_novelty_stats.json"
    params:
        outdir = "{results}/05_novelty",
        sample = "{sample}",
        databases = config["novelty_detection"]["homology"]["databases"],
        evalue = config["novelty_detection"]["homology"]["diamond"]["evalue"],
        threads = config["novelty_detection"]["homology"]["diamond"]["threads"],
        high_threshold = config["novelty_detection"]["homology"]["thresholds"]["high_novelty"],
        medium_threshold = config["novelty_detection"]["homology"]["thresholds"]["medium_novelty"],
        low_threshold = config["novelty_detection"]["homology"]["thresholds"]["low_novelty"],
        plots = "--plots" if config.get("advanced", {}).get("detailed_logs", False) else ""
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        runtime = 360
    log:
        "{logs}/novelty/{sample}_homology.log"
    conda:
        "environment.yml"
    shell:
        """
        python scripts/analyze_homology_novelty.py \
            --assembly {input.assembly} \
            --output {params.outdir} \
            --databases {params.databases} \
            --evalue {params.evalue} \
            --threads {params.threads} \
            --sample-name {params.sample} \
            --high-novelty-threshold {params.high_threshold} \
            --medium-novelty-threshold {params.medium_threshold} \
            --low-novelty-threshold {params.low_threshold} \
            {params.plots} \
            > {log} 2>&1
        """

rule ml_novelty:
    """Detect novelty using ML-based methods"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta"
    output:
        contigs = "{results}/05_novelty/{sample}_ml_novelty_contigs.tsv",
        detailed = "{results}/05_novelty/{sample}_ml_novelty_detailed.tsv",
        stats = "{results}/05_novelty/{sample}_ml_novelty_stats.json"
    params:
        outdir = "{results}/05_novelty",
        sample = "{sample}",
        model_name = config["novelty_detection"]["machine_learning"]["dnabert"]["model_name"],
        max_seq_length = config["novelty_detection"]["machine_learning"]["dnabert"]["max_sequence_length"],
        batch_size = config["novelty_detection"]["machine_learning"]["dnabert"]["batch_size"],
        device = config["novelty_detection"]["machine_learning"]["dnabert"]["device"],
        contamination = config["novelty_detection"]["machine_learning"]["isolation_forest"]["contamination"],
        n_estimators = config["novelty_detection"]["machine_learning"]["isolation_forest"]["n_estimators"],
        min_length = config["novelty_detection"]["machine_learning"]["preprocessing"]["min_length"],
        max_length = config["novelty_detection"]["machine_learning"]["preprocessing"]["max_length"],
        overlap_size = config["novelty_detection"]["machine_learning"]["preprocessing"]["overlap_size"],
        plots = "--plots" if config.get("advanced", {}).get("detailed_logs", False) else ""
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        runtime = 480,
        gpus = 1 if config["novelty_detection"]["machine_learning"]["dnabert"]["device"] == "cuda" else 0
    log:
        "{logs}/novelty/{sample}_ml.log"
    conda:
        "environment.yml"
    shell:
        """
        python scripts/ml_novelty_detection.py \
            --assembly {input.assembly} \
            --output {params.outdir} \
            --sample-name {params.sample} \
            --model-name {params.model_name} \
            --max-seq-length {params.max_seq_length} \
            --batch-size {params.batch_size} \
            --device {params.device} \
            --contamination {params.contamination} \
            --n-estimators {params.n_estimators} \
            --min-length {params.min_length} \
            --max-length {params.max_length} \
            --overlap-size {params.overlap_size} \
            {params.plots} \
            > {log} 2>&1
        """

rule combine_novelty:
    """Combine homology and ML novelty results"""
    input:
        homology = "{results}/05_novelty/{sample}_homology_novelty.tsv",
        ml = "{results}/05_novelty/{sample}_ml_novelty_contigs.tsv"
    output:
        combined = "{results}/05_novelty/{sample}_combined_novelty.tsv",
        novel = "{results}/05_novelty/{sample}_high_confidence_novel.tsv",
        stats = "{results}/05_novelty/{sample}_combined_novelty_stats.json"
    params:
        outdir = "{results}/05_novelty",
        sample = "{sample}",
        homology_weight = config["novelty_detection"]["integration"]["homology_weight"],
        ml_weight = config["novelty_detection"]["integration"]["ml_weight"],
        confidence_boost = config["novelty_detection"]["integration"]["confidence_boost"],
        agreement_threshold = config["novelty_detection"]["integration"]["agreement_threshold"],
        plots = "--plots" if config.get("advanced", {}).get("detailed_logs", False) else ""
    resources:
        mem_mb = 16000,
        runtime = 60
    log:
        "{logs}/novelty/{sample}_combine.log"
    conda:
        "environment.yml"
    shell:
        """
        python scripts/combine_novelty_results.py \
            --homology-results {input.homology} \
            --ml-results {input.ml} \
            --output {params.outdir} \
            --sample-name {params.sample} \
            --homology-weight {params.homology_weight} \
            --ml-weight {params.ml_weight} \
            --confidence-boost {params.confidence_boost} \
            --agreement-threshold {params.agreement_threshold} \
            {params.plots} \
            > {log} 2>&1
        """

# =============================================================================
# Abundance Calculation Rules
# =============================================================================

rule contig_abundance:
    """Calculate contig-level abundance"""
    input:
        assembly = "{results}/02_assembly/{sample}_contigs_filtered.fasta",
        r1 = "{results}/01_quality_control/{sample}_trimmed_R1.fastq.gz",
        r2 = "{results}/01_quality_control/{sample}_trimmed_R2.fastq.gz"
    output:
        abundance = "{results}/06_abundance/{sample}_contig_abundance.tsv",
        matrix = "{results}/06_abundance/{sample}_abundance_matrix.tsv",
        stats = "{results}/06_abundance/{sample}_abundance_stats.json"
    params:
        outdir = "{results}/06_abundance",
        sample = "{sample}",
        mapper = "bwa",
        normalization = config["abundance"]["normalization"]["method"],
        min_mapq = config["abundance"]["coverage"]["min_mapq"],
        min_baseq = config["abundance"]["coverage"]["min_baseq"],
        threads = config["abundance"]["mapping"]["bwa"]["threads"],
        plots = "--plots" if config.get("advanced", {}).get("detailed_logs", False) else ""
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        runtime = 240
    log:
        "{logs}/abundance/{sample}_contig.log"
    conda:
        "environment.yml"
    shell:
        """
        python scripts/calculate_contig_abundance.py \
            --assembly {input.assembly} \
            --reads1 {input.r1} \
            --reads2 {input.r2} \
            --output {params.outdir} \
            --sample-name {params.sample} \
            --mapper {params.mapper} \
            --normalization {params.normalization} \
            --min-mapq {params.min_mapq} \
            --min-baseq {params.min_baseq} \
            --threads {params.threads} \
            {params.plots} \
            > {log} 2>&1
        """

rule mag_abundance:
    """Calculate MAG-level abundance"""
    input:
        contig_abundance = "{results}/06_abundance/{sample}_contig_abundance.tsv",
        bins_dir = "{results}/03_binning/{sample}/final_bins",
        bin_quality = "{results}/03_binning/{sample}_bin_quality.tsv",
        taxonomy = "{results}/04_taxonomy/{sample}_gtdbtk_summary.tsv"
    output:
        abundance = "{results}/06_abundance/{sample}_mag_abundance.tsv",
        high_quality = "{results}/06_abundance/{sample}_high_quality_mag_abundance.tsv",
        matrix = "{results}/06_abundance/{sample}_mag_abundance_matrix.tsv",
        stats = "{results}/06_abundance/{sample}_mag_abundance_stats.json"
    params:
        outdir = "{results}/06_abundance",
        sample = "{sample}",
        min_completeness = config["binning"]["quality_thresholds"]["min_completeness"],
        max_contamination = config["binning"]["quality_thresholds"]["max_contamination"],
        normalization = config["abundance"]["normalization"]["method"],
        plots = "--plots" if config.get("advanced", {}).get("detailed_logs", False) else ""
    resources:
        mem_mb = 16000,
        runtime = 60
    log:
        "{logs}/abundance/{sample}_mag.log"
    conda:
        "environment.yml"
    shell:
        """
        python scripts/calculate_mag_abundance.py \
            --contig-abundance {input.contig_abundance} \
            --bins-dir {input.bins_dir} \
            --bin-quality {input.bin_quality} \
            --taxonomy-file {input.taxonomy} \
            --output {params.outdir} \
            --sample-name {params.sample} \
            --min-completeness {params.min_completeness} \
            --max-contamination {params.max_contamination} \
            --normalization {params.normalization} \
            --quality-weighting \
            {params.plots} \
            > {log} 2>&1
        """

# =============================================================================
# Reporting Rules
# =============================================================================

rule phase1_report:
    """Generate comprehensive Phase 1 report"""
    input:
        assembly_stats = "{results}/02_assembly/{sample}_assembly_stats.json",
        bin_quality = "{results}/03_binning/{sample}_bin_quality.tsv",
        taxonomy = "{results}/04_taxonomy/{sample}_gtdbtk_summary.tsv",
        novelty = "{results}/05_novelty/{sample}_combined_novelty.tsv",
        contig_abundance = "{results}/06_abundance/{sample}_contig_abundance.tsv",
        mag_abundance = "{results}/06_abundance/{sample}_mag_abundance.tsv"
    output:
        report = "{results}/07_reports/{sample}_phase1_summary.html",
        json = "{results}/07_reports/{sample}_phase1_summary.json"
    params:
        outdir = "{results}/07_reports",
        sample = "{sample}",
        results_dir = "{results}"
    resources:
        mem_mb = 8000,
        runtime = 30
    log:
        "{logs}/reports/{sample}_report.log"
    conda:
        "environment.yml"
    shell:
        """
        python scripts/generate_phase1_report.py \
            --sample-name {params.sample} \
            --results-dir {params.results_dir} \
            --output {params.outdir} \
            > {log} 2>&1
        """

# =============================================================================
# Utility Rules
# =============================================================================

rule clean:
    """Clean temporary files"""
    shell:
        """
        rm -rf {TEMP_DIR}/*
        find {RESULTS_DIR} -name "*.tmp" -delete
        find {RESULTS_DIR} -name "*.temp" -delete
        """

rule clean_intermediate:
    """Clean intermediate files"""
    shell:
        """
        find {RESULTS_DIR} -name "*_mapped.bam" -delete
        find {RESULTS_DIR} -name "*_mapped.bam.bai" -delete
        find {RESULTS_DIR} -path "*/metaspades/*" -delete
        """

# =============================================================================
# Error Handling and Logging
# =============================================================================

onerror:
    """Handle pipeline errors"""
    print(f"Pipeline failed. Check logs in {LOGS_DIR}")
    shell("find {LOGS_DIR} -name '*.log' -exec tail -20 {{}} \\;")

onsuccess:
    """Handle successful completion"""
    print("🧬 Phase 1 pipeline completed successfully!")
    print(f"Results available in: {RESULTS_DIR}")
    
    if config.get("advanced", {}).get("cleanup_temp_files", False):
        shell("snakemake clean_intermediate")

# =============================================================================
# Configuration Validation
# =============================================================================

def validate_config():
    """Validate configuration parameters"""
    required_sections = ["samples", "directories", "assembly", "binning", 
                        "taxonomy", "novelty_detection", "abundance"]
    
    for section in required_sections:
        if section not in config:
            raise ValueError(f"Missing required config section: {section}")
    
    # Validate database paths
    databases = config.get("databases", {})
    for db_name, db_path in databases.items():
        if db_path and not os.path.exists(db_path):
            print(f"Warning: Database not found: {db_name} at {db_path}")

# Run validation
validate_config()


