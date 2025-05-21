
SAMPLES = ["sample1", "sample2", "sample3"]
RAW_DIR = "raw/"
QC_DIR = "qc/"
ID_DIR = "identification/"
QUANT_DIR = "quantification/"
NORMAL_DIR = "normalized/"
RESULTS_DIR = "results/"

rule all:
    input:
        expand(f"{RESULTS_DIR}{{sample}}_de_results.csv", sample=SAMPLES)

# 1. Convert raw data (if needed)
rule convert_raw_to_mzML:
    input:
        f"{RAW_DIR}{{sample}}.raw"
    output:
        f"{RAW_DIR}{{sample}}.mzML"
    shell:
        "msconvert {input} --mzML -o {RAW_DIR}"

# 2. Quality Control
rule run_rawtools_qc:
    input:
        f"{RAW_DIR}{{sample}}.mzML"
    output:
        f"{QC_DIR}{{sample}}_qc.html"
    shell:
        "RawToolsQC --input {input} --output {output}"

# 3. Identification using MaxQuant
rule run_maxquant:
    input:
        mzml=expand(f"{RAW_DIR}{{sample}}.mzML", sample=SAMPLES),
        fasta="db/uniprot.fasta"
    output:
        f"{ID_DIR}/combined/txt/proteinGroups.txt"
    shell:
        "MaxQuantCmd.exe mqpar.xml"

# 4. Protein Quantification (from MaxQuant output)
rule extract_protein_matrix:
    input:
        f"{ID_DIR}/combined/txt/proteinGroups.txt"
    output:
        f"{QUANT_DIR}/protein_abundance_matrix.csv"
    script:
        "scripts/extract_quant.py"

# 5. Normalization and Filtering
rule normalize_and_filter:
    input:
        f"{QUANT_DIR}/protein_abundance_matrix.csv"
    output:
        f"{NORMAL_DIR}/normalized_filtered.csv"
    script:
        "scripts/normalize_filter.py"

# 6. Statistical and ML Analysis
rule differential_expression:
    input:
        f"{NORMAL_DIR}/normalized_filtered.csv"
    output:
        f"{RESULTS_DIR}{{sample}}_de_results.csv"
    script:
        "scripts/differential_expression.R"

# 7. Outliers detection
rule detect_outliers:
    input:
        f"{NORMAL_DIR}/normalized_filtered.csv"
    output:
        f"{RESULTS_DIR}{{sample}}_outliers.csv"
    script:
        "scripts/detect_outliers.py"


