# Metagenomics Shotgun Pipeline - Single Sample Processing

A comprehensive Snakemake pipeline for metagenomics shotgun data analysis focusing on **community structure** and **novelty detection** using an assembly-first approach.

## Key Features

- **Assembly-First Approach**: Prioritizes de novo assembly for comprehensive novelty detection
- **Dual Novelty Detection**: Combines traditional homology-based methods with ML-based approaches (DNABert + Isolation Forest)
- **Comprehensive Binning**: Uses multiple binning algorithms with DAS Tool refinement
- **Quality Control**: Extensive QC throughout the pipeline
- **Contig-Based Analysis**: All diversity and abundance calculations based on assembled contigs

## Pipeline Overview

### Single Sample Processing Workflow

1. **Quality Control & Preprocessing**
   - FastQC quality assessment
   - Adapter trimming with fastp
   - Host decontamination

2. **Assembly & Binning**
   - MetaSPAdes de novo assembly
   - Multi-algorithm binning (MetaBAT2, CONCOCT, MaxBin2)
   - DAS Tool bin refinement
   - CheckM2 quality assessment

3. **Taxonomic Profiling** (Assembly-Based)
   - GTDB-Tk classification of MAGs
   - CAT/BAT contig-level classification
   - Phylogenetic placement for novel lineages

4. **Enhanced Novelty Detection**
   - **Homology-based**: DIAMOND BLAST vs NCBI database
   - **ML-based**: DNABert embeddings + Isolation Forest outlier detection
   - **Integrated scoring**: Combined novelty assessment

5. **Abundance Estimation**
   - Read mapping to contigs
   - Coverage-based abundance calculation
   - MAG-level abundance aggregation

## Installation & Setup

### Requirements

```bash
# Core dependencies
snakemake>=7.0
conda/mamba
python>=3.8

# Bioinformatics tools
metaspades
fastqc
fastp
bowtie2
samtools
diamond
gtdbtk
checkm2
metabat2
concoct
maxbin2
dastool

# Python packages
biopython
pandas
numpy
scikit-learn
torch
transformers
pysam
```

### Database Setup

1. **GTDB Database**
```bash
# Download GTDB-Tk database
wget https://data.gtdb.ecogenomic.org/releases/release214/auxillary_files/gtdbtk_package/full_package/gtdbtk_r214_data.tar.gz
tar -xzf gtdbtk_r214_data.tar.gz
```

2. **NCBI NR Database**
```bash
# Download and prepare DIAMOND database
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in nr.gz --db nr
```

3. **CheckM2 Database**
```bash
# Download CheckM2 database
checkm2 database --download --path /path/to/checkm2_db
```

4. **Host Genome Index**
```bash
# Build Bowtie2 index for host genome (e.g., human)
bowtie2-build human_genome.fasta human_genome
```

### Configuration

1. **Edit `config/config.yaml`**:
```yaml
# Update sample paths
samples:
  sample1:
    R1: "data/sample1_R1.fastq.gz"
    R2: "data/sample1_R2.fastq.gz"

# Update database paths
databases:
  gtdb: "/path/to/gtdb/release214"
  host_genome: "/path/to/host/genome/index"
  ncbi: "/path/to/ncbi_nr/nr.dmnd"
  checkm2: "/path/to/checkm2/CheckM2_database"
```

## Usage

### Quick Start

```bash
# Clone the pipeline
git clone https://github.com/your-repo/metagenomics-pipeline
cd metagenomics-pipeline

# Edit configuration
vim config/config.yaml

# Dry run to check pipeline
./run_pipeline.sh --dry-run

# Run pipeline
./run_pipeline.sh --cores 32 --memory 200GB
```

### Advanced Usage

```bash
# Run on SLURM cluster
./run_pipeline.sh --cluster --cores 128

# Custom configuration
./run_pipeline.sh --config custom_config.yaml

# Unlock if needed
./run_pipeline.sh --unlock
```

### Individual Steps

```bash
# Run specific rules
snakemake --cores 8 results/01_quality_control/sample1/sample1_trimmed_R1.fastq.gz
snakemake --cores 32 results/02_assembly/sample1/sample1_contigs.fasta
snakemake --cores 8 results/05_novelty/sample1/combined_novelty_report.tsv
```

## Output Structure

```
results/
├── 01_quality_control/
│   └── sample1/
│       ├── sample1_R1_fastqc.html
│       ├── sample1_trimmed_R1.fastq.gz
│       └── sample1_decontaminated_R1.fastq.gz
├── 02_assembly/
│   └── sample1/
│       ├── sample1_contigs.fasta
│       ├── sample1_assembly_stats.txt
│       └── sample1_mapped.bam
├── 03_binning/
│   └── sample1/
│       ├── refined_bins/
│       └── checkm2_results/
├── 04_taxonomy/
│   └── sample1/
│       ├── gtdbtk/
│       └── cat_bat/
├── 05_novelty/
│   └── sample1/
│       ├── homology_novelty.tsv
│       ├── ml_novelty_scores.tsv
│       └── combined_novelty_report.tsv
└── 06_abundance/
    └── sample1/
        ├── contig_abundance.tsv
        └── mag_abundance.tsv
```

## Key Output Files

### Novelty Detection Results

**`combined_novelty_report.tsv`** - Main novelty detection output:
- `integrated_novelty_score`: Combined score from both methods
- `integrated_novelty_class`: High/moderate/low confidence classifications
- `method_agreement`: Agreement between homology and ML methods
- Individual method scores and classifications

### Community Structure

**`contig_abundance.tsv`** - Contig-level abundance:
- Coverage statistics, relative abundance, TPM values

**`mag_abundance.tsv`** - MAG-level abundance:
- Aggregated coverage, genome completeness, relative abundance

### Assembly Quality

**`assembly_stats.txt`** - Assembly metrics:
- N50/N90, total length, contig distribution, GC content

## Novelty Detection Methods

### 1. Homology-Based Detection
- DIAMOND BLAST against NCBI NR database
- Identifies sequences with low identity/coverage to known sequences
- Classifies: Known, Moderate novelty, High novelty, No significant hit

### 2. ML-Based Detection (DNABert + Isolation Forest)
- Generates sequence embeddings using pre-trained DNABert model
- Applies Isolation Forest for unsupervised outlier detection
- Captures novel sequence patterns independent of homology

### 3. Integrated Scoring
- Combines both methods with weighted scoring
- Provides confidence levels based on method agreement
- Enables discovery of different types of novelty

## Performance Considerations

### Computational Requirements
- **Memory**: 200+ GB recommended for large datasets
- **CPU**: 32+ cores for parallel processing
- **GPU**: Optional for DNABert (speeds up embedding generation)
- **Storage**: ~50-100GB per sample for intermediate files

### Runtime Estimates
- Small sample (1-5GB): 4-8 hours
- Medium sample (5-20GB): 8-24 hours  
- Large sample (20-50GB): 24-48 hours

### Optimization Tips
- Use SSD storage for faster I/O
- Enable GPU for ML novelty detection
- Adjust memory allocation per tool in config
- Use cluster mode for multiple samples

## Troubleshooting

### Common Issues

1. **Memory errors during assembly**
   - Increase memory allocation in config
   - Use `--careful` mode for MetaSPAdes

2. **GTDB-Tk classification fails**
   - Check GTDB database path and version
   - Ensure sufficient memory (100+ GB)

3. **DNABert out of memory**
   - Reduce chunk size in ML novelty detection
   - Use GPU with more memory

4. **No bins generated**
   - Check assembly quality (N50, total length)
   - Verify read mapping success
   - Adjust binning parameters

### Getting Help

- Check Snakemake logs in `logs/` directory
- Review individual tool logs
- Ensure all databases are properly installed
- Verify input file formats and paths

## Next Steps

After single sample processing:

1. **Cross-Sample Analysis**: Compare novelty patterns across samples
2. **Functional Annotation**: Annotate novel sequences with function
3. **Phylogenetic Analysis**: Place novel MAGs in phylogenetic context
4. **Comparative Genomics**: Analyze novel gene content and pathways

## Citation

If you use this pipeline, please cite:
- Snakemake workflow management system
- Individual tools used (MetaSPAdes, GTDB-Tk, etc.)
- DNABert model for sequence embeddings
- Relevant databases (GTDB, NCBI NR)

## License

This pipeline is released under the MIT License. See LICENSE file for details.
