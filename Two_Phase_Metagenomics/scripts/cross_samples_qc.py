#!/bin/bash

# Cross-Sample Metagenomics Analysis Runner
# Execute comprehensive cross-sample analysis pipeline

set -euo pipefail

# Default values
CONFIG_FILE="config/cross_sample_config.yaml"
SNAKEFILE="cross_sample_Snakefile"
CORES=8
DRY_RUN=false
UNLOCK=false
PROFILE=""
CLUSTER_CONFIG=""
TARGET="all"
VERBOSE=false
RESUME=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

# Help function
show_help() {
    cat << EOF
Cross-Sample Metagenomics Analysis Pipeline

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -c, --config FILE           Configuration file (default: $CONFIG_FILE)
    -j, --cores NUM            Number of cores to use (default: $CORES)
    -n, --dry-run              Perform a dry run (don't execute)
    -u, --unlock               Unlock working directory
    -p, --profile PROFILE      Use cluster profile
    --cluster-config FILE      Cluster configuration file
    -t, --target TARGET        Specific target to run (default: all)
    -v, --verbose              Verbose output
    -r, --resume               Resume interrupted analysis
    -h, --help                 Show this help message

TARGETS:
    all                        Complete cross-sample analysis
    abundance-only             Abundance integration only
    ml-only                    ML analysis only
    novelty-only               Novelty analysis only
    temporal-only              Temporal analysis only
    pathogen-only              Pathogen detection only

EXAMPLES:
    # Basic analysis with 16 cores
    $0 --cores 16

    # Dry run to check workflow
    $0 --dry-run

    # Run only ML analysis
    $0 --target ml-only

    # Run on cluster with SLURM
    $0 --profile slurm --cluster-config config/cluster.yaml

    # Resume interrupted analysis
    $0 --resume --cores 32

    # Novelty analysis with custom config
    $0 --config custom_config.yaml --target novelty-only

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -j|--cores)
            CORES="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -u|--unlock)
            UNLOCK=true
            shift
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        --cluster-config)
            CLUSTER_CONFIG="$2"
            shift 2
            ;;
        -t|--target)
            TARGET="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -r|--resume)
            RESUME=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Function to check dependencies
check_dependencies() {
    print_step "Checking dependencies..."
    
    local deps=("snakemake" "python3" "conda")
    local missing=()
    
    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            missing+=("$dep")
        fi
    done
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        print_error "Missing dependencies: ${missing[*]}"
        print_error "Please install missing dependencies and try again"
        exit 1
    fi
    
    print_status "All dependencies found"
}

# Function to validate configuration
validate_config() {
    print_step "Validating configuration..."
    
    if [[ ! -f "$CONFIG_FILE" ]]; then
        print_error "Configuration file not found: $CONFIG_FILE"
        exit 1
    fi
    
    if [[ ! -f "$SNAKEFILE" ]]; then
        print_error "Snakefile not found: $SNAKEFILE"
        exit 1
    fi
    
    # Check if required directories exist
    local input_dir
    input_dir=$(python3 -c "
import yaml
with open('$CONFIG_FILE', 'r') as f:
    config = yaml.safe_load(f)
    print(config.get('input_directory', ''))
")
    
    if [[ ! -d "$input_dir" ]]; then
        print_error "Input directory not found: $input_dir"
        print_error "Please ensure single-sample analysis has been completed"
        exit 1
    fi
    
    print_status "Configuration validated"
}

# Function to check input data
check_input_data() {
    print_step "Checking input data availability..."
    
    local python_check=$(cat << 'EOF'
import yaml
import os
import sys

with open(sys.argv[1], 'r') as f:
    config = yaml.safe_load(f)

input_dir = config.get('input_directory', '')
metadata_file = config.get('sample_metadata', '')

# Check metadata
if not os.path.exists(metadata_file):
    print(f"ERROR: Metadata file not found: {metadata_file}")
    sys.exit(1)

# Read sample list from metadata
import pandas as pd
try:
    metadata_df = pd.read_csv(metadata_file, sep='\t')
    samples = metadata_df['sample_id'].tolist()
    print(f"Found {len(samples)} samples in metadata")
except Exception as e:
    print(f"ERROR: Could not read metadata: {e}")
    sys.exit(1)

# Check for required files
missing_files = []
required_patterns = [
    'abundance/contig_abundance.tsv',
    'novelty/combined_novelty_results.tsv',
    'assembly/contigs.fasta'
]

for sample in samples[:5]:  # Check first 5 samples
    for pattern in required_patterns:
        file_path = os.path.join(input_dir, sample, pattern)
        if not os.path.exists(file_path):
            missing_files.append(file_path)

if missing_files:
    print(f"WARNING: Some expected files are missing (showing first 10):")
    for f in missing_files[:10]:
        print(f"  - {f}")
    if len(missing_files) > 10:
        print(f"  ... and {len(missing_files) - 10} more")
else:
    print("Input data check passed")
EOF
    )
    
    python3 -c "$python_check" "$CONFIG_FILE"
}

# Function to setup environment
setup_environment() {
    print_step "Setting up analysis environment..."
    
    # Activate conda environment if specified
    if [[ -f "environment.yml" ]]; then
        local env_name
        env_name=$(grep "name:" environment.yml | cut -d':' -f2 | xargs)
        
        if conda info --envs | grep -q "$env_name"; then
            print_status "Activating conda environment: $env_name"
            # Note: In practice, this should be sourced by the user
            print_warning "Please ensure conda environment '$env_name' is activated"
        else
            print_warning "Conda environment '$env_name' not found"
            print_warning "Please create environment: conda env create -f environment.yml"
        fi
    fi
    
    # Create output directories
    local output_dir
    output_dir=$(python3 -c "
import yaml
with open('$CONFIG_FILE', 'r') as f:
    config = yaml.safe_load(f)
    print(config.get('output_directory', 'results/cross_sample'))
")
    
    mkdir -p "$output_dir"/{abundance,mag_dereplication,ml_analysis,novelty_analysis,temporal_analysis,pathogen_detection,reports,logs}
    
    print_status "Environment setup complete"
}

# Function to build Snakemake command
build_snakemake_cmd() {
    local cmd="snakemake"
    
    # Add basic parameters
    cmd+=" --snakefile $SNAKEFILE"
    cmd+=" --configfile $CONFIG_FILE"
    cmd+=" --cores $CORES"
    cmd+=" --target $TARGET"
    
    # Add conditional parameters
    if [[ "$DRY_RUN" == true ]]; then
        cmd+=" --dry-run"
    fi
    
    if [[ "$VERBOSE" == true ]]; then
        cmd+=" --verbose"
    fi
    
    if [[ "$RESUME" == true ]]; then
        cmd+=" --rerun-incomplete"
    fi
    
    if [[ -n "$PROFILE" ]]; then
        cmd+=" --profile $PROFILE"
    fi
    
    if [[ -n "$CLUSTER_CONFIG" ]]; then
        cmd+=" --cluster-config $CLUSTER_CONFIG"
    fi
    
    # Add useful flags
    cmd+=" --use-conda"
    cmd+=" --conda-frontend conda"
    cmd+=" --printshellcmds"
    cmd+=" --reason"
    cmd+=" --stats stats.json"
    
    echo "$cmd"
}

# Function to run analysis
run_analysis() {
    print_step "Starting cross-sample analysis..."
    
    local cmd
    cmd=$(build_snakemake_cmd)
    
    print_status "Snakemake command: $cmd"
    
    if [[ "$DRY_RUN" == true ]]; then
        print_warning "DRY RUN MODE - No files will be created"
    fi
    
    # Execute command
    if eval "$cmd"; then
        if [[ "$DRY_RUN" == false ]]; then
            print_status "Analysis completed successfully!"
            show_results_summary
        else
            print_status "Dry run completed successfully!"
        fi
    else
        print_error "Analysis failed!"
        show_troubleshooting_tips
        exit 1
    fi
}

# Function to show results summary
show_results_summary() {
    print_step "Analysis Results Summary"
    
    local output_dir
    output_dir=$(python3 -c "
import yaml
with open('$CONFIG_FILE', 'r') as f:
    config = yaml.safe_load(f)
    print(config.get('output_directory', 'results/cross_sample'))
")
    
    echo "📊 Results available in: $output_dir"
    echo ""
    echo "📋 Key output files:"
    
    local key_files=(
        "reports/comprehensive_analysis_report.html:📄 Comprehensive HTML Report"
        "reports/executive_summary.pdf:📊 Executive Summary PDF"
        "pathogen_detection/high_risk_alerts.tsv:🚨 High-Risk Pathogen Alerts"
        "novelty_analysis/pan_novelty_results.tsv:🧬 Novel Features Detected"
        "ml_analysis/isolation_forest_results.tsv:🤖 ML Anomaly Detection"
        "temporal_analysis/temporal_trends.tsv:⏰ Temporal Trends"
    )
    
    for file_desc in "${key_files[@]}"; do
        local file_path="${file_desc%%:*}"
        local description="${file_desc##*:}"
        local full_path="$output_dir/$file_path"
        
        if [[ -f "$full_path" ]]; then
            local size
            size=$(du -h "$full_path" | cut -f1)
            echo "  ✅ $description ($size)"
            echo "     $full_path"
        else
            echo "  ❌ $description (not found)"
        fi
    done
    
    echo ""
    echo "🔗 Quick access:"
    echo "  View main report: open $output_dir/reports/comprehensive_analysis_report.html"
    echo "  Check high-risk alerts: cat $output_dir/pathogen_detection/high_risk_alerts.tsv"
    echo ""
}

# Function to show troubleshooting tips
show_troubleshooting_tips() {
    print_step "Troubleshooting Tips"
    
    echo "❓ Common issues and solutions:"
    echo ""
    echo "1. 🔒 Locked working directory:"
    echo "   Solution: $0 --unlock"
    echo ""
    echo "2. 📁 Missing input files:"
    echo "   Solution: Ensure single-sample analysis completed successfully"
    echo "   Check: ls -la results/single_sample/*/abundance/"
    echo ""
    echo "3. 💾 Insufficient memory:"
    echo "   Solution: Reduce cores or increase memory allocation"
    echo "   Example: $0 --cores 4"
    echo ""
    echo "4. 🐍 Conda environment issues:"
    echo "   Solution: conda env create -f environment.yml"
    echo "   Then: conda activate metagenomics-pipeline"
    echo ""
    echo "5. 📊 View detailed logs:"
    echo "   Check: tail -f results/cross_sample/logs/*.log"
    echo ""
    echo "6. 🔄 Resume interrupted analysis:"
    echo "   Solution: $0 --resume"
    echo ""
}

# Function to unlock directory
unlock_directory() {
    print_step "Unlocking working directory..."
    
    local cmd="snakemake --snakefile $SNAKEFILE --configfile $CONFIG_FILE --unlock"
    
    if eval "$cmd"; then
        print_status "Directory unlocked successfully"
    else
        print_error "Failed to unlock directory"
        exit 1
    fi
}

# Function to display banner
show_banner() {
    cat << 'EOF'
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║               🧬 Cross-Sample Metagenomics Analysis Pipeline 🧬              ║
║                                                                              ║
║                    Assembly-First with ML-Enhanced Detection                 ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝

EOF
}

# Main execution function
main() {
    show_banner
    
    print_status "Starting cross-sample metagenomics analysis"
    print_status "Target: $TARGET | Cores: $CORES | Config: $CONFIG_FILE"
    echo ""
    
    # Handle special cases
    if [[ "$UNLOCK" == true ]]; then
        unlock_directory
        exit 0
    fi
    
    # Run pre-flight checks
    check_dependencies
    validate_config
    check_input_data
    setup_environment
    
    # Run analysis
    run_analysis
    
    print_status "Cross-sample analysis pipeline completed!"
}

# Execute main function
main "$@"

