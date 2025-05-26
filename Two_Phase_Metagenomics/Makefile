# Metagenomics Pipeline Makefile
# Provides convenient commands for setup, validation, and execution

.PHONY: help setup validate test run clean plot report all
.DEFAULT_GOAL := help

# Configuration
PIPELINE_DIR := $(shell pwd)
CONDA_ENV := metagenomics-pipeline
CONFIG_FILE := config/config.yaml
CORES := 32
MEMORY := 200GB

# Colors for output
RED := \033[0;31m
GREEN := \033[0;32m
YELLOW := \033[1;33m
BLUE := \033[0;34m
NC := \033[0m # No Color

# Help target
help: ## Show this help message
	@echo "$(BLUE)Metagenomics Pipeline - Make Commands$(NC)"
	@echo "======================================"
	@echo ""
	@echo "$(GREEN)Setup Commands:$(NC)"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; /^[a-zA-Z_-]+:.*?## .*$$/ {printf "  $(YELLOW)%-15s$(NC) %s\n", $$1, $$2}' | \
		grep -E "(setup|install|environment|database)"
	@echo ""
	@echo "$(GREEN)Validation Commands:$(NC)"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; /^[a-zA-Z_-]+:.*?## .*$$/ {printf "  $(YELLOW)%-15s$(NC) %s\n", $$1, $$2}' | \
		grep -E "(validate|test|check)"
	@echo ""
	@echo "$(GREEN)Execution Commands:$(NC)"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; /^[a-zA-Z_-]+:.*?## .*$$/ {printf "  $(YELLOW)%-15s$(NC) %s\n", $$1, $$2}' | \
		grep -E "(run|dry-run|cluster)"
	@echo ""
	@echo "$(GREEN)Analysis Commands:$(NC)"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; /^[a-zA-Z_-]+:.*?## .*$$/ {printf "  $(YELLOW)%-15s$(NC) %s\n", $$1, $$2}' | \
		grep -E "(plot|report|clean)"
	@echo ""
	@echo "$(GREEN)Combined Commands:$(NC)"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; /^[a-zA-Z_-]+:.*?## .*$$/ {printf "  $(YELLOW)%-15s$(NC) %s\n", $$1, $$2}' | \
		grep -E "(all|quick)"

# Setup commands
setup: ## Complete pipeline setup (environment + databases)
	@echo "$(BLUE)Starting complete pipeline setup...$(NC)"
	@chmod +x setup.sh
	@./setup.sh
	@echo "$(GREEN)Setup completed!$(NC)"

environment: ## Create conda environment only
	@echo "$(BLUE)Creating conda environment...$(NC)"
	@chmod +x setup.sh
	@./setup.sh --env-only
	@echo "$(GREEN)Environment created!$(NC)"

databases: ## Download databases only
	@echo "$(BLUE)Downloading databases...$(NC)"
	@chmod +x setup.sh
	@./setup.sh --db-only
	@echo "$(GREEN)Databases downloaded!$(NC)"

install-deps: ## Install additional dependencies
	@echo "$(BLUE)Installing additional dependencies...$(NC)"
	@conda activate $(CONDA_ENV) && \
		pip install --upgrade torch transformers && \
		pip install checkm2 cat-bat
	@echo "$(GREEN)Dependencies installed!$(NC)"

# Validation commands
validate: ## Validate pipeline setup and dependencies
	@echo "$(BLUE)Validating pipeline setup...$(NC)"
	@python validate_pipeline.py --pipeline-dir $(PIPELINE_DIR)
	@echo "$(GREEN)Validation completed! Check validation_report.txt$(NC)"

test: ## Run pipeline test with synthetic data
	@echo "$(BLUE)Running pipeline test...$(NC)"
	@python validate_pipeline.py --include-test --pipeline-dir $(PIPELINE_DIR)
	@echo "$(GREEN)Test completed!$(NC)"

check-syntax: ## Check Snakemake syntax only
	@echo "$(BLUE)Checking pipeline syntax...$(NC)"
	@snakemake --configfile $(CONFIG_FILE) --dry-run --quiet
	@echo "$(GREEN)Syntax check passed!$(NC)"

check-config: ## Validate configuration file
	@echo "$(BLUE)Checking configuration...$(NC)"
	@python -c "import yaml; yaml.safe_load(open('$(CONFIG_FILE)'))"
	@echo "$(GREEN)Configuration is valid!$(NC)"

# Execution commands
dry-run: ## Show what the pipeline would do without executing
	@echo "$(BLUE)Running pipeline dry-run...$(NC)"
	@chmod +x run_pipeline.sh
	@./run_pipeline.sh --dry-run --cores $(CORES)

run: ## Run the complete pipeline
	@echo "$(BLUE)Running complete pipeline...$(NC)"
	@chmod +x run_pipeline.sh
	@./run_pipeline.sh --cores $(CORES) --memory $(MEMORY)
	@echo "$(GREEN)Pipeline completed!$(NC)"

run-local: ## Run pipeline locally (alias for run)
	@$(MAKE) run

run-cluster: ## Run pipeline on SLURM cluster
	@echo "$(BLUE)Running pipeline on cluster...$(NC)"
	@chmod +x run_pipeline.sh
	@./run_pipeline.sh --cluster --cores 128
	@echo "$(GREEN)Pipeline submitted to cluster!$(NC)"

unlock: ## Unlock Snakemake working directory
	@echo "$(BLUE)Unlocking working directory...$(NC)"
	@./run_pipeline.sh --unlock
	@echo "$(GREEN)Directory unlocked!$(NC)"

# Individual pipeline steps
qc: ## Run quality control only
	@echo "$(BLUE)Running quality control...$(NC)"
	@snakemake --configfile $(CONFIG_FILE) --cores $(CORES) \
		results/01_quality_control/*/

assembly: ## Run assembly only
	@echo "$(BLUE)Running assembly...$(NC)"
	@snakemake --configfile $(CONFIG_FILE) --cores $(CORES) \
		results/02_assembly/*/

binning: ## Run binning only
	@echo "$(BLUE)Running binning...$(NC)"
	@snakemake --configfile $(CONFIG_FILE) --cores $(CORES) \
		results/03_binning/*/

taxonomy: ## Run taxonomic classification only
	@echo "$(BLUE)Running taxonomic classification...$(NC)"
	@snakemake --configfile $(CONFIG_FILE) --cores $(CORES) \
		results/04_taxonomy/*/

novelty: ## Run novelty detection only
	@echo "$(BLUE)Running novelty detection...$(NC)"
	@snakemake --configfile $(CONFIG_FILE) --cores $(CORES) \
		results/05_novelty/*/

abundance: ## Run abundance calculation only
	@echo "$(BLUE)Running abundance calculation...$(NC)"
	@snakemake --configfile $(CONFIG_FILE) --cores $(CORES) \
		results/06_abundance/*/

# Analysis commands
plot: ## Generate visualization plots
	@echo "$(BLUE)Generating plots...$(NC)"
	@conda activate $(CONDA_ENV) && python visualize_results.py
	@echo "$(GREEN)Plots generated! Check plots/ directory$(NC)"

report: ## Generate Snakemake HTML report
	@echo "$(BLUE)Generating pipeline report...$(NC)"
	@snakemake --configfile $(CONFIG_FILE) --report results/pipeline_report.html
	@echo "$(GREEN)Report generated: results/pipeline_report.html$(NC)"

summary: ## Generate summary statistics
	@echo "$(BLUE)Generating summary statistics...$(NC)"
	@conda activate $(CONDA_ENV) && python visualize_results.py --plot-type all
	@echo "$(GREEN)Summary completed! Check plots/summary_report.md$(NC)"

# Maintenance commands
clean: ## Clean intermediate files (keep final results)
	@echo "$(BLUE)Cleaning intermediate files...$(NC)"
	@rm -rf .snakemake/
	@rm -rf logs/slurm-*.out logs/slurm-*.err
	@find results/ -name "*.tmp" -delete
	@find results/ -name "*.temp" -delete
	@echo "$(GREEN)Cleaned intermediate files!$(NC)"

clean-all: ## Clean all results and outputs
	@echo "$(RED)Warning: This will delete ALL results!$(NC)"
	@read -p "Are you sure? [y/N] " -n 1 -r; \
	if [[ $$REPLY =~ ^[Yy]$$ ]]; then \
		echo ""; \
		echo "$(BLUE)Cleaning all results...$(NC)"; \
		rm -rf results/ plots/ test_results/ test_data/; \
		rm -f validation_report.txt pipeline_stats.json; \
		echo "$(GREEN)All results cleaned!$(NC)"; \
	else \
		echo ""; \
		echo "$(YELLOW)Cancelled.$(NC)"; \
	fi

clean-logs: ## Clean log files
	@echo "$(BLUE)Cleaning log files...$(NC)"
	@rm -rf logs/
	@mkdir -p logs
	@echo "$(GREEN)Log files cleaned!$(NC)"

# Combined workflows
quick-start: ## Quick setup and validation for new users
	@echo "$(BLUE)Quick start setup...$(NC)"
	@$(MAKE) environment
	@$(MAKE) validate
	@echo "$(GREEN)Quick start completed!$(NC)"
	@echo "$(YELLOW)Next steps:$(NC)"
	@echo "1. Edit config/config.yaml with your sample paths"
	@echo "2. Run: make dry-run"
	@echo "3. Run: make run"

full-setup: ## Complete setup with test run
	@echo "$(BLUE)Full setup with testing...$(NC)"
	@$(MAKE) setup
	@$(MAKE) test
	@echo "$(GREEN)Full setup completed!$(NC)"

all: ## Run complete pipeline workflow
	@echo "$(BLUE)Running complete workflow...$(NC)"
	@$(MAKE) validate
	@$(MAKE) run
	@$(MAKE) plot
	@$(MAKE) report
	@echo "$(GREEN)Complete workflow finished!$(NC)"

# Development commands
dev-setup: ## Setup for development
	@echo "$(BLUE)Setting up development environment...$(NC)"
	@$(MAKE) environment
	@conda activate $(CONDA_ENV) && \
		pip install pytest black flake8 mypy pre-commit
	@echo "$(GREEN)Development environment ready!$(NC)"

lint: ## Run code linting
	@echo "$(BLUE)Running code linting...$(NC)"
	@conda activate $(CONDA_ENV) && \
		black scripts/*.py *.py && \
		flake8 scripts/*.py *.py
	@echo "$(GREEN)Linting completed!$(NC)"

# Status and monitoring
status: ## Show pipeline status
	@echo "$(BLUE)Pipeline Status$(NC)"
	@echo "==============="
	@echo "Pipeline directory: $(PIPELINE_DIR)"
	@echo "Configuration file: $(CONFIG_FILE)"
	@echo "Conda environment: $(CONDA_ENV)"
	@echo ""
	@echo "$(GREEN)Available samples:$(NC)"
	@if [ -f $(CONFIG_FILE) ]; then \
		python -c "import yaml; config=yaml.safe_load(open('$(CONFIG_FILE)')); [print(f'  - {s}') for s in config.get('samples', {}).keys()]"; \
	else \
		echo "  No configuration file found"; \
	fi
	@echo ""
	@echo "$(GREEN)Results status:$(NC)"
	@if [ -d results ]; then \
		echo "  Results directory: ✓"; \
		find results -name "*.tsv" | wc -l | xargs echo "  Output files: "; \
	else \
		echo "  Results directory: ✗"; \
	fi

monitor: ## Monitor running pipeline
	@echo "$(BLUE)Monitoring pipeline...$(NC)"
	@watch -n 30 'ls -la results/*/; echo ""; tail -5 logs/slurm-*.out 2>/dev/null || echo "No cluster logs found"'

# Information commands
info: ## Show pipeline information
	@echo "$(BLUE)Metagenomics Pipeline Information$(NC)"
	@echo "================================="
	@echo ""
	@echo "$(GREEN)Pipeline Features:$(NC)"
	@echo "  • Assembly-first approach for comprehensive analysis"
	@echo "  • Dual novelty detection (homology + ML-based)"
	@echo "  • Multi-algorithm binning with refinement"
	@echo "  • Contig-based community structure analysis"
	@echo "  • DNABert + Isolation Forest for pattern detection"
	@echo ""
	@echo "$(GREEN)Pipeline Steps:$(NC)"
	@echo "  1. Quality control & preprocessing"
	@echo "  2. De novo assembly & binning"
	@echo "  3. Taxonomic classification"
	@echo "  4. Novelty detection (homology + ML)"
	@echo "  5. Abundance estimation"
	@echo ""
	@echo "$(GREEN)Output Structure:$(NC)"
	@echo "  results/"
	@echo "    ├── 01_quality_control/"
	@echo "    ├── 02_assembly/"
	@echo "    ├── 03_binning/"
	@echo "    ├── 04_taxonomy/"
	@echo "    ├── 05_novelty/"
	@echo "    └── 06_abundance/"

docs: ## Open documentation
	@echo "$(BLUE)Opening documentation...$(NC)"
	@if command -v xdg-open >/dev/null; then \
		xdg-open README.md; \
	elif command -v open >/dev/null; then \
		open README.md; \
	else \
		echo "Please open README.md manually"; \
	fi

# Example configurations
example-config: ## Generate example configuration file
	@echo "$(BLUE)Generating example configuration...$(NC)"
	@mkdir -p config
	@cat > config/example_config.yaml << 'EOF'
# Example Configuration for Metagenomics Pipeline

samples:
  sample_A:
    R1: "data/sample_A_R1.fastq.gz"
    R2: "data/sample_A_R2.fastq.gz"
  sample_B:
    R1: "data/sample_B_R1.fastq.gz"
    R2: "data/sample_B_R2.fastq.gz"

output_dir: "results"
threads: 32
memory: 200

databases:
  gtdb: "databases/gtdb/release214"
  host_genome: "databases/host_genomes/human_genome"
  ncbi: "databases/ncbi/nr.dmnd"
  checkm2: "databases/checkm2/CheckM2_database"

# See config.yaml for full configuration options
EOF
	@echo "$(GREEN)Example configuration created: config/example_config.yaml$(NC)"

# Version information
version: ## Show version information
	@echo "$(BLUE)Version Information$(NC)"
	@echo "==================="
	@echo "Pipeline version: 1.0.0"
	@echo "Snakemake: $$(snakemake --version 2>/dev/null || echo 'Not installed')"
	@echo "Python: $$(python --version 2>/dev/null || echo 'Not installed')"
	@echo "Conda: $$(conda --version 2>/dev/null || echo 'Not installed')"
	@echo ""
	@echo "For detailed dependency versions, run: make validate"

