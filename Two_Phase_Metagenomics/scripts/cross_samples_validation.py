#!/usr/bin/env python3
"""
Cross-Sample Pipeline Validation Script
Comprehensive validation of pipeline setup and readiness
"""

import sys
import os
import yaml
import pandas as pd
import json
import subprocess
import importlib
from pathlib import Path
from collections import defaultdict
import argparse

# ANSI color codes
class Colors:
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'

def print_status(message, status="INFO"):
    """Print colored status message"""
    color_map = {
        "INFO": Colors.BLUE,
        "SUCCESS": Colors.GREEN,
        "WARNING": Colors.YELLOW,
        "ERROR": Colors.RED
    }
    color = color_map.get(status, Colors.BLUE)
    print(f"{color}[{status}]{Colors.END} {message}")

def print_header(title):
    """Print section header"""
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{title.center(60)}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.END}")

class CrossSampleValidator:
    def __init__(self, config_file="config/cross_sample_config.yaml"):
        self.config_file = config_file
        self.config = None
        self.validation_results = defaultdict(list)
        self.overall_status = True
        
    def load_config(self):
        """Load and validate configuration file"""
        print_header("CONFIGURATION VALIDATION")
        
        try:
            if not os.path.exists(self.config_file):
                self.add_error("Configuration file not found", f"Missing: {self.config_file}")
                return False
            
            with open(self.config_file, 'r') as f:
                self.config = yaml.safe_load(f)
            
            print_status(f"Configuration loaded: {self.config_file}", "SUCCESS")
            return True
            
        except yaml.YAMLError as e:
            self.add_error("Invalid YAML format", str(e))
            return False
        except Exception as e:
            self.add_error("Configuration loading failed", str(e))
            return False
    
    def validate_dependencies(self):
        """Validate required dependencies"""
        print_header("DEPENDENCY VALIDATION")
        
        # Check system dependencies
        system_deps = {
            'snakemake': 'Workflow management',
            'python3': 'Python interpreter',
            'conda': 'Package manager',
            'makeblastdb': 'BLAST database creation',
            'blastn': 'BLAST nucleotide search',
            'fastANI': 'Average nucleotide identity'
        }
        
        for dep, description in system_deps.items():
            if self.check_command(dep):
                print_status(f"{dep}: {description}", "SUCCESS")
            else:
                self.add_warning("Missing system dependency", f"{dep} ({description})")
        
        # Check Python packages
        python_packages = {
            'pandas': 'Data manipulation',
            'numpy': 'Numerical computing',
            'scipy': 'Scientific computing',
            'sklearn': 'Machine learning',
            'torch': 'Deep learning',
            'matplotlib': 'Plotting',
            'seaborn': 'Statistical plotting',
            'yaml': 'YAML parsing',
            'networkx': 'Graph analysis'
        }
        
        for package, description in python_packages.items():
            if self.check_python_package(package):
                print_status(f"Python: {package} ({description})", "SUCCESS")
            else:
                self.add_error("Missing Python package", f"{package} ({description})")
    
    def validate_file_structure(self):
        """Validate pipeline file structure"""
        print_header("FILE STRUCTURE VALIDATION")
        
        required_files = {
            'cross_sample_Snakefile': 'Main workflow file',
            'run_cross_sample_analysis.sh': 'Pipeline runner script',
            'Makefile': 'Make commands',
            'environment.yml': 'Conda environment'
        }
        
        for file_path, description in required_files.items():
            if os.path.exists(file_path):
                print_status(f"{file_path}: {description}", "SUCCESS")
            else:
                self.add_error("Missing pipeline file", f"{file_path} ({description})")
        
        # Check script directory
        scripts_dir = "scripts"
        if os.path.exists(scripts_dir):
            script_files = [
                'integrate_abundance_data.py',
                'mag_dereplication.py',
                'vae_community_embedding.py',
                'isolation_forest_cross_sample.py',
                'pan_novelty_analysis.py',
                'temporal_analysis.py',
                'pathogen_detection.py',
                'generate_comprehensive_report.py',
                'quality_control_cross_sample.py'
            ]
            
            for script in script_files:
                script_path = os.path.join(scripts_dir, script)
                if os.path.exists(script_path):
                    print_status(f"Script: {script}", "SUCCESS")
                else:
                    self.add_error("Missing script", script)
        else:
            self.add_error("Missing directory", "scripts/")
    
    def validate_configuration_content(self):
        """Validate configuration content"""
        print_header("CONFIGURATION CONTENT VALIDATION")
        
        if not self.config:
            self.add_error("Configuration validation skipped", "Config not loaded")
            return
        
        # Required configuration sections
        required_sections = {
            'input_directory': str,
            'output_directory': str,
            'sample_metadata': str,
            'mag_dereplication': dict,
            'vae_settings': dict,
            'isolation_forest': dict,
            'pathogen_detection': dict,
            'resources': dict
        }
        
        for section, expected_type in required_sections.items():
            if section in self.config:
                if isinstance(self.config[section], expected_type):
                    print_status(f"Config section: {section}", "SUCCESS")
                else:
                    self.add_warning("Config type mismatch", f"{section} should be {expected_type.__name__}")
            else:
                self.add_error("Missing config section", section)
        
        # Validate specific configuration values
        self.validate_config_values()
    
    def validate_config_values(self):
        """Validate specific configuration values"""
        
        # Check numeric ranges
        numeric_checks = [
            ('mag_dereplication.ani_threshold', 0.8, 1.0),
            ('vae_settings.latent_dim', 1, 1000),
            ('isolation_forest.contamination', 0.01, 0.5),
            ('resources.default_threads', 1, 128)
        ]
        
        for path, min_val, max_val in numeric_checks:
            value = self.get_nested_config(path)
            if value is not None:
                if min_val <= value <= max_val:
                    print_status(f"Config value: {path} = {value}", "SUCCESS")
                else:
                    self.add_warning("Config value out of range", f"{path} = {value} (expected {min_val}-{max_val})")
    
    def validate_input_data(self):
        """Validate input data availability"""
        print_header("INPUT DATA VALIDATION")
        
        if not self.config:
            self.add_error("Input validation skipped", "Config not loaded")
            return
        
        # Check metadata file
        metadata_file = self.config.get('sample_metadata')
        if metadata_file and os.path.exists(metadata_file):
            try:
                metadata_df = pd.read_csv(metadata_file, sep='\t')
                print_status(f"Metadata file: {len(metadata_df)} samples", "SUCCESS")
                
                # Validate metadata columns
                required_cols = ['sample_id']
                optional_cols = ['collection_date', 'site_id', 'treatment']
                
                for col in required_cols:
                    if col in metadata_df.columns:
                        print_status(f"Metadata column: {col}", "SUCCESS")
                    else:
                        self.add_error("Missing metadata column", col)
                
                for col in optional_cols:
                    if col in metadata_df.columns:
                        print_status(f"Optional metadata column: {col}", "SUCCESS")
                    else:
                        self.add_warning("Missing optional column", col)
                
                # Check for duplicates
                if metadata_df['sample_id'].duplicated().any():
                    self.add_error("Duplicate sample IDs", "in metadata file")
                
                # Validate single-sample results
                self.validate_single_sample_results(metadata_df['sample_id'].tolist())
                
            except Exception as e:
                self.add_error("Metadata validation failed", str(e))
        else:
            self.add_error("Metadata file not found", metadata_file or "not specified")
    
    def validate_single_sample_results(self, sample_ids):
        """Validate single-sample analysis results"""
        
        input_dir = self.config.get('input_directory', 'results/single_sample')
        
        if not os.path.exists(input_dir):
            self.add_error("Input directory not found", input_dir)
            return
        
        # Check a subset of samples
        check_samples = sample_ids[:min(5, len(sample_ids))]
        
        required_files = [
            'abundance/contig_abundance.tsv',
            'novelty/combined_novelty_results.tsv',
            'assembly/contigs.fasta'
        ]
        
        missing_files = []
        
        for sample_id in check_samples:
            for file_pattern in required_files:
                file_path = os.path.join(input_dir, sample_id, file_pattern)
                if os.path.exists(file_path):
                    continue
                else:
                    missing_files.append(f"{sample_id}/{file_pattern}")
        
        if missing_files:
            self.add_warning("Missing single-sample files", f"{len(missing_files)} files missing")
            for missing_file in missing_files[:10]:  # Show first 10
                self.add_warning("Missing file", missing_file)
        else:
            print_status(f"Single-sample results: {len(check_samples)} samples checked", "SUCCESS")
    
    def validate_databases(self):
        """Validate pathogen databases"""
        print_header("DATABASE VALIDATION")
        
        if not self.config:
            self.add_error("Database validation skipped", "Config not loaded")
            return
        
        databases = self.config.get('pathogen_detection', {}).get('databases', {})
        
        for db_name, db_path in databases.items():
            if db_path and os.path.exists(db_path):
                file_size = os.path.getsize(db_path) / (1024 * 1024)  # MB
                print_status(f"Database: {db_name} ({file_size:.1f} MB)", "SUCCESS")
            elif db_path:
                self.add_warning("Database not found", f"{db_name}: {db_path}")
            else:
                self.add_warning("Database path not set", db_name)
    
    def validate_output_permissions(self):
        """Validate output directory permissions"""
        print_header("OUTPUT VALIDATION")
        
        if not self.config:
            self.add_error("Output validation skipped", "Config not loaded")
            return
        
        output_dir = self.config.get('output_directory', 'results/cross_sample')
        
        # Try to create output directory
        try:
            os.makedirs(output_dir, exist_ok=True)
            print_status(f"Output directory: {output_dir}", "SUCCESS")
            
            # Test write permissions
            test_file = os.path.join(output_dir, '.test_write')
            with open(test_file, 'w') as f:
                f.write("test")
            os.remove(test_file)
            print_status("Write permissions: OK", "SUCCESS")
            
        except Exception as e:
            self.add_error("Output directory error", str(e))
    
    def validate_memory_requirements(self):
        """Validate memory requirements"""
        print_header("RESOURCE VALIDATION")
        
        if not self.config:
            self.add_error("Resource validation skipped", "Config not loaded")
            return
        
        memory_limits = self.config.get('resources', {}).get('memory_limits', {})
        
        # Check if memory requirements are reasonable
        total_memory_gb = 0
        for component, memory_gb in memory_limits.items():
            print_status(f"Memory requirement: {component} = {memory_gb} GB", "INFO")
            total_memory_gb = max(total_memory_gb, memory_gb)
        
        # Check system memory if possible
        try:
            if sys.platform.startswith('linux'):
                with open('/proc/meminfo', 'r') as f:
                    for line in f:
                        if line.startswith('MemTotal:'):
                            total_mem_kb = int(line.split()[1])
                            system_memory_gb = total_mem_kb / (1024 * 1024)
                            
                            if total_memory_gb <= system_memory_gb * 0.8:
                                print_status(f"System memory: {system_memory_gb:.1f} GB available", "SUCCESS")
                            else:
                                self.add_warning("Memory concern", f"Required: {total_memory_gb} GB, Available: {system_memory_gb:.1f} GB")
                            break
        except:
            print_status("System memory check: Skipped", "INFO")
    
    def check_command(self, command):
        """Check if command is available"""
        try:
            subprocess.run([command, '--version'], capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            try:
                subprocess.run([command, '-h'], capture_output=True, check=True)
                return True
            except (subprocess.CalledProcessError, FileNotFoundError):
                return False
    
    def check_python_package(self, package):
        """Check if Python package is available"""
        try:
            importlib.import_module(package)
            return True
        except ImportError:
            return False
    
    def get_nested_config(self, path):
        """Get nested configuration value"""
        keys = path.split('.')
        value = self.config
        
        for key in keys:
            if isinstance(value, dict) and key in value:
                value = value[key]
            else:
                return None
        
        return value
    
    def add_error(self, category, message):
        """Add error to validation results"""
        self.validation_results['errors'].append(f"{category}: {message}")
        self.overall_status = False
        print_status(f"{category}: {message}", "ERROR")
    
    def add_warning(self, category, message):
        """Add warning to validation results"""
        self.validation_results['warnings'].append(f"{category}: {message}")
        print_status(f"{category}: {message}", "WARNING")
    
    def add_success(self, category, message):
        """Add success to validation results"""
        self.validation_results['successes'].append(f"{category}: {message}")
        print_status(f"{category}: {message}", "SUCCESS")
    
    def generate_report(self):
        """Generate validation report"""
        print_header("VALIDATION SUMMARY")
        
        # Count results
        error_count = len(self.validation_results['errors'])
        warning_count = len(self.validation_results['warnings'])
        success_count = len(self.validation_results['successes'])
        
        # Overall status
        if self.overall_status:
            print_status("Overall Status: READY", "SUCCESS")
        else:
            print_status("Overall Status: ISSUES FOUND", "ERROR")
        
        print(f"\n📊 Results Summary:")
        print(f"   ✅ Successes: {success_count}")
        print(f"   ⚠️  Warnings:  {warning_count}")
        print(f"   ❌ Errors:    {error_count}")
        
        # Show errors if any
        if error_count > 0:
            print(f"\n{Colors.RED}❌ Critical Issues:{Colors.END}")
            for error in self.validation_results['errors']:
                print(f"   • {error}")
        
        # Show warnings if any
        if warning_count > 0:
            print(f"\n{Colors.YELLOW}⚠️  Warnings:{Colors.END}")
            for warning in self.validation_results['warnings'][:10]:  # Show first 10
                print(f"   • {warning}")
            if warning_count > 10:
                print(f"   • ... and {warning_count - 10} more warnings")
        
        # Recommendations
        print(f"\n💡 Recommendations:")
        if error_count > 0:
            print("   1. Fix critical errors before running pipeline")
            print("   2. Run setup script: ./setup_cross_sample.sh")
            print("   3. Check configuration: make validate-config")
        elif warning_count > 0:
            print("   1. Review warnings for optimal performance")
            print("   2. Consider running: make setup-env")
        else:
            print("   1. Pipeline is ready to run!")
            print("   2. Start analysis: make run")
        
        return self.overall_status
    
    def run_validation(self):
        """Run complete validation"""
        print_status("Starting cross-sample pipeline validation...", "INFO")
        
        # Run validation steps
        if self.load_config():
            self.validate_configuration_content()
        
        self.validate_dependencies()
        self.validate_file_structure()
        self.validate_input_data()
        self.validate_databases()
        self.validate_output_permissions()
        self.validate_memory_requirements()
        
        # Generate final report
        return self.generate_report()

def main():
    parser = argparse.ArgumentParser(
        description="Validate cross-sample metagenomics pipeline setup"
    )
    parser.add_argument(
        '--config', 
        default='config/cross_sample_config.yaml',
        help='Configuration file to validate'
    )
    parser.add_argument(
        '--output', 
        help='Save validation report to JSON file'
    )
    parser.add_argument(
        '--quiet', 
        action='store_true',
        help='Only show errors and warnings'
    )
    
    args = parser.parse_args()
    
    # Create validator
    validator = CrossSampleValidator(args.config)
    
    # Run validation
    success = validator.run_validation()
    
    # Save report if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump({
                'overall_status': success,
                'validation_results': dict(validator.validation_results),
                'timestamp': pd.Timestamp.now().isoformat()
            }, f, indent=2)
        print_status(f"Validation report saved: {args.output}", "INFO")
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()

