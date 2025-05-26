#!/usr/bin/env python3
"""
🔗 DAS Tool Input Preparation Utility
Part of the Metagenomics Pipeline - Phase 1

Prepares input files for DAS Tool bin refinement from multiple binning tools.
Creates scaffolds2bin files and validates bin quality.
"""

import argparse
import os
import sys
import logging
from pathlib import Path
from Bio import SeqIO
import json

def setup_logging(log_level="INFO"):
    """Setup logging configuration"""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Prepare input files for DAS Tool from binning results"
    )
    parser.add_argument(
        "--metabat2",
        help="MetaBAT2 bins list file"
    )
    parser.add_argument(
        "--maxbin2",
        help="MaxBin2 bins list file"
    )
    parser.add_argument(
        "--concoct",
        help="CONCOCT bins list file"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for DAS Tool input files"
    )
    parser.add_argument(
        "--min-contigs",
        type=int,
        default=5,
        help="Minimum number of contigs per bin (default: 5)"
    )
    parser.add_argument(
        "--min-bin-size",
        type=int,
        default=50000,
        help="Minimum total bin size in bp (default: 50,000)"
    )
    parser.add_argument(
        "--prefix",
        default="DAS_Tool_input",
        help="Prefix for output files (default: DAS_Tool_input)"
    )
    parser.add_argument(
        "--validate-contigs",
        help="Path to assembly file to validate contig IDs"
    )
    parser.add_argument(
        "--ignore-errors",
        action="store_true",
        help="Continue processing even with errors"
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )
    
    return parser.parse_args()

def read_bin_files(bin_list_file):
    """Read list of bin files and return valid ones"""
    if not bin_list_file or not os.path.exists(bin_list_file):
        return []
    
    bin_files = []
    try:
        with open(bin_list_file, 'r') as f:
            for line in f:
                bin_file = line.strip()
                if bin_file and os.path.exists(bin_file):
                    bin_files.append(bin_file)
                elif bin_file:
                    logging.warning(f"Bin file not found: {bin_file}")
        
        logging.info(f"Found {len(bin_files)} valid bin files in {bin_list_file}")
        return bin_files
        
    except Exception as e:
        logging.error(f"Error reading bin list file {bin_list_file}: {e}")
        return []

def extract_contigs_from_bin(bin_file):
    """Extract contig IDs from a bin FASTA file"""
    contigs = []
    total_length = 0
    
    try:
        with open(bin_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # Clean contig ID (remove description after first space)
                contig_id = record.id.split()[0]
                contigs.append(contig_id)
                total_length += len(record.seq)
        
        logging.debug(f"Extracted {len(contigs)} contigs from {bin_file}, total length: {total_length:,} bp")
        return contigs, total_length
        
    except Exception as e:
        logging.error(f"Error reading bin file {bin_file}: {e}")
        return [], 0

def validate_bin_file(bin_file, min_contigs=1, min_total_length=1000):
    """Validate if bin file meets minimum requirements"""
    if not os.path.exists(bin_file):
        return False, "File does not exist"
    
    if os.path.getsize(bin_file) == 0:
        return False, "File is empty"
    
    try:
        contigs, total_length = extract_contigs_from_bin(bin_file)
        
        if len(contigs) < min_contigs:
            return False, f"Too few contigs ({len(contigs)} < {min_contigs})"
        
        if total_length < min_total_length:
            return False, f"Total length too small ({total_length} < {min_total_length})"
        
        return True, f"Valid bin: {len(contigs)} contigs, {total_length:,} bp"
        
    except Exception as e:
        return False, f"Error validating file: {e}"

def create_scaffolds2bin_file(bin_files, tool_name, output_dir, min_contigs=5, min_bin_size=50000):
    """Create scaffolds2bin file for DAS Tool"""
    if not bin_files:
        logging.warning(f"No bin files found for {tool_name}")
        return None
    
    output_file = Path(output_dir) / f"{tool_name}.scaffolds2bin.tsv"
    
    valid_bins = 0
    total_contigs = 0
    bin_stats = []
    
    try:
        with open(output_file, 'w') as f:
            for i, bin_file in enumerate(bin_files):
                # Validate bin file
                is_valid, validation_msg = validate_bin_file(
                    bin_file, min_contigs=min_contigs, min_total_length=min_bin_size
                )
                
                if not is_valid:
                    logging.warning(f"Skipping invalid bin {bin_file}: {validation_msg}")
                    continue
                
                bin_name = f"{tool_name}.{i+1}"
                contigs, total_length = extract_contigs_from_bin(bin_file)
                
                if contigs:
                    contig_count = 0
                    for contig in contigs:
                        # Write contig-bin mapping
                        f.write(f"{contig}\t{bin_name}\n")
                        contig_count += 1
                        total_contigs += 1
                    
                    valid_bins += 1
                    bin_stats.append({
                        'bin_file': bin_file,
                        'bin_name': bin_name,
                        'contigs': contig_count,
                        'total_length': total_length
                    })
                    
                    logging.info(f"Processed {tool_name} bin {i+1}: {contig_count} contigs, {total_length:,} bp")
        
        if valid_bins > 0:
            logging.info(f"Created {tool_name} scaffolds2bin file: {output_file}")
            logging.info(f"  Valid bins: {valid_bins}")
            logging.info(f"  Total contigs: {total_contigs:,}")
            
            # Save bin statistics
            stats_file = Path(output_dir) / f"{tool_name}_bin_stats.json"
            with open(stats_file, 'w') as f:
                json.dump({
                    'tool_name': tool_name,
                    'valid_bins': valid_bins,
                    'total_contigs': total_contigs,
                    'bin_details': bin_stats
                }, f, indent=2)
            
            return str(output_file)
        else:
            logging.error(f"No valid bins found for {tool_name}")
            return None
        
    except Exception as e:
        logging.error(f"Error creating scaffolds2bin file for {tool_name}: {e}")
        return None

def validate_scaffolds2bin_file(file_path):
    """Validate scaffolds2bin file format"""
    if not os.path.exists(file_path):
        return False
    
    try:
        line_count = 0
        unique_contigs = set()
        unique_bins = set()
        
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line:
                    parts = line.split('\t')
                    if len(parts) != 2:
                        logging.error(f"Invalid line format in {file_path} at line {line_num}: {line}")
                        return False
                    
                    contig_id, bin_id = parts
                    unique_contigs.add(contig_id)
                    unique_bins.add(bin_id)
                    line_count += 1
        
        logging.info(f"Validated {file_path}:")
        logging.info(f"  {line_count:,} contig-bin mappings")
        logging.info(f"  {len(unique_contigs):,} unique contigs")
        logging.info(f"  {len(unique_bins):,} unique bins")
        
        return line_count > 0
        
    except Exception as e:
        logging.error(f"Error validating {file_path}: {e}")
        return False

def validate_contig_ids(created_files, assembly_file):
    """Validate that contig IDs in bins exist in assembly"""
    if not assembly_file or not os.path.exists(assembly_file):
        logging.warning("No assembly file provided for validation")
        return True
    
    logging.info("Validating contig IDs against assembly")
    
    # Load assembly contig IDs
    assembly_contigs = set()
    try:
        with open(assembly_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                assembly_contigs.add(record.id.split()[0])
        
        logging.info(f"Loaded {len(assembly_contigs):,} contig IDs from assembly")
        
    except Exception as e:
        logging.error(f"Error reading assembly file: {e}")
        return False
    
    # Check each scaffolds2bin file
    validation_results = {}
    
    for tool_name, file_path in created_files.items():
        if not file_path or not os.path.exists(file_path):
            continue
        
        missing_contigs = set()
        total_contigs = 0
        
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if line.strip():
                        contig_id = line.split('\t')[0]
                        total_contigs += 1
                        
                        if contig_id not in assembly_contigs:
                            missing_contigs.add(contig_id)
            
            validation_results[tool_name] = {
                'total_contigs': total_contigs,
                'missing_contigs': len(missing_contigs),
                'valid_percentage': ((total_contigs - len(missing_contigs)) / total_contigs * 100) if total_contigs > 0 else 0
            }
            
            if missing_contigs:
                logging.warning(f"{tool_name}: {len(missing_contigs)} contig IDs not found in assembly")
                if len(missing_contigs) <= 10:  # Show first 10 missing IDs
                    logging.warning(f"  Missing IDs: {', '.join(list(missing_contigs)[:10])}")
            else:
                logging.info(f"{tool_name}: All contig IDs validated successfully")
                
        except Exception as e:
            logging.error(f"Error validating {tool_name}: {e}")
            validation_results[tool_name] = {'error': str(e)}
    
    return validation_results

def create_summary_file(output_dir, created_files):
    """Create a summary of created DAS Tool input files"""
    summary_file = Path(output_dir) / "dastool_input_summary.txt"
    
    try:
        with open(summary_file, 'w') as f:
            f.write("DAS Tool Input Preparation Summary\n")
            f.write("=" * 40 + "\n\n")
            
            total_files = 0
            total_mappings = 0
            
            for tool_name, file_path in created_files.items():
                if file_path and os.path.exists(file_path):
                    total_files += 1
                    
                    # Count mappings
                    mappings = 0
                    unique_bins = set()
                    try:
                        with open(file_path, 'r') as input_f:
                            for line in input_f:
                                if line.strip():
                                    parts = line.strip().split('\t')
                                    if len(parts) == 2:
                                        mappings += 1
                                        unique_bins.add(parts[1])
                    except Exception as e:
                        mappings = 0
                        logging.warning(f"Error counting mappings for {tool_name}: {e}")
                    
                    total_mappings += mappings
                    
                    f.write(f"{tool_name}:\n")
                    f.write(f"  File: {file_path}\n")
                    f.write(f"  Contig-bin mappings: {mappings:,}\n")
                    f.write(f"  Unique bins: {len(unique_bins)}\n")
                    f.write(f"  Status: ✓ Ready for DAS Tool\n\n")
                else:
                    f.write(f"{tool_name}:\n")
                    f.write(f"  Status: ✗ No valid bins found\n\n")
            
            f.write(f"Summary:\n")
            f.write(f"  Tools with valid input: {total_files}\n")
            f.write(f"  Total contig-bin mappings: {total_mappings:,}\n")
            
            if total_files >= 2:
                f.write(f"  Ready for DAS Tool refinement: ✓\n")
                
                # Add DAS Tool command
                f.write(f"\nSuggested DAS Tool command:\n")
                valid_files = [f for f in created_files.values() if f is not None]
                tool_names = [tool for tool, file_path in created_files.items() if file_path is not None]
                
                f.write(f"DAS_Tool \\\n")
                f.write(f"  -i {','.join(valid_files)} \\\n")
                f.write(f"  -l {','.join(tool_names)} \\\n")
                f.write(f"  -c ASSEMBLY.fasta \\\n")
                f.write(f"  -o DAS_Tool_output \\\n")
                f.write(f"  --write_bins \\\n")
                f.write(f"  --threads 8\n")
            else:
                f.write(f"  Ready for DAS Tool refinement: ✗ (need at least 2 tools)\n")
        
        logging.info(f"Summary saved to {summary_file}")
        
    except Exception as e:
        logging.error(f"Error creating summary file: {e}")

def main():
    """Main function"""
    args = parse_arguments()
    setup_logging(args.log_level)
    
    logging.info("Starting DAS Tool input preparation")
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each binning tool
    tools = {
        'metabat2': args.metabat2,
        'maxbin2': args.maxbin2,
        'concoct': args.concoct
    }
    
    created_files = {}
    
    for tool_name, bin_list_file in tools.items():
        if bin_list_file:
            logging.info(f"Processing {tool_name} bins")
            
            # Read bin files
            bin_files = read_bin_files(bin_list_file)
            
            if not bin_files:
                logging.warning(f"No valid bin files found for {tool_name}")
                created_files[tool_name] = None
                continue
            
            # Create scaffolds2bin file with validation
            scaffolds2bin_file = create_scaffolds2bin_file(
                bin_files, tool_name, output_dir, 
                min_contigs=args.min_contigs,
                min_bin_size=args.min_bin_size
            )
            
            if scaffolds2bin_file and validate_scaffolds2bin_file(scaffolds2bin_file):
                created_files[tool_name] = scaffolds2bin_file
            else:
                logging.warning(f"Failed to create valid scaffolds2bin file for {tool_name}")
                if not args.ignore_errors:
                    created_files[tool_name] = None
                else:
                    logging.info(f"Ignoring errors for {tool_name} (--ignore-errors flag set)")
                    created_files[tool_name] = scaffolds2bin_file
        else:
            logging.info(f"No bin list file provided for {tool_name}")
            created_files[tool_name] = None
    
    # Validate contig IDs if assembly provided
    if args.validate_contigs:
        validation_results = validate_contig_ids(created_files, args.validate_contigs)
        
        # Save validation results
        validation_file = output_dir / "contig_validation_results.json"
        with open(validation_file, 'w') as f:
            json.dump(validation_results, f, indent=2)
        
        logging.info(f"Contig validation results saved to {validation_file}")
    
    # Create summary
    create_summary_file(output_dir, created_files)
    
    # Check if we have enough tools for DAS Tool
    valid_tools = sum(1 for file_path in created_files.values() if file_path is not None)
    
    if valid_tools >= 2:
        logging.info(f"✓ Successfully prepared input for {valid_tools} tools")
        logging.info("Ready for DAS Tool bin refinement")
        
        # Print DAS Tool command hint
        valid_files = [f for f in created_files.values() if f is not None]
        tool_names = [tool for tool, file_path in created_files.items() if file_path is not None]
        
        print(f"\n🔧 DAS Tool Input Preparation Complete")
        print("=" * 40)
        print(f"Valid input files: {valid_tools}")
        print(f"Tools: {', '.join(tool_names)}")
        print("\nSuggested DAS Tool command:")
        print(f"DAS_Tool \\")
        print(f"  -i {','.join(valid_files)} \\")
        print(f"  -l {','.join(tool_names)} \\")
        print(f"  -c ASSEMBLY.fasta \\")
        print(f"  -o {args.prefix} \\")
        print(f"  --write_bins \\")
        print(f"  --threads 8")
        
    else:
        logging.warning(f"Only {valid_tools} valid tool(s) found")
        logging.warning("DAS Tool requires at least 2 binning tools for refinement")
        
        if valid_tools == 1:
            logging.info("Consider using the single binning result directly, or run additional binning tools")
        elif valid_tools == 0:
            logging.error("No valid binning results found. Check your input files and try again")
        
        if not args.ignore_errors:
            sys.exit(1)
    
    logging.info("DAS Tool input preparation completed successfully")

if __name__ == "__main__":
    main()

