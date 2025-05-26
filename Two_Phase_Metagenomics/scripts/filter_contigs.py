#!/usr/bin/env python3
"""
🔧 Contig Filtering Utility
Part of the Metagenomics Pipeline - Phase 1

Filters assembly contigs based on length, coverage, GC content, complexity, and quality criteria.
"""

import argparse
import os
import sys
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import re
import numpy as np

def setup_logging(log_level="INFO"):
    """Setup logging configuration"""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Filter assembly contigs by length, coverage, and quality"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input assembly FASTA file"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output filtered FASTA file"
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=1000,
        help="Minimum contig length (default: 1000)"
    )
    parser.add_argument(
        "--max-length",
        type=int,
        help="Maximum contig length (optional)"
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=3.0,
        help="Minimum coverage depth (default: 3.0)"
    )
    parser.add_argument(
        "--max-coverage",
        type=float,
        help="Maximum coverage depth (optional)"
    )
    parser.add_argument(
        "--min-gc",
        type=float,
        default=20.0,
        help="Minimum GC content percentage (default: 20.0)"
    )
    parser.add_argument(
        "--max-gc",
        type=float,
        default=80.0,
        help="Maximum GC content percentage (default: 80.0)"
    )
    parser.add_argument(
        "--max-n-content",
        type=float,
        default=5.0,
        help="Maximum N content percentage (default: 5.0)"
    )
    parser.add_argument(
        "--exclude-patterns",
        nargs="*",
        default=["plasmid", "phage", "virus"],
        help="Exclude contigs with these patterns in headers"
    )
    parser.add_argument(
        "--min-complexity",
        type=float,
        default=1.0,
        help="Minimum sequence complexity (Shannon entropy, default: 1.0)"
    )
    parser.add_argument(
        "--max-run-fraction",
        type=float,
        default=0.8,
        help="Maximum fraction of sequence as nucleotide runs (default: 0.8)"
    )
    parser.add_argument(
        "--require-orf",
        action="store_true",
        help="Require potential open reading frames"
    )
    parser.add_argument(
        "--min-orf-length",
        type=int,
        default=300,
        help="Minimum ORF length in bp (default: 300)"
    )
    parser.add_argument(
        "--remove-duplicates",
        action="store_true",
        help="Remove duplicate sequences"
    )
    parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=0.95,
        help="Similarity threshold for duplicate removal (default: 0.95)"
    )
    parser.add_argument(
        "--stats-output",
        help="Output file for filtering statistics"
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )
    
    return parser.parse_args()

def extract_coverage_from_header(header):
    """Extract coverage information from contig header"""
    coverage = None
    
    # Common SPAdes coverage patterns
    patterns = [
        r'cov_(\d+\.?\d*)',           # cov_X.X
        r'coverage_(\d+\.?\d*)',      # coverage_X.X
        r'_(\d+\.?\d*)x',             # _X.Xx
        r'cov=(\d+\.?\d*)',           # cov=X.X
        r'depth=(\d+\.?\d*)',         # depth=X.X
        r'multi=(\d+\.?\d*)',         # multi=X.X (SPAdes multiplicity)
        r'len_\d+_cov_(\d+\.?\d*)',   # SPAdes format: len_1234_cov_5.67
        r'_length_\d+_cov_(\d+\.?\d*)', # Alternative SPAdes format
        r'NODE_\d+_length_\d+_cov_(\d+\.?\d*)', # Full SPAdes format
    ]
    
    for pattern in patterns:
        match = re.search(pattern, header, re.IGNORECASE)
        if match:
            try:
                coverage = float(match.group(1))
                break
            except ValueError:
                continue
    
    return coverage

def extract_length_from_header(header):
    """Extract sequence length from contig header if available"""
    length = None
    
    # Common length patterns
    patterns = [
        r'length_(\d+)',              # length_1234
        r'len_(\d+)',                 # len_1234
        r'size_(\d+)',                # size_1234
        r'NODE_\d+_length_(\d+)',     # SPAdes format
    ]
    
    for pattern in patterns:
        match = re.search(pattern, header, re.IGNORECASE)
        if match:
            try:
                length = int(match.group(1))
                break
            except ValueError:
                continue
    
    return length

def calculate_sequence_complexity(sequence):
    """Calculate sequence complexity metrics"""
    sequence = sequence.upper()
    length = len(sequence)
    
    if length == 0:
        return 0.0, 0.0, 0.0
    
    # Calculate base composition
    base_counts = {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'G': sequence.count('G'),
        'C': sequence.count('C'),
        'N': sequence.count('N')
    }
    
    # Shannon entropy (complexity measure)
    total_bases = sum(base_counts.values())
    entropy = 0.0
    if total_bases > 0:
        for count in base_counts.values():
            if count > 0:
                p = count / total_bases
                entropy -= p * np.log2(p)
    
    # Low complexity regions (simple repeats)
    # Count runs of the same nucleotide
    max_run = 0
    current_run = 1
    for i in range(1, length):
        if sequence[i] == sequence[i-1]:
            current_run += 1
        else:
            max_run = max(max_run, current_run)
            current_run = 1
    max_run = max(max_run, current_run)
    
    # Dinucleotide complexity
    dinuc_counts = {}
    for i in range(length - 1):
        dinuc = sequence[i:i+2]
        dinuc_counts[dinuc] = dinuc_counts.get(dinuc, 0) + 1
    
    dinuc_entropy = 0.0
    total_dinucs = len(dinuc_counts)
    if total_dinucs > 0:
        for count in dinuc_counts.values():
            p = count / total_dinucs
            dinuc_entropy -= p * np.log2(p)
    
    return entropy, max_run / length, dinuc_entropy

def has_potential_orf(sequence, min_orf_length=300):
    """Check if sequence has potential open reading frames"""
    sequence = sequence.upper()
    
    # Start and stop codons
    start_codons = ['ATG', 'GTG', 'TTG']
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    # Check all 6 reading frames
    for frame in range(3):
        # Forward strand
        for start_pos in range(frame, len(sequence) - 2, 3):
            codon = sequence[start_pos:start_pos + 3]
            if len(codon) == 3 and codon in start_codons:
                # Look for stop codon
                for stop_pos in range(start_pos + 3, len(sequence) - 2, 3):
                    stop_codon = sequence[stop_pos:stop_pos + 3]
                    if len(stop_codon) == 3 and stop_codon in stop_codons:
                        orf_length = stop_pos - start_pos
                        if orf_length >= min_orf_length:
                            return True
                        break
        
        # Reverse strand
        rev_sequence = reverse_complement(sequence)
        for start_pos in range(frame, len(rev_sequence) - 2, 3):
            codon = rev_sequence[start_pos:start_pos + 3]
            if len(codon) == 3 and codon in start_codons:
                # Look for stop codon
                for stop_pos in range(start_pos + 3, len(rev_sequence) - 2, 3):
                    stop_codon = rev_sequence[stop_pos:stop_pos + 3]
                    if len(stop_codon) == 3 and stop_codon in stop_codons:
                        orf_length = stop_pos - start_pos
                        if orf_length >= min_orf_length:
                            return True
                        break
    
    return False

def reverse_complement(sequence):
    """Get reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))

def calculate_sequence_similarity(seq1, seq2):
    """Calculate sequence similarity using simple alignment"""
    if len(seq1) != len(seq2):
        return 0.0
    
    matches = sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a == b)
    return matches / len(seq1)

def passes_filters(record, args):
    """Check if a contig passes all filtering criteria"""
    sequence = str(record.seq)
    length = len(sequence)
    header = record.description
    
    # Length filters
    if length < args.min_length:
        return False, "length_too_short"
    
    if args.max_length and length > args.max_length:
        return False, "length_too_long"
    
    # Coverage filters
    coverage = extract_coverage_from_header(header)
    if coverage is not None:
        if coverage < args.min_coverage:
            return False, "coverage_too_low"
        
        if args.max_coverage and coverage > args.max_coverage:
            return False, "coverage_too_high"
    
    # GC content filters
    gc_content = gc_fraction(record.seq) * 100
    if gc_content < args.min_gc:
        return False, "gc_too_low"
    
    if gc_content > args.max_gc:
        return False, "gc_too_high"
    
    # N content filter
    n_count = sequence.upper().count('N')
    n_content = (n_count / length) * 100
    if n_content > args.max_n_content:
        return False, "n_content_too_high"
    
    # Sequence complexity filters
    if hasattr(args, 'min_complexity') and args.min_complexity:
        entropy, max_run_fraction, dinuc_entropy = calculate_sequence_complexity(sequence)
        
        # Filter out low complexity sequences
        if entropy < args.min_complexity:
            return False, "low_complexity"
        
        # Filter out sequences with long runs of same nucleotide
        if hasattr(args, 'max_run_fraction') and args.max_run_fraction:
            if max_run_fraction > args.max_run_fraction:
                return False, "long_nucleotide_runs"
    
    # Pattern exclusion filters
    if args.exclude_patterns:
        header_lower = header.lower()
        for pattern in args.exclude_patterns:
            if pattern.lower() in header_lower:
                return False, f"excluded_pattern_{pattern}"
    
    # Additional quality filters
    if hasattr(args, 'require_orf') and args.require_orf:
        if not has_potential_orf(sequence, args.min_orf_length):
            return False, "no_potential_orf"
    
    return True, "passed"

def filter_contigs(input_file, output_file, args):
    """Filter contigs and write results"""
    logging.info(f"Filtering contigs from {input_file}")
    
    # Statistics tracking
    stats = {
        "total_contigs": 0,
        "passed_contigs": 0,
        "total_length": 0,
        "passed_length": 0,
        "filter_reasons": {},
        "duplicates_removed": 0
    }
    
    # First pass: collect all sequences and filter by criteria
    filtered_records = []
    seen_sequences = set()
    sequence_hashes = {}
    
    try:
        with open(input_file, 'r') as in_handle:
            for record in SeqIO.parse(in_handle, "fasta"):
                stats["total_contigs"] += 1
                stats["total_length"] += len(record.seq)
                
                passed, reason = passes_filters(record, args)
                
                if passed:
                    # Check for duplicates if requested
                    if args.remove_duplicates:
                        seq_hash = hash(str(record.seq).upper())
                        if seq_hash in sequence_hashes:
                            # Check similarity more precisely
                            existing_seq = sequence_hashes[seq_hash]
                            similarity = calculate_sequence_similarity(
                                str(record.seq), str(existing_seq.seq)
                            )
                            
                            if similarity >= args.similarity_threshold:
                                stats["duplicates_removed"] += 1
                                if "duplicate_sequence" not in stats["filter_reasons"]:
                                    stats["filter_reasons"]["duplicate_sequence"] = 0
                                stats["filter_reasons"]["duplicate_sequence"] += 1
                                continue
                        
                        sequence_hashes[seq_hash] = record
                    
                    filtered_records.append(record)
                    stats["passed_contigs"] += 1
                    stats["passed_length"] += len(record.seq)
                else:
                    if reason not in stats["filter_reasons"]:
                        stats["filter_reasons"][reason] = 0
                    stats["filter_reasons"][reason] += 1
        
        # Second pass: write filtered sequences
        with open(output_file, 'w') as out_handle:
            for record in filtered_records:
                SeqIO.write(record, out_handle, "fasta")
        
        # Log statistics
        logging.info(f"Filtering complete:")
        logging.info(f"  Input contigs: {stats['total_contigs']:,}")
        logging.info(f"  Passed contigs: {stats['passed_contigs']:,} ({stats['passed_contigs']/stats['total_contigs']*100:.1f}%)")
        logging.info(f"  Input length: {stats['total_length']:,} bp")
        logging.info(f"  Passed length: {stats['passed_length']:,} bp ({stats['passed_length']/stats['total_length']*100:.1f}%)")
        
        if args.remove_duplicates:
            logging.info(f"  Duplicates removed: {stats['duplicates_removed']:,}")
        
        if stats["filter_reasons"]:
            logging.info("Filter reasons:")
            for reason, count in sorted(stats["filter_reasons"].items()):
                logging.info(f"    {reason}: {count:,} contigs")
        
        # Save statistics if requested
        if args.stats_output:
            save_statistics(stats, args.stats_output)
        
        return stats
        
    except Exception as e:
        logging.error(f"Error filtering contigs: {e}")
        sys.exit(1)

def save_statistics(stats, stats_file):
    """Save filtering statistics to file"""
    try:
        with open(stats_file, 'w') as f:
            f.write("Contig Filtering Statistics\n")
            f.write("=" * 30 + "\n\n")
            
            f.write(f"Total input contigs: {stats['total_contigs']:,}\n")
            f.write(f"Passed contigs: {stats['passed_contigs']:,}\n")
            f.write(f"Retention rate: {stats['passed_contigs']/stats['total_contigs']*100:.1f}%\n\n")
            
            f.write(f"Total input length: {stats['total_length']:,} bp\n")
            f.write(f"Passed length: {stats['passed_length']:,} bp\n")
            f.write(f"Length retention: {stats['passed_length']/stats['total_length']*100:.1f}%\n\n")
            
            if stats.get("duplicates_removed", 0) > 0:
                f.write(f"Duplicates removed: {stats['duplicates_removed']:,}\n\n")
            
            if stats["filter_reasons"]:
                f.write("Filtering reasons:\n")
                for reason, count in sorted(stats["filter_reasons"].items()):
                    percentage = count / stats['total_contigs'] * 100
                    f.write(f"  {reason}: {count:,} contigs ({percentage:.1f}%)\n")
        
        logging.info(f"Statistics saved to {stats_file}")
        
    except Exception as e:
        logging.error(f"Error saving statistics: {e}")

def main():
    """Main function"""
    args = parse_arguments()
    setup_logging(args.log_level)
    
    logging.info("Starting contig filtering")
    
    # Validate inputs
    if not os.path.exists(args.input):
        logging.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Validate parameters
    if args.min_length <= 0:
        logging.error("Minimum length must be positive")
        sys.exit(1)
    
    if args.max_length and args.max_length <= args.min_length:
        logging.error("Maximum length must be greater than minimum length")
        sys.exit(1)
    
    if args.min_coverage < 0:
        logging.error("Minimum coverage must be non-negative")
        sys.exit(1)
    
    if args.max_coverage and args.max_coverage <= args.min_coverage:
        logging.error("Maximum coverage must be greater than minimum coverage")
        sys.exit(1)
    
    if not (0 <= args.min_gc <= 100) or not (0 <= args.max_gc <= 100):
        logging.error("GC content percentages must be between 0 and 100")
        sys.exit(1)
    
    if args.min_gc >= args.max_gc:
        logging.error("Minimum GC must be less than maximum GC")
        sys.exit(1)
    
    # Filter contigs
    stats = filter_contigs(args.input, args.output, args)
    
    # Print summary
    print(f"\n🔧 Contig Filtering Summary")
    print("=" * 30)
    print(f"Input contigs: {stats['total_contigs']:,}")
    print(f"Passed contigs: {stats['passed_contigs']:,}")
    print(f"Retention rate: {stats['passed_contigs']/stats['total_contigs']*100:.1f}%")
    print(f"Output file: {args.output}")
    
    if args.remove_duplicates and stats.get('duplicates_removed', 0) > 0:
        print(f"Duplicates removed: {stats['duplicates_removed']:,}")
    
    logging.info("Contig filtering completed successfully")

if __name__ == "__main__":
    main()

