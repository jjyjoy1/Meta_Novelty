#!/usr/bin/env python3
"""
🧬 Assembly Statistics Calculator
Part of the Metagenomics Pipeline - Phase 1

Calculates comprehensive assembly statistics including N50, coverage depth,
GC content, and quality metrics for metagenomic assemblies.
"""

import argparse
import os
import sys
import logging
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
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
        description="Calculate comprehensive assembly statistics"
    )
    parser.add_argument(
        "--assembly", "-a",
        required=True,
        help="Path to assembly FASTA file"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "--coverage-file", "-c",
        help="Coverage file from mapping (optional)"
    )
    parser.add_argument(
        "--sample-name", "-s",
        default="sample",
        help="Sample name for output files"
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=1000,
        help="Minimum contig length to include (default: 1000)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use"
    )
    parser.add_argument(
        "--plots",
        action="store_true",
        help="Generate visualization plots"
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )
    
    return parser.parse_args()

class AssemblyAnalyzer:
    """Main class for assembly analysis"""
    
    def __init__(self, assembly_file, min_length=1000):
        self.assembly_file = assembly_file
        self.min_length = min_length
        self.contigs = {}
        self.stats = {}
        
    def load_assembly(self):
        """Load and parse assembly file"""
        logging.info(f"Loading assembly from {self.assembly_file}")
        
        try:
            with open(self.assembly_file, 'r') as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if len(record.seq) >= self.min_length:
                        self.contigs[record.id] = {
                            'sequence': str(record.seq),
                            'length': len(record.seq),
                            'gc_content': gc_fraction(record.seq) * 100,
                            'n_content': (str(record.seq).upper().count('N') / len(record.seq)) * 100
                        }
            
            logging.info(f"Loaded {len(self.contigs)} contigs >= {self.min_length} bp")
            
        except Exception as e:
            logging.error(f"Error loading assembly: {e}")
            sys.exit(1)
    
    def calculate_basic_stats(self):
        """Calculate basic assembly statistics"""
        logging.info("Calculating basic assembly statistics")
        
        lengths = [contig['length'] for contig in self.contigs.values()]
        gc_contents = [contig['gc_content'] for contig in self.contigs.values()]
        n_contents = [contig['n_content'] for contig in self.contigs.values()]
        
        # Sort lengths in descending order for N50 calculation
        lengths_sorted = sorted(lengths, reverse=True)
        total_length = sum(lengths)
        
        # Calculate N50, N75, N90
        def calculate_nx(x, lengths_sorted, total_length):
            target = total_length * (x / 100)
            cumulative = 0
            for length in lengths_sorted:
                cumulative += length
                if cumulative >= target:
                    return length
            return 0
        
        n50 = calculate_nx(50, lengths_sorted, total_length)
        n75 = calculate_nx(75, lengths_sorted, total_length)
        n90 = calculate_nx(90, lengths_sorted, total_length)
        
        # Calculate L50 (number of contigs needed to reach N50)
        l50 = 0
        cumulative = 0
        for length in lengths_sorted:
            l50 += 1
            cumulative += length
            if cumulative >= total_length * 0.5:
                break
        
        self.stats.update({
            'total_contigs': len(self.contigs),
            'total_length': total_length,
            'max_contig_length': max(lengths) if lengths else 0,
            'min_contig_length': min(lengths) if lengths else 0,
            'mean_contig_length': np.mean(lengths) if lengths else 0,
            'median_contig_length': np.median(lengths) if lengths else 0,
            'n50': n50,
            'n75': n75,
            'n90': n90,
            'l50': l50,
            'mean_gc_content': np.mean(gc_contents) if gc_contents else 0,
            'median_gc_content': np.median(gc_contents) if gc_contents else 0,
            'std_gc_content': np.std(gc_contents) if gc_contents else 0,
            'mean_n_content': np.mean(n_contents) if n_contents else 0,
            'contigs_over_1kb': sum(1 for l in lengths if l >= 1000),
            'contigs_over_10kb': sum(1 for l in lengths if l >= 10000),
            'contigs_over_100kb': sum(1 for l in lengths if l >= 100000),
            'contigs_over_1mb': sum(1 for l in lengths if l >= 1000000)
        })
    
    def load_coverage_data(self, coverage_file):
        """Load coverage data if available"""
        if not coverage_file or not os.path.exists(coverage_file):
            logging.warning("Coverage file not provided or doesn't exist")
            return
        
        logging.info(f"Loading coverage data from {coverage_file}")
        
        try:
            # Try to read coverage file (assuming samtools depth format)
            coverage_data = {}
            with open(coverage_file, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            contig_id = parts[0]
                            coverage = float(parts[2])
                            if contig_id not in coverage_data:
                                coverage_data[contig_id] = []
                            coverage_data[contig_id].append(coverage)
            
            # Calculate mean coverage per contig
            for contig_id in self.contigs:
                if contig_id in coverage_data:
                    self.contigs[contig_id]['coverage'] = np.mean(coverage_data[contig_id])
                else:
                    self.contigs[contig_id]['coverage'] = 0.0
            
            # Calculate coverage statistics
            coverages = [contig['coverage'] for contig in self.contigs.values()]
            self.stats.update({
                'mean_coverage': np.mean(coverages) if coverages else 0,
                'median_coverage': np.median(coverages) if coverages else 0,
                'std_coverage': np.std(coverages) if coverages else 0,
                'max_coverage': max(coverages) if coverages else 0,
                'min_coverage': min(coverages) if coverages else 0
            })
            
        except Exception as e:
            logging.error(f"Error loading coverage data: {e}")
    
    def calculate_quality_metrics(self):
        """Calculate assembly quality metrics"""
        logging.info("Calculating quality metrics")
        
        lengths = [contig['length'] for contig in self.contigs.values()]
        
        if not lengths:
            return
        
        # Assembly contiguity metrics
        total_length = sum(lengths)
        
        # Effective genome size (sum of squared contig lengths / total length)
        effective_genome_size = sum(l**2 for l in lengths) / total_length
        
        # Assembly fragmentation
        # Lower values indicate more fragmented assemblies
        auN = sum(l**2 for l in lengths) / sum(lengths)  # Area under N curve
        
        self.stats.update({
            'effective_genome_size': effective_genome_size,
            'auN': auN,
            'assembly_fragmentation': len(lengths) / (total_length / 1000000),  # contigs per Mb
        })
        
        # GC content distribution analysis
        gc_contents = [contig['gc_content'] for contig in self.contigs.values()]
        if gc_contents:
            # Calculate GC skew and identify potential contamination
            gc_q25 = np.percentile(gc_contents, 25)
            gc_q75 = np.percentile(gc_contents, 75)
            gc_iqr = gc_q75 - gc_q25
            
            # Identify outlier contigs (potential contamination)
            outlier_threshold_low = gc_q25 - 1.5 * gc_iqr
            outlier_threshold_high = gc_q75 + 1.5 * gc_iqr
            
            outlier_contigs = [
                contig_id for contig_id, data in self.contigs.items()
                if data['gc_content'] < outlier_threshold_low or 
                   data['gc_content'] > outlier_threshold_high
            ]
            
            self.stats.update({
                'gc_q25': gc_q25,
                'gc_q75': gc_q75,
                'gc_iqr': gc_iqr,
                'potential_contamination_contigs': len(outlier_contigs),
                'contamination_percentage': (len(outlier_contigs) / len(self.contigs)) * 100
            })
    
    def generate_contig_table(self):
        """Generate detailed contig information table"""
        logging.info("Generating contig information table")
        
        contig_data = []
        for contig_id, data in self.contigs.items():
            row = {
                'contig_id': contig_id,
                'length': data['length'],
                'gc_content': data['gc_content'],
                'n_content': data['n_content']
            }
            
            if 'coverage' in data:
                row['coverage'] = data['coverage']
            
            contig_data.append(row)
        
        return pd.DataFrame(contig_data)

def create_visualization_plots(analyzer, output_dir, sample_name):
    """Create visualization plots for assembly statistics"""
    logging.info("Creating visualization plots")
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create plots directory
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Extract data for plotting
    lengths = [contig['length'] for contig in analyzer.contigs.values()]
    gc_contents = [contig['gc_content'] for contig in analyzer.contigs.values()]
    
    # 1. Contig length distribution
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Length histogram
    ax1.hist(lengths, bins=50, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Contig Length (bp)')
    ax1.set_ylabel('Frequency')
    ax1.set_title(f'{sample_name}: Contig Length Distribution')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    
    # Length cumulative plot
    lengths_sorted = sorted(lengths, reverse=True)
    cumulative_lengths = np.cumsum(lengths_sorted)
    cumulative_percentages = (cumulative_lengths / cumulative_lengths[-1]) * 100
    
    ax2.plot(range(1, len(lengths_sorted) + 1), cumulative_percentages)
    ax2.axhline(y=50, color='red', linestyle='--', label='N50')
    ax2.axhline(y=75, color='orange', linestyle='--', label='N75')
    ax2.axhline(y=90, color='green', linestyle='--', label='N90')
    ax2.set_xlabel('Contig Rank')
    ax2.set_ylabel('Cumulative Assembly (%)')
    ax2.set_title(f'{sample_name}: Cumulative Assembly Coverage')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # GC content distribution
    ax3.hist(gc_contents, bins=30, edgecolor='black', alpha=0.7)
    ax3.axvline(analyzer.stats['mean_gc_content'], color='red', linestyle='--', 
                label=f"Mean: {analyzer.stats['mean_gc_content']:.1f}%")
    ax3.set_xlabel('GC Content (%)')
    ax3.set_ylabel('Frequency')
    ax3.set_title(f'{sample_name}: GC Content Distribution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Length vs GC content scatter
    if 'coverage' in list(analyzer.contigs.values())[0]:
        coverages = [contig['coverage'] for contig in analyzer.contigs.values()]
        scatter = ax4.scatter(lengths, gc_contents, c=coverages, 
                            alpha=0.6, cmap='viridis', s=20)
        plt.colorbar(scatter, ax=ax4, label='Coverage')
    else:
        ax4.scatter(lengths, gc_contents, alpha=0.6, s=20)
    
    ax4.set_xlabel('Contig Length (bp)')
    ax4.set_ylabel('GC Content (%)')
    ax4.set_title(f'{sample_name}: Length vs GC Content')
    ax4.set_xscale('log')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(plots_dir / f"{sample_name}_assembly_overview.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Assembly quality metrics plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Contig size categories
    categories = ['≥1kb', '≥10kb', '≥100kb', '≥1Mb']
    counts = [
        analyzer.stats['contigs_over_1kb'],
        analyzer.stats['contigs_over_10kb'],
        analyzer.stats['contigs_over_100kb'],
        analyzer.stats['contigs_over_1mb']
    ]
    
    ax1.bar(categories, counts, color='skyblue', edgecolor='black')
    ax1.set_ylabel('Number of Contigs')
    ax1.set_title(f'{sample_name}: Contigs by Size Category')
    ax1.grid(True, alpha=0.3)
    
    # Add counts on top of bars
    for i, count in enumerate(counts):
        ax1.text(i, count + max(counts) * 0.01, str(count), 
                ha='center', va='bottom')
    
    # N-statistics
    n_stats = ['N50', 'N75', 'N90']
    n_values = [analyzer.stats['n50'], analyzer.stats['n75'], analyzer.stats['n90']]
    
    ax2.bar(n_stats, n_values, color='lightcoral', edgecolor='black')
    ax2.set_ylabel('Contig Length (bp)')
    ax2.set_title(f'{sample_name}: N-Statistics')
    ax2.grid(True, alpha=0.3)
    
    # Add values on top of bars
    for i, value in enumerate(n_values):
        ax2.text(i, value + max(n_values) * 0.01, f'{value:,}', 
                ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(plots_dir / f"{sample_name}_quality_metrics.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Plots saved to {plots_dir}")

def save_results(analyzer, contig_df, output_dir, sample_name):
    """Save analysis results to files"""
    logging.info("Saving results")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save summary statistics
    stats_file = output_path / f"{sample_name}_assembly_stats.json"
    with open(stats_file, 'w') as f:
        json.dump(analyzer.stats, f, indent=2)
    
    # Save detailed contig table
    contig_file = output_path / f"{sample_name}_contig_stats.tsv"
    contig_df.to_csv(contig_file, sep='\t', index=False)
    
    # Save summary table
    summary_data = {
        'Sample': [sample_name],
        'Total_Contigs': [analyzer.stats['total_contigs']],
        'Total_Length': [analyzer.stats['total_length']],
        'N50': [analyzer.stats['n50']],
        'Mean_GC': [analyzer.stats['mean_gc_content']],
        'Max_Contig': [analyzer.stats['max_contig_length']],
        'Contigs_1kb+': [analyzer.stats['contigs_over_1kb']],
        'Contigs_10kb+': [analyzer.stats['contigs_over_10kb']]
    }
    
    if 'mean_coverage' in analyzer.stats:
        summary_data['Mean_Coverage'] = [analyzer.stats['mean_coverage']]
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_path / f"{sample_name}_assembly_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    
    logging.info(f"Results saved to {output_path}")
    
    # Print summary to console
    print(f"\n🧬 Assembly Statistics Summary for {sample_name}")
    print("=" * 50)
    print(f"Total contigs: {analyzer.stats['total_contigs']:,}")
    print(f"Total length: {analyzer.stats['total_length']:,} bp")
    print(f"N50: {analyzer.stats['n50']:,} bp")
    print(f"L50: {analyzer.stats['l50']:,} contigs")
    print(f"Max contig: {analyzer.stats['max_contig_length']:,} bp")
    print(f"Mean GC content: {analyzer.stats['mean_gc_content']:.1f}%")
    print(f"Contigs ≥1kb: {analyzer.stats['contigs_over_1kb']:,}")
    print(f"Contigs ≥10kb: {analyzer.stats['contigs_over_10kb']:,}")
    
    if 'mean_coverage' in analyzer.stats:
        print(f"Mean coverage: {analyzer.stats['mean_coverage']:.1f}x")

def main():
    """Main function"""
    args = parse_arguments()
    setup_logging(args.log_level)
    
    logging.info("Starting assembly statistics analysis")
    
    # Initialize analyzer
    analyzer = AssemblyAnalyzer(args.assembly, args.min_length)
    
    # Load and analyze assembly
    analyzer.load_assembly()
    analyzer.calculate_basic_stats()
    analyzer.load_coverage_data(args.coverage_file)
    analyzer.calculate_quality_metrics()
    
    # Generate contig table
    contig_df = analyzer.generate_contig_table()
    
    # Save results
    save_results(analyzer, contig_df, args.output, args.sample_name)
    
    # Generate plots if requested
    if args.plots:
        create_visualization_plots(analyzer, args.output, args.sample_name)
    
    logging.info("Assembly statistics analysis completed successfully")

if __name__ == "__main__":
    main()


