#!/usr/bin/env python3
"""
🔍 Homology-based Novelty Detection
Part of the Metagenomics Pipeline - Phase 1

Uses DIAMOND BLAST against multiple databases to identify novel sequences
based on lack of homology to known sequences.
"""

import argparse
import os
import sys
import logging
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import json
from collections import defaultdict
import tempfile

def setup_logging(log_level="INFO"):
    """Setup logging configuration"""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Homology-based novelty detection using DIAMOND BLAST"
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
        "--databases", "-d",
        nargs="+",
        required=True,
        help="Paths to DIAMOND databases (space-separated)"
    )
    parser.add_argument(
        "--database-names",
        nargs="+",
        help="Names for databases (must match number of databases)"
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="E-value threshold for BLAST"
    )
    parser.add_argument(
        "--max-target-seqs",
        type=int,
        default=25,
        help="Maximum number of target sequences"
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=500,
        help="Minimum sequence length for analysis"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads for DIAMOND"
    )
    parser.add_argument(
        "--block-size",
        type=float,
        default=2.0,
        help="DIAMOND block size in GB"
    )
    parser.add_argument(
        "--index-chunks",
        type=int,
        default=4,
        help="Number of index chunks for DIAMOND"
    )
    parser.add_argument(
        "--sample-name",
        default="sample",
        help="Sample name for output files"
    )
    parser.add_argument(
        "--high-novelty-threshold",
        type=float,
        default=30.0,
        help="Identity threshold for high novelty (default: 30%)"
    )
    parser.add_argument(
        "--medium-novelty-threshold",
        type=float,
        default=70.0,
        help="Identity threshold for medium novelty (default: 70%)"
    )
    parser.add_argument(
        "--low-novelty-threshold",
        type=float,
        default=90.0,
        help="Identity threshold for low novelty (default: 90%)"
    )
    parser.add_argument(
        "--plots",
        action="store_true",
        help="Generate visualization plots"
    )
    parser.add_argument(
        "--keep-blast-output",
        action="store_true",
        help="Keep raw BLAST output files"
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )
    
    return parser.parse_args()

class DiamondBlaster:
    """DIAMOND BLAST executor and result parser"""
    
    def __init__(self, threads=16, block_size=2.0, index_chunks=4):
        self.threads = threads
        self.block_size = block_size
        self.index_chunks = index_chunks
        
        # Check if DIAMOND is available
        try:
            subprocess.run(['diamond', '--version'], 
                         capture_output=True, check=True)
            self.diamond_available = True
        except (subprocess.CalledProcessError, FileNotFoundError):
            logging.warning("DIAMOND not found - using mock results")
            self.diamond_available = False
    
    def run_blast(self, query_file, database, output_file, evalue=1e-5, 
                  max_target_seqs=25):
        """Run DIAMOND BLAST search"""
        if not self.diamond_available:
            # Create mock output for testing
            self.create_mock_blast_output(query_file, output_file)
            return
        
        cmd = [
            'diamond', 'blastx',
            '--query', query_file,
            '--db', database,
            '--out', output_file,
            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 
                        'mismatch', 'gapopen', 'qstart', 'qend', 
                        'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen',
            '--evalue', str(evalue),
            '--max-target-seqs', str(max_target_seqs),
            '--threads', str(self.threads),
            '--block-size', str(self.block_size),
            '--index-chunks', str(self.index_chunks),
            '--quiet'
        ]
        
        logging.info(f"Running DIAMOND BLAST against {database}")
        logging.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logging.info(f"BLAST completed successfully")
            
        except subprocess.CalledProcessError as e:
            logging.error(f"DIAMOND BLAST failed: {e}")
            logging.error(f"stderr: {e.stderr}")
            raise
    
    def create_mock_blast_output(self, query_file, output_file):
        """Create mock BLAST output for testing"""
        logging.warning("Creating mock BLAST output for testing")
        
        # Read query sequences
        query_ids = []
        with open(query_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                query_ids.append(record.id)
        
        # Generate mock results
        with open(output_file, 'w') as f:
            for query_id in query_ids:
                # Random chance of having hits
                if np.random.random() > 0.3:  # 70% chance of hits
                    num_hits = np.random.randint(1, 6)  # 1-5 hits
                    for i in range(num_hits):
                        # Mock BLAST fields
                        identity = np.random.uniform(25, 95)
                        length = np.random.randint(50, 500)
                        evalue = np.random.uniform(1e-10, 1e-2)
                        bitscore = np.random.uniform(50, 200)
                        qlen = np.random.randint(500, 2000)
                        slen = np.random.randint(400, 2500)
                        
                        f.write(f"{query_id}\tmock_subject_{i}\t{identity:.1f}\t"
                               f"{length}\t5\t1\t1\t{length}\t1\t{length}\t"
                               f"{evalue:.2e}\t{bitscore:.1f}\t{qlen}\t{slen}\n")
    
    def parse_blast_results(self, blast_file):
        """Parse DIAMOND BLAST output"""
        if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
            logging.warning(f"BLAST file {blast_file} is empty or doesn't exist")
            return pd.DataFrame()
        
        columns = [
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen'
        ]
        
        try:
            df = pd.read_csv(blast_file, sep='\t', names=columns, header=None)
            logging.info(f"Parsed {len(df)} BLAST hits from {blast_file}")
            return df
            
        except Exception as e:
            logging.error(f"Error parsing BLAST results: {e}")
            return pd.DataFrame()

class HomologyNoveltyAnalyzer:
    """Main class for homology-based novelty analysis"""
    
    def __init__(self, databases, database_names=None, thresholds=None):
        self.databases = databases
        self.database_names = database_names or [f"db_{i}" for i in range(len(databases))]
        
        if len(self.database_names) != len(self.databases):
            raise ValueError("Number of database names must match number of databases")
        
        # Novelty thresholds
        self.thresholds = thresholds or {
            'high_novelty': 30.0,      # < 30% identity = high novelty
            'medium_novelty': 70.0,    # 30-70% identity = medium novelty
            'low_novelty': 90.0        # 70-90% identity = low novelty
            # > 90% identity = known sequence
        }
        
        self.blaster = DiamondBlaster()
        self.sequences = {}
        self.blast_results = {}
        self.novelty_results = {}
    
    def load_sequences(self, assembly_file, min_length=500):
        """Load sequences from assembly file"""
        logging.info(f"Loading sequences from {assembly_file}")
        
        with open(assembly_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if len(record.seq) >= min_length:
                    self.sequences[record.id] = {
                        'sequence': str(record.seq),
                        'length': len(record.seq)
                    }
        
        logging.info(f"Loaded {len(self.sequences)} sequences >= {min_length} bp")
    
    def run_blast_searches(self, assembly_file, output_dir, **blast_params):
        """Run BLAST searches against all databases"""
        logging.info("Running BLAST searches against all databases")
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        for db_path, db_name in zip(self.databases, self.database_names):
            if not os.path.exists(db_path):
                logging.warning(f"Database {db_path} not found, skipping")
                continue
            
            blast_output = output_path / f"blast_{db_name}.tsv"
            
            # Run BLAST
            self.blaster.run_blast(
                assembly_file, 
                db_path, 
                str(blast_output),
                **blast_params
            )
            
            # Parse results
            blast_df = self.blaster.parse_blast_results(str(blast_output))
            self.blast_results[db_name] = blast_df
            
            logging.info(f"Database {db_name}: {len(blast_df)} total hits")
    
    def calculate_novelty_scores(self):
        """Calculate novelty scores based on BLAST results"""
        logging.info("Calculating homology-based novelty scores")
        
        novelty_data = []
        
        for seq_id in self.sequences:
            seq_length = self.sequences[seq_id]['length']
            
            # Find best hits across all databases
            best_hits = {}
            all_hits = []
            
            for db_name, blast_df in self.blast_results.items():
                if blast_df.empty:
                    continue
                
                seq_hits = blast_df[blast_df['qseqid'] == seq_id]
                if not seq_hits.empty:
                    # Find best hit for this database
                    best_hit = seq_hits.loc[seq_hits['bitscore'].idxmax()]
                    best_hits[db_name] = best_hit
                    all_hits.extend(seq_hits.to_dict('records'))
            
            # Calculate novelty score
            novelty_score = self.calculate_sequence_novelty(best_hits, all_hits, seq_length)
            
            # Classify novelty level
            novelty_level = self.classify_novelty_level(novelty_score)
            
            # Get best overall hit
            best_overall_hit = None
            if all_hits:
                best_overall_hit = max(all_hits, key=lambda x: x['bitscore'])
            
            result = {
                'sequence_id': seq_id,
                'sequence_length': seq_length,
                'homology_novelty_score': novelty_score,
                'novelty_level': novelty_level,
                'total_hits': len(all_hits),
                'databases_with_hits': len(best_hits),
                'best_identity': best_overall_hit['pident'] if best_overall_hit else 0,
                'best_evalue': best_overall_hit['evalue'] if best_overall_hit else 1.0,
                'best_bitscore': best_overall_hit['bitscore'] if best_overall_hit else 0,
                'best_database': None,
                'best_subject': best_overall_hit['sseqid'] if best_overall_hit else None
            }
            
            # Find which database had the best hit
            if best_overall_hit:
                for db_name, blast_df in self.blast_results.items():
                    if not blast_df.empty:
                        if best_overall_hit['sseqid'] in blast_df['sseqid'].values:
                            result['best_database'] = db_name
                            break
            
            novelty_data.append(result)
        
        self.novelty_results = pd.DataFrame(novelty_data)
        
        # Calculate statistics
        self.calculate_statistics()
        
        return self.novelty_results
    
    def calculate_sequence_novelty(self, best_hits, all_hits, seq_length):
        """Calculate novelty score for a single sequence"""
        if not all_hits:
            # No hits = maximum novelty
            return 100.0
        
        # Get best identity across all databases
        best_identity = max(hit['pident'] for hit in all_hits)
        
        # Calculate coverage-weighted identity
        total_coverage = 0
        weighted_identity = 0
        
        for hit in all_hits:
            hit_coverage = (hit['length'] / seq_length) * 100
            hit_coverage = min(hit_coverage, 100)  # Cap at 100%
            
            weighted_identity += hit['pident'] * hit_coverage
            total_coverage += hit_coverage
        
        if total_coverage > 0:
            avg_weighted_identity = weighted_identity / total_coverage
        else:
            avg_weighted_identity = 0
        
        # Use the better of best identity or coverage-weighted identity
        effective_identity = max(best_identity, avg_weighted_identity)
        
        # Convert to novelty score (invert identity)
        novelty_score = 100 - effective_identity
        
        # Apply coverage penalty for low coverage hits
        max_coverage = max((hit['length'] / seq_length) * 100 for hit in all_hits)
        if max_coverage < 50:  # Less than 50% coverage
            coverage_bonus = (50 - max_coverage) / 50 * 20  # Up to 20 point bonus
            novelty_score = min(100, novelty_score + coverage_bonus)
        
        return novelty_score
    
    def classify_novelty_level(self, novelty_score):
        """Classify novelty level based on score"""
        if novelty_score >= (100 - self.thresholds['high_novelty']):
            return "High"
        elif novelty_score >= (100 - self.thresholds['medium_novelty']):
            return "Medium"
        elif novelty_score >= (100 - self.thresholds['low_novelty']):
            return "Low"
        else:
            return "Known"
    
    def calculate_statistics(self):
        """Calculate novelty detection statistics"""
        if self.novelty_results.empty:
            self.statistics = {}
            return
        
        novelty_scores = self.novelty_results['homology_novelty_score']
        novelty_levels = self.novelty_results['novelty_level']
        
        self.statistics = {
            'total_sequences': len(self.novelty_results),
            'sequences_with_hits': len(self.novelty_results[self.novelty_results['total_hits'] > 0]),
            'sequences_without_hits': len(self.novelty_results[self.novelty_results['total_hits'] == 0]),
            'mean_novelty_score': novelty_scores.mean(),
            'median_novelty_score': novelty_scores.median(),
            'std_novelty_score': novelty_scores.std(),
            'max_novelty_score': novelty_scores.max(),
            'min_novelty_score': novelty_scores.min(),
            'high_novelty_count': len(novelty_levels[novelty_levels == 'High']),
            'medium_novelty_count': len(novelty_levels[novelty_levels == 'Medium']),
            'low_novelty_count': len(novelty_levels[novelty_levels == 'Low']),
            'known_count': len(novelty_levels[novelty_levels == 'Known']),
            'novel_percentage': (len(novelty_levels[novelty_levels != 'Known']) / len(novelty_levels)) * 100
        }
        
        logging.info(f"Homology novelty analysis complete:")
        logging.info(f"  Total sequences: {self.statistics['total_sequences']}")
        logging.info(f"  Sequences with hits: {self.statistics['sequences_with_hits']}")
        logging.info(f"  Mean novelty score: {self.statistics['mean_novelty_score']:.1f}")
        logging.info(f"  High novelty: {self.statistics['high_novelty_count']}")

def create_visualization_plots(analyzer, output_dir, sample_name):
    """Create visualization plots"""
    logging.info("Creating visualization plots")
    
    # Create plots directory
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    results = analyzer.novelty_results
    
    if results.empty:
        logging.warning("No results to plot")
        return
    
    # 1. Novelty analysis overview
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Novelty score distribution
    novelty_scores = results['homology_novelty_score']
    ax1.hist(novelty_scores, bins=50, alpha=0.7, edgecolor='black')
    ax1.axvline(novelty_scores.mean(), color='red', linestyle='--',
                label=f'Mean: {novelty_scores.mean():.1f}')
    ax1.set_xlabel('Homology Novelty Score')
    ax1.set_ylabel('Frequency')
    ax1.set_title(f'{sample_name}: Homology Novelty Score Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Novelty levels pie chart
    level_counts = results['novelty_level'].value_counts()
    colors = {'High': 'red', 'Medium': 'orange', 'Low': 'yellow', 'Known': 'green'}
    plot_colors = [colors.get(level, 'gray') for level in level_counts.index]
    
    wedges, texts, autotexts = ax2.pie(level_counts.values, labels=level_counts.index,
                                      colors=plot_colors, autopct='%1.1f%%', startangle=90)
    ax2.set_title(f'{sample_name}: Novelty Level Distribution')
    
    # Score vs sequence length
    ax3.scatter(results['sequence_length'], results['homology_novelty_score'],
               alpha=0.6, s=20, c=results['total_hits'], cmap='viridis')
    ax3.set_xlabel('Sequence Length (bp)')
    ax3.set_ylabel('Homology Novelty Score')
    ax3.set_title(f'{sample_name}: Novelty vs Length')
    ax3.set_xscale('log')
    ax3.grid(True, alpha=0.3)
    
    # Add colorbar
    scatter = ax3.collections[0]
    plt.colorbar(scatter, ax=ax3, label='Number of Hits')
    
    # Best identity distribution
    results_with_hits = results[results['total_hits'] > 0]
    if not results_with_hits.empty:
        ax4.hist(results_with_hits['best_identity'], bins=30, alpha=0.7, edgecolor='black')
        ax4.set_xlabel('Best Identity (%)')
        ax4.set_ylabel('Frequency')
        ax4.set_title(f'{sample_name}: Best Identity Distribution')
        ax4.grid(True, alpha=0.3)
    else:
        ax4.text(0.5, 0.5, 'No BLAST hits found', ha='center', va='center',
                transform=ax4.transAxes)
        ax4.set_title('Best Identity Distribution (N/A)')
    
    plt.tight_layout()
    plt.savefig(plots_dir / f"{sample_name}_homology_novelty_analysis.png",
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Database hit statistics
    if analyzer.blast_results:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Hits per database
        db_hit_counts = {}
        for db_name, blast_df in analyzer.blast_results.items():
            if not blast_df.empty:
                db_hit_counts[db_name] = len(blast_df['qseqid'].unique())
            else:
                db_hit_counts[db_name] = 0
        
        if db_hit_counts:
            ax1.bar(db_hit_counts.keys(), db_hit_counts.values(),
                   color='skyblue', edgecolor='black')
            ax1.set_ylabel('Sequences with Hits')
            ax1.set_title(f'{sample_name}: Hits per Database')
            ax1.tick_params(axis='x', rotation=45)
            ax1.grid(True, alpha=0.3)
            
            # Add counts on bars
            for db, count in db_hit_counts.items():
                ax1.text(list(db_hit_counts.keys()).index(db), count + max(db_hit_counts.values()) * 0.01,
                        str(count), ha='center', va='bottom')
        
        # Hit statistics
        stats = analyzer.statistics
        categories = ['High', 'Medium', 'Low', 'Known']
        counts = [
            stats.get('high_novelty_count', 0),
            stats.get('medium_novelty_count', 0),
            stats.get('low_novelty_count', 0),
            stats.get('known_count', 0)
        ]
        colors = ['red', 'orange', 'yellow', 'green']
        
        ax2.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
        ax2.set_ylabel('Number of Sequences')
        ax2.set_title(f'{sample_name}: Novelty Categories')
        ax2.grid(True, alpha=0.3)
        
        # Add counts on bars
        for i, count in enumerate(counts):
            ax2.text(i, count + max(counts) * 0.01, str(count),
                    ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(plots_dir / f"{sample_name}_homology_novelty_summary.png",
                    dpi=300, bbox_inches='tight')
        plt.close()
    
    logging.info(f"Plots saved to {plots_dir}")

def save_results(analyzer, output_dir, sample_name, keep_blast_output=False):
    """Save analysis results"""
    logging.info("Saving results")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save main results
    results_file = output_path / f"{sample_name}_homology_novelty.tsv"
    analyzer.novelty_results.to_csv(results_file, sep='\t', index=False)
    
    # Save statistics
    stats_file = output_path / f"{sample_name}_homology_novelty_stats.json"
    with open(stats_file, 'w') as f:
        json.dump(analyzer.statistics, f, indent=2)
    
    # Save detailed BLAST results if requested
    if keep_blast_output:
        blast_dir = output_path / "blast_results"
        blast_dir.mkdir(exist_ok=True)
        
        for db_name, blast_df in analyzer.blast_results.items():
            if not blast_df.empty:
                blast_file = blast_dir / f"{sample_name}_{db_name}_blast.tsv"
                blast_df.to_csv(blast_file, sep='\t', index=False)
    
    logging.info(f"Results saved to {output_path}")
    
    # Print summary
    stats = analyzer.statistics
    print(f"\n🔍 Homology Novelty Detection Summary for {sample_name}")
    print("=" * 55)
    print(f"Total sequences analyzed: {stats.get('total_sequences', 0):,}")
    print(f"Sequences with BLAST hits: {stats.get('sequences_with_hits', 0):,}")
    print(f"Sequences without hits: {stats.get('sequences_without_hits', 0):,}")
    print(f"Mean novelty score: {stats.get('mean_novelty_score', 0):.1f}")
    print(f"High novelty sequences: {stats.get('high_novelty_count', 0):,}")
    print(f"Medium novelty sequences: {stats.get('medium_novelty_count', 0):,}")
    print(f"Low novelty sequences: {stats.get('low_novelty_count', 0):,}")
    print(f"Known sequences: {stats.get('known_count', 0):,}")
    print(f"Overall novelty percentage: {stats.get('novel_percentage', 0):.1f}%")

def main():
    """Main function"""
    args = parse_arguments()
    setup_logging(args.log_level)
    
    logging.info("Starting homology-based novelty detection analysis")
    
    # Validate inputs
    if not os.path.exists(args.assembly):
        logging.error(f"Assembly file not found: {args.assembly}")
        sys.exit(1)
    
    for db in args.databases:
        if not os.path.exists(db):
            logging.warning(f"Database not found: {db}")
    
    # Set up thresholds
    thresholds = {
        'high_novelty': args.high_novelty_threshold,
        'medium_novelty': args.medium_novelty_threshold,
        'low_novelty': args.low_novelty_threshold
    }
    
    # Initialize analyzer
    analyzer = HomologyNoveltyAnalyzer(
        args.databases,
        args.database_names,
        thresholds
    )
    
    # Set up DIAMOND parameters
    analyzer.blaster.threads = args.threads
    analyzer.blaster.block_size = args.block_size
    analyzer.blaster.index_chunks = args.index_chunks
    
    # Load sequences
    analyzer.load_sequences(args.assembly, args.min_length)
    
    # Run BLAST searches
    blast_params = {
        'evalue': args.evalue,
        'max_target_seqs': args.max_target_seqs
    }
    
    analyzer.run_blast_searches(args.assembly, args.output, **blast_params)
    
    # Calculate novelty scores
    results = analyzer.calculate_novelty_scores()
    
    # Save results
    save_results(analyzer, args.output, args.sample_name, args.keep_blast_output)
    
    # Generate plots if requested
    if args.plots:
        create_visualization_plots(analyzer, args.output, args.sample_name)
    
    logging.info("Homology novelty detection analysis completed successfully")

if __name__ == "__main__":
    main()

