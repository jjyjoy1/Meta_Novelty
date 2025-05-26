#!/usr/bin/env python3
"""
📊 Phase 1 Summary Report Generator
Part of the Metagenomics Pipeline - Phase 1

Generates a comprehensive HTML report summarizing all Phase 1 analysis results.
"""

import argparse
import os
import sys
import logging
import json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import base64
from io import BytesIO
import warnings

warnings.filterwarnings('ignore')

def setup_logging(log_level="INFO"):
    """Setup logging configuration"""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Generate comprehensive Phase 1 summary report"
    )
    parser.add_argument(
        "--sample-name",
        required=True,
        help="Sample name"
    )
    parser.add_argument(
        "--results-dir",
        required=True,
        help="Results directory containing all Phase 1 outputs"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for report"
    )
    parser.add_argument(
        "--template",
        help="Custom HTML template file"
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )
    
    return parser.parse_args()

class Phase1ReportGenerator:
    """Main class for generating Phase 1 reports"""
    
    def __init__(self, sample_name, results_dir, output_dir):
        self.sample_name = sample_name
        self.results_dir = Path(results_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Data containers
        self.data = {}
        self.plots = {}
        self.summary_stats = {}
        
        # Set plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
    def load_all_data(self):
        """Load all available data from Phase 1 results"""
        logging.info("Loading Phase 1 analysis results")
        
        # Assembly statistics
        self.load_assembly_data()
        
        # Binning results
        self.load_binning_data()
        
        # Taxonomy data
        self.load_taxonomy_data()
        
        # Novelty detection results
        self.load_novelty_data()
        
        # Abundance data
        self.load_abundance_data()
        
        # Quality control data
        self.load_qc_data()
    
    def load_assembly_data(self):
        """Load assembly statistics"""
        try:
            stats_file = self.results_dir / "02_assembly" / f"{self.sample_name}_assembly_stats.json"
            if stats_file.exists():
                with open(stats_file, 'r') as f:
                    self.data['assembly'] = json.load(f)
                logging.info("Loaded assembly statistics")
            else:
                logging.warning("Assembly statistics file not found")
                self.data['assembly'] = {}
                
        except Exception as e:
            logging.error(f"Error loading assembly data: {e}")
            self.data['assembly'] = {}
    
    def load_binning_data(self):
        """Load binning and bin quality data"""
        try:
            # Bin quality
            quality_file = self.results_dir / "03_binning" / f"{self.sample_name}_bin_quality.tsv"
            if quality_file.exists():
                self.data['bin_quality'] = pd.read_csv(quality_file, sep='\t')
                logging.info(f"Loaded quality data for {len(self.data['bin_quality'])} bins")
            else:
                logging.warning("Bin quality file not found")
                self.data['bin_quality'] = pd.DataFrame()
                
        except Exception as e:
            logging.error(f"Error loading binning data: {e}")
            self.data['bin_quality'] = pd.DataFrame()
    
    def load_taxonomy_data(self):
        """Load taxonomy classification data"""
        try:
            # GTDB-Tk results
            gtdb_file = self.results_dir / "04_taxonomy" / f"{self.sample_name}_gtdbtk_summary.tsv"
            if gtdb_file.exists():
                self.data['taxonomy'] = pd.read_csv(gtdb_file, sep='\t')
                logging.info(f"Loaded taxonomy for {len(self.data['taxonomy'])} MAGs")
            else:
                logging.warning("GTDB-Tk taxonomy file not found")
                self.data['taxonomy'] = pd.DataFrame()
                
        except Exception as e:
            logging.error(f"Error loading taxonomy data: {e}")
            self.data['taxonomy'] = pd.DataFrame()
    
    def load_novelty_data(self):
        """Load novelty detection results"""
        try:
            # Combined novelty results
            novelty_file = self.results_dir / "05_novelty" / f"{self.sample_name}_combined_novelty.tsv"
            if novelty_file.exists():
                self.data['novelty'] = pd.read_csv(novelty_file, sep='\t')
                logging.info(f"Loaded novelty data for {len(self.data['novelty'])} contigs")
            else:
                logging.warning("Combined novelty file not found")
                self.data['novelty'] = pd.DataFrame()
            
            # Novelty statistics
            novelty_stats_file = self.results_dir / "05_novelty" / f"{self.sample_name}_combined_novelty_stats.json"
            if novelty_stats_file.exists():
                with open(novelty_stats_file, 'r') as f:
                    self.data['novelty_stats'] = json.load(f)
            else:
                self.data['novelty_stats'] = {}
                
        except Exception as e:
            logging.error(f"Error loading novelty data: {e}")
            self.data['novelty'] = pd.DataFrame()
            self.data['novelty_stats'] = {}
    
    def load_abundance_data(self):
        """Load abundance calculation results"""
        try:
            # Contig abundance
            contig_abundance_file = self.results_dir / "06_abundance" / f"{self.sample_name}_contig_abundance.tsv"
            if contig_abundance_file.exists():
                self.data['contig_abundance'] = pd.read_csv(contig_abundance_file, sep='\t')
                logging.info(f"Loaded contig abundance for {len(self.data['contig_abundance'])} contigs")
            else:
                logging.warning("Contig abundance file not found")
                self.data['contig_abundance'] = pd.DataFrame()
            
            # MAG abundance
            mag_abundance_file = self.results_dir / "06_abundance" / f"{self.sample_name}_mag_abundance.tsv"
            if mag_abundance_file.exists():
                self.data['mag_abundance'] = pd.read_csv(mag_abundance_file, sep='\t')
                logging.info(f"Loaded MAG abundance for {len(self.data['mag_abundance'])} MAGs")
            else:
                logging.warning("MAG abundance file not found")
                self.data['mag_abundance'] = pd.DataFrame()
                
        except Exception as e:
            logging.error(f"Error loading abundance data: {e}")
            self.data['contig_abundance'] = pd.DataFrame()
            self.data['mag_abundance'] = pd.DataFrame()
    
    def load_qc_data(self):
        """Load quality control data"""
        try:
            # FastP statistics (if available)
            fastp_file = self.results_dir / "01_quality_control" / f"{self.sample_name}_fastp_report.json"
            if fastp_file.exists():
                with open(fastp_file, 'r') as f:
                    self.data['qc'] = json.load(f)
                logging.info("Loaded QC statistics")
            else:
                logging.warning("QC statistics file not found")
                self.data['qc'] = {}
                
        except Exception as e:
            logging.error(f"Error loading QC data: {e}")
            self.data['qc'] = {}
    
    def calculate_summary_statistics(self):
        """Calculate overall summary statistics"""
        logging.info("Calculating summary statistics")
        
        self.summary_stats = {
            'sample_name': self.sample_name,
            'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'pipeline_version': 'Phase 1 v1.0'
        }
        
        # Assembly stats
        assembly_data = self.data.get('assembly', {})
        self.summary_stats.update({
            'total_contigs': assembly_data.get('total_contigs', 0),
            'total_assembly_length': assembly_data.get('total_length', 0),
            'n50': assembly_data.get('n50', 0),
            'max_contig_length': assembly_data.get('max_contig_length', 0),
            'mean_gc_content': assembly_data.get('mean_gc_content', 0)
        })
        
        # Binning stats
        bin_quality = self.data.get('bin_quality', pd.DataFrame())
        if not bin_quality.empty:
            high_quality_bins = len(bin_quality[
                (bin_quality.get('Completeness', 0) >= 90) & 
                (bin_quality.get('Contamination', 100) < 5)
            ])
            medium_quality_bins = len(bin_quality[
                (bin_quality.get('Completeness', 0) >= 50) & 
                (bin_quality.get('Contamination', 100) <= 10)
            ]) - high_quality_bins
        else:
            high_quality_bins = medium_quality_bins = 0
        
        self.summary_stats.update({
            'total_bins': len(bin_quality),
            'high_quality_bins': high_quality_bins,
            'medium_quality_bins': medium_quality_bins
        })
        
        # Taxonomy stats
        taxonomy_data = self.data.get('taxonomy', pd.DataFrame())
        self.summary_stats['classified_mags'] = len(taxonomy_data)
        
        # Novelty stats
        novelty_data = self.data.get('novelty', pd.DataFrame())
        if not novelty_data.empty and 'combined_novelty_level' in novelty_data.columns:
            high_novelty = len(novelty_data[novelty_data['combined_novelty_level'] == 'High'])
            self.summary_stats['high_novelty_contigs'] = high_novelty
        else:
            self.summary_stats['high_novelty_contigs'] = 0
        
        # Abundance stats
        mag_abundance = self.data.get('mag_abundance', pd.DataFrame())
        if not mag_abundance.empty and 'relative_abundance' in mag_abundance.columns:
            self.summary_stats['most_abundant_mag'] = mag_abundance.iloc[0].get('mag_id', 'Unknown')
            self.summary_stats['max_mag_abundance'] = mag_abundance['relative_abundance'].max()
        else:
            self.summary_stats['most_abundant_mag'] = 'Unknown'
            self.summary_stats['max_mag_abundance'] = 0
    
    def generate_plots(self):
        """Generate all summary plots"""
        logging.info("Generating summary plots")
        
        # Assembly overview plot
        self.create_assembly_plot()
        
        # Bin quality plot
        self.create_bin_quality_plot()
        
        # Taxonomy composition plot
        self.create_taxonomy_plot()
        
        # Novelty distribution plot
        self.create_novelty_plot()
        
        # Abundance overview plot
        self.create_abundance_plot()
    
    def create_assembly_plot(self):
        """Create assembly statistics plot"""
        assembly_data = self.data.get('assembly', {})
        
        if not assembly_data:
            self.plots['assembly'] = None
            return
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # Contig size categories
        categories = ['≥1kb', '≥10kb', '≥100kb', '≥1Mb']
        counts = [
            assembly_data.get('contigs_over_1kb', 0),
            assembly_data.get('contigs_over_10kb', 0),
            assembly_data.get('contigs_over_100kb', 0),
            assembly_data.get('contigs_over_1mb', 0)
        ]
        
        ax1.bar(categories, counts, color='skyblue', edgecolor='black')
        ax1.set_ylabel('Number of Contigs')
        ax1.set_title('Assembly Contiguity')
        ax1.grid(True, alpha=0.3)
        
        # Add counts on bars
        for i, count in enumerate(counts):
            if count > 0:
                ax1.text(i, count + max(counts) * 0.01, str(count), 
                        ha='center', va='bottom')
        
        # N-statistics
        n_stats = ['N50', 'N75', 'N90']
        n_values = [
            assembly_data.get('n50', 0),
            assembly_data.get('n75', 0),
            assembly_data.get('n90', 0)
        ]
        
        ax2.bar(n_stats, n_values, color='lightcoral', edgecolor='black')
        ax2.set_ylabel('Contig Length (bp)')
        ax2.set_title('N-Statistics')
        ax2.grid(True, alpha=0.3)
        
        # Assembly metrics
        metrics = ['Total Length\n(Mb)', 'Total Contigs', 'Max Contig\n(kb)', 'Mean GC\n(%)']
        values = [
            assembly_data.get('total_length', 0) / 1e6,
            assembly_data.get('total_contigs', 0),
            assembly_data.get('max_contig_length', 0) / 1000,
            assembly_data.get('mean_gc_content', 0)
        ]
        
        colors = ['green', 'blue', 'orange', 'purple']
        bars = ax3.bar(metrics, values, color=colors, alpha=0.7, edgecolor='black')
        ax3.set_title('Assembly Overview')
        ax3.grid(True, alpha=0.3)
        
        # Add values on bars
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height + max(values) * 0.01,
                    f'{value:.1f}', ha='center', va='bottom')
        
        # Assembly quality assessment
        quality_score = self.calculate_assembly_quality_score(assembly_data)
        colors = ['red' if quality_score < 50 else 'orange' if quality_score < 75 else 'green']
        
        ax4.pie([quality_score, 100-quality_score], labels=['Quality Score', ''], 
               colors=colors + ['lightgray'], startangle=90,
               wedgeprops=dict(width=0.3))
        ax4.set_title(f'Assembly Quality\n{quality_score:.1f}/100')
        
        plt.tight_layout()
        self.plots['assembly'] = self.fig_to_base64(fig)
        plt.close(fig)
    
    def calculate_assembly_quality_score(self, assembly_data):
        """Calculate assembly quality score (0-100)"""
        score = 0
        
        # N50 contribution (40 points max)
        n50 = assembly_data.get('n50', 0)
        if n50 > 50000:
            score += 40
        elif n50 > 10000:
            score += 30
        elif n50 > 1000:
            score += 20
        else:
            score += 10
        
        # Total length contribution (20 points max)
        total_length = assembly_data.get('total_length', 0)
        if total_length > 50e6:  # 50 Mb
            score += 20
        elif total_length > 10e6:  # 10 Mb
            score += 15
        elif total_length > 1e6:   # 1 Mb
            score += 10
        else:
            score += 5
        
        # Contiguity contribution (20 points max)
        long_contigs = assembly_data.get('contigs_over_10kb', 0)
        total_contigs = assembly_data.get('total_contigs', 1)
        if long_contigs / total_contigs > 0.1:
            score += 20
        elif long_contigs / total_contigs > 0.05:
            score += 15
        else:
            score += 10
        
        # GC content contribution (20 points max)
        gc_content = assembly_data.get('mean_gc_content', 50)
        if 30 <= gc_content <= 70:  # Reasonable range
            score += 20
        elif 20 <= gc_content <= 80:
            score += 15
        else:
            score += 10
        
        return min(score, 100)
    
    def create_bin_quality_plot(self):
        """Create bin quality assessment plot"""
        bin_quality = self.data.get('bin_quality', pd.DataFrame())
        
        if bin_quality.empty:
            self.plots['bin_quality'] = None
            return
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # Quality scatter plot
        if 'Completeness' in bin_quality.columns and 'Contamination' in bin_quality.columns:
            completeness = bin_quality['Completeness']
            contamination = bin_quality['Contamination']
            
            # Color by quality
            colors = []
            for c, cont in zip(completeness, contamination):
                if c >= 90 and cont < 5:
                    colors.append('green')  # High quality
                elif c >= 50 and cont <= 10:
                    colors.append('orange')  # Medium quality
                else:
                    colors.append('red')    # Low quality
            
            ax1.scatter(completeness, contamination, c=colors, alpha=0.7, s=50)
            ax1.set_xlabel('Completeness (%)')
            ax1.set_ylabel('Contamination (%)')
            ax1.set_title('Bin Quality Distribution')
            ax1.grid(True, alpha=0.3)
            
            # Add quality thresholds
            ax1.axhline(y=10, color='red', linestyle='--', alpha=0.5, label='10% contamination')
            ax1.axhline(y=5, color='orange', linestyle='--', alpha=0.5, label='5% contamination')
            ax1.axvline(x=50, color='red', linestyle='--', alpha=0.5, label='50% completeness')
            ax1.axvline(x=90, color='orange', linestyle='--', alpha=0.5, label='90% completeness')
            ax1.legend()
        
        # Quality categories
        high_quality = len(bin_quality[
            (bin_quality.get('Completeness', 0) >= 90) & 
            (bin_quality.get('Contamination', 100) < 5)
        ])
        medium_quality = len(bin_quality[
            (bin_quality.get('Completeness', 0) >= 50) & 
            (bin_quality.get('Contamination', 100) <= 10)
        ]) - high_quality
        low_quality = len(bin_quality) - high_quality - medium_quality
        
        categories = ['High\n(>90% complete,\n<5% contamination)', 
                     'Medium\n(>50% complete,\n<10% contamination)', 
                     'Low\n(other)']
        counts = [high_quality, medium_quality, low_quality]
        colors = ['green', 'orange', 'red']
        
        ax2.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
        ax2.set_ylabel('Number of Bins')
        ax2.set_title('Quality Categories')
        ax2.grid(True, alpha=0.3)
        
        # Add counts on bars
        for i, count in enumerate(counts):
            if count > 0:
                ax2.text(i, count + max(counts) * 0.01, str(count), 
                        ha='center', va='bottom')
        
        # Completeness distribution
        if 'Completeness' in bin_quality.columns:
            ax3.hist(bin_quality['Completeness'], bins=20, alpha=0.7, 
                    color='skyblue', edgecolor='black')
            ax3.axvline(bin_quality['Completeness'].mean(), color='red', 
                       linestyle='--', label=f"Mean: {bin_quality['Completeness'].mean():.1f}%")
            ax3.set_xlabel('Completeness (%)')
            ax3.set_ylabel('Frequency')
            ax3.set_title('Completeness Distribution')
            ax3.legend()
            ax3.grid(True, alpha=0.3)
        
        # Contamination distribution
        if 'Contamination' in bin_quality.columns:
            ax4.hist(bin_quality['Contamination'], bins=20, alpha=0.7, 
                    color='lightcoral', edgecolor='black')
            ax4.axvline(bin_quality['Contamination'].mean(), color='red', 
                       linestyle='--', label=f"Mean: {bin_quality['Contamination'].mean():.1f}%")
            ax4.set_xlabel('Contamination (%)')
            ax4.set_ylabel('Frequency')
            ax4.set_title('Contamination Distribution')
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        self.plots['bin_quality'] = self.fig_to_base64(fig)
        plt.close(fig)
    
    def create_taxonomy_plot(self):
        """Create taxonomy composition plot"""
        taxonomy_data = self.data.get('taxonomy', pd.DataFrame())
        
        if taxonomy_data.empty:
            self.plots['taxonomy'] = None
            return
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        
        # Extract phylum information (assuming GTDB format)
        if 'classification' in taxonomy_data.columns:
            phyla = []
            for classification in taxonomy_data['classification']:
                if pd.isna(classification):
                    phyla.append('Unclassified')
                    continue
                
                # Parse GTDB classification
                parts = str(classification).split(';')
                phylum = 'Unclassified'
                for part in parts:
                    if part.startswith('p__'):
                        phylum = part[3:].strip()
                        if phylum == '':
                            phylum = 'Unclassified'
                        break
                phyla.append(phylum)
            
            # Count phyla
            phylum_counts = pd.Series(phyla).value_counts()
            
            # Plot top 10 phyla
            top_phyla = phylum_counts.head(10)
            if len(top_phyla) > 0:
                wedges, texts, autotexts = ax1.pie(top_phyla.values, labels=top_phyla.index, 
                                                  autopct='%1.1f%%', startangle=90)
                ax1.set_title('Taxonomic Composition (Phylum)')
        
        # Classification success rate
        if 'classification' in taxonomy_data.columns:
            classified = len(taxonomy_data.dropna(subset=['classification']))
            unclassified = len(taxonomy_data) - classified
            
            ax2.pie([classified, unclassified], labels=['Classified', 'Unclassified'],
                   colors=['green', 'red'], autopct='%1.1f%%', startangle=90)
            ax2.set_title('Classification Success')
        
        plt.tight_layout()
        self.plots['taxonomy'] = self.fig_to_base64(fig)
        plt.close(fig)
    
    def create_novelty_plot(self):
        """Create novelty detection summary plot"""
        novelty_data = self.data.get('novelty', pd.DataFrame())
        novelty_stats = self.data.get('novelty_stats', {})
        
        if novelty_data.empty:
            self.plots['novelty'] = None
            return
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # Combined novelty score distribution
        if 'combined_novelty_score' in novelty_data.columns:
            scores = novelty_data['combined_novelty_score']
            ax1.hist(scores, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
            ax1.axvline(scores.mean(), color='red', linestyle='--',
                       label=f'Mean: {scores.mean():.1f}')
            ax1.set_xlabel('Combined Novelty Score')
            ax1.set_ylabel('Frequency')
            ax1.set_title('Novelty Score Distribution')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
        
        # Novelty level categories
        if 'combined_novelty_level' in novelty_data.columns:
            level_counts = novelty_data['combined_novelty_level'].value_counts()
            colors = {'High': 'red', 'Medium': 'orange', 'Low': 'green'}
            plot_colors = [colors.get(level, 'gray') for level in level_counts.index]
            
            wedges, texts, autotexts = ax2.pie(level_counts.values, labels=level_counts.index,
                                              colors=plot_colors, autopct='%1.1f%%', startangle=90)
            ax2.set_title('Novelty Level Distribution')
        
        # Method agreement
        if 'method_agreement' in novelty_data.columns:
            agreement_counts = novelty_data['method_agreement'].value_counts()
            colors = ['lightcoral', 'lightgreen']
            labels = ['Disagreement', 'Agreement']
            
            ax3.pie(agreement_counts.values, labels=labels, colors=colors, 
                   autopct='%1.1f%%', startangle=90)
            ax3.set_title('Method Agreement')
        
        # High-confidence novel sequences
        if 'high_confidence_novel' in novelty_data.columns:
            high_conf_novel = novelty_data['high_confidence_novel'].sum()
            total_sequences = len(novelty_data)
            other_sequences = total_sequences - high_conf_novel
            
            ax4.pie([high_conf_novel, other_sequences], 
                   labels=['High-Confidence Novel', 'Other'],
                   colors=['red', 'lightblue'], autopct='%1.1f%%', startangle=90)
            ax4.set_title('High-Confidence Novel Sequences')
        
        plt.tight_layout()
        self.plots['novelty'] = self.fig_to_base64(fig)
        plt.close(fig)
    
    def create_abundance_plot(self):
        """Create abundance overview plot"""
        mag_abundance = self.data.get('mag_abundance', pd.DataFrame())
        
        if mag_abundance.empty:
            self.plots['abundance'] = None
            return
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # Relative abundance distribution
        if 'relative_abundance' in mag_abundance.columns:
            abundances = mag_abundance['relative_abundance']
            abundances = abundances[abundances > 0]  # Remove zeros for log scale
            
            if len(abundances) > 0:
                ax1.hist(abundances, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
                ax1.set_xlabel('Relative Abundance (%)')
                ax1.set_ylabel('Frequency')
                ax1.set_title('MAG Relative Abundance Distribution')
                ax1.set_yscale('log')
                ax1.grid(True, alpha=0.3)
        
        # Top 10 most abundant MAGs
        if 'relative_abundance' in mag_abundance.columns:
            top_mags = mag_abundance.nlargest(10, 'relative_abundance')
            
            # Create labels
            labels = []
            for _, row in top_mags.iterrows():
                mag_id = str(row['mag_id'])[:10] + "..." if len(str(row['mag_id'])) > 10 else str(row['mag_id'])
                labels.append(mag_id)
            
            y_pos = np.arange(len(labels))
            ax2.barh(y_pos, top_mags['relative_abundance'], alpha=0.7, color='lightcoral')
            ax2.set_yticks(y_pos)
            ax2.set_yticklabels(labels)
            ax2.set_xlabel('Relative Abundance (%)')
            ax2.set_title('Top 10 Most Abundant MAGs')
            ax2.grid(True, alpha=0.3)
        
        # Quality vs abundance scatter
        if 'relative_abundance' in mag_abundance.columns and 'quality_score' in mag_abundance.columns:
            ax3.scatter(mag_abundance['quality_score'], mag_abundance['relative_abundance'],
                       alpha=0.6, s=30)
            ax3.set_xlabel('Quality Score')
            ax3.set_ylabel('Relative Abundance (%)')
            ax3.set_title('Quality vs Abundance')
            ax3.set_yscale('log')
            ax3.grid(True, alpha=0.3)
        
        # Genome size distribution
        if 'total_length' in mag_abundance.columns:
            genome_sizes = mag_abundance['total_length'] / 1e6  # Convert to Mb
            
            ax4.hist(genome_sizes, bins=20, alpha=0.7, color='lightgreen', edgecolor='black')
            ax4.axvline(genome_sizes.mean(), color='red', linestyle='--',
                       label=f'Mean: {genome_sizes.mean():.1f} Mb')
            ax4.set_xlabel('Genome Size (Mb)')
            ax4.set_ylabel('Frequency')
            ax4.set_title('MAG Genome Size Distribution')
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        self.plots['abundance'] = self.fig_to_base64(fig)
        plt.close(fig)
    
    def fig_to_base64(self, fig):
        """Convert matplotlib figure to base64 string"""
        buffer = BytesIO()
        fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
        buffer.seek(0)
        image_png = buffer.getvalue()
        buffer.close()
        graphic = base64.b64encode(image_png).decode()
        return graphic
    
    def generate_html_report(self, template_file=None):
        """Generate HTML report"""
        logging.info("Generating HTML report")
        
        if template_file and os.path.exists(template_file):
            with open(template_file, 'r') as f:
                template = f.read()
        else:
            template = self.get_default_template()
        
        # Format template with data
        html_content = template.format(
            sample_name=self.sample_name,
            analysis_date=self.summary_stats.get('analysis_date', 'Unknown'),
            pipeline_version=self.summary_stats.get('pipeline_version', 'Unknown'),
            
            # Summary stats
            total_contigs=f"{self.summary_stats.get('total_contigs', 0):,}",
            total_assembly_length=f"{self.summary_stats.get('total_assembly_length', 0):,}",
            n50=f"{self.summary_stats.get('n50', 0):,}",
            max_contig_length=f"{self.summary_stats.get('max_contig_length', 0):,}",
            mean_gc_content=f"{self.summary_stats.get('mean_gc_content', 0):.1f}",
            
            total_bins=self.summary_stats.get('total_bins', 0),
            high_quality_bins=self.summary_stats.get('high_quality_bins', 0),
            medium_quality_bins=self.summary_stats.get('medium_quality_bins', 0),
            
            classified_mags=self.summary_stats.get('classified_mags', 0),
            high_novelty_contigs=self.summary_stats.get('high_novelty_contigs', 0),
            
            most_abundant_mag=self.summary_stats.get('most_abundant_mag', 'Unknown'),
            max_mag_abundance=f"{self.summary_stats.get('max_mag_abundance', 0):.2f}",
            
            # Plots
            assembly_plot=self.plots.get('assembly', ''),
            bin_quality_plot=self.plots.get('bin_quality', ''),
            taxonomy_plot=self.plots.get('taxonomy', ''),
            novelty_plot=self.plots.get('novelty', ''),
            abundance_plot=self.plots.get('abundance', '')
        )
        
        # Save HTML report
        html_file = self.output_dir / f"{self.sample_name}_phase1_summary.html"
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        # Save JSON summary
        json_file = self.output_dir / f"{self.sample_name}_phase1_summary.json"
        with open(json_file, 'w') as f:
            json.dump(self.summary_stats, f, indent=2)
        
        logging.info(f"Report saved to {html_file}")
        return str(html_file)
    
    def get_default_template(self):
        """Get default HTML template"""
        return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Phase 1 Analysis Report - {sample_name}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f4f4f4;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
        }}
        .header {{
            text-align: center;
            border-bottom: 3px solid #333;
            padding-bottom: 20px;
            margin-bottom: 30px;
        }}
        .header h1 {{
            color: #333;
            margin-bottom: 10px;
        }}
        .section {{
            margin-bottom: 40px;
        }}
        .section h2 {{
            color: #333;
            border-bottom: 2px solid #ddd;
            padding-bottom: 10px;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .stat-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #007bff;
        }}
        .stat-value {{
            font-size: 24px;
            font-weight: bold;
            color: #007bff;
        }}
        .stat-label {{
            color: #666;
            font-size: 14px;
        }}
        .plot-container {{
            text-align: center;
            margin: 30px 0;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 8px;
        }}
        .plot-container h3 {{
            margin-top: 0;
            color: #333;
        }}
        .plot-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 5px;
        }}
        .footer {{
            text-align: center;
            margin-top: 50px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #666;
        }}
        .quality-indicator {{
            display: inline-block;
            padding: 5px 10px;
            border-radius: 15px;
            color: white;
            font-weight: bold;
        }}
        .high-quality {{ background-color: #28a745; }}
        .medium-quality {{ background-color: #ffc107; color: #333; }}
        .low-quality {{ background-color: #dc3545; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Metagenomics Analysis Report</h1>
            <h2>Phase 1: Single-Sample Processing</h2>
            <p><strong>Sample:</strong> {sample_name}</p>
            <p><strong>Analysis Date:</strong> {analysis_date}</p>
            <p><strong>Pipeline Version:</strong> {pipeline_version}</p>
        </div>

        <div class="section">
            <h2>📊 Summary Statistics</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{total_contigs}</div>
                    <div class="stat-label">Total Contigs</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{total_assembly_length} bp</div>
                    <div class="stat-label">Total Assembly Length</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{n50} bp</div>
                    <div class="stat-label">N50</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{mean_gc_content}%</div>
                    <div class="stat-label">Mean GC Content</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{total_bins}</div>
                    <div class="stat-label">Total MAGs</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{high_quality_bins}</div>
                    <div class="stat-label">High-Quality MAGs</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{classified_mags}</div>
                    <div class="stat-label">Classified MAGs</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{high_novelty_contigs}</div>
                    <div class="stat-label">High-Novelty Contigs</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>🔬 Assembly Analysis</h2>
            <div class="plot-container">
                <h3>Assembly Statistics Overview</h3>
                {assembly_plot_img}
            </div>
        </div>

        <div class="section">
            <h2>🧪 MAG Quality Assessment</h2>
            <div class="plot-container">
                <h3>Bin Quality Distribution</h3>
                {bin_quality_plot_img}
            </div>
        </div>

        <div class="section">
            <h2>🌳 Taxonomic Classification</h2>
            <div class="plot-container">
                <h3>Taxonomic Composition</h3>
                {taxonomy_plot_img}
            </div>
        </div>

        <div class="section">
            <h2>🔍 Novelty Detection</h2>
            <div class="plot-container">
                <h3>Novelty Analysis Results</h3>
                {novelty_plot_img}
            </div>
        </div>

        <div class="section">
            <h2>📈 Abundance Analysis</h2>
            <div class="plot-container">
                <h3>MAG Abundance Overview</h3>
                {abundance_plot_img}
            </div>
            <p><strong>Most Abundant MAG:</strong> {most_abundant_mag} ({max_mag_abundance}%)</p>
        </div>

        <div class="footer">
            <p>Generated by the Metagenomics Pipeline Phase 1</p>
            <p>For questions or support, please refer to the pipeline documentation</p>
        </div>
    </div>

    <script>
        // Add plot images
        document.addEventListener('DOMContentLoaded', function() {{
            const plots = {{
                'assembly_plot_img': '{assembly_plot}',
                'bin_quality_plot_img': '{bin_quality_plot}',
                'taxonomy_plot_img': '{taxonomy_plot}',
                'novelty_plot_img': '{novelty_plot}',
                'abundance_plot_img': '{abundance_plot}'
            }};
            
            for (const [placeholder, base64_data] of Object.entries(plots)) {{
                const element = document.querySelector('.' + placeholder.replace('_img', '') + ' h3');
                if (element && base64_data) {{
                    const img = document.createElement('img');
                    img.src = 'data:image/png;base64,' + base64_data;
                    img.style.maxWidth = '100%';
                    element.parentNode.appendChild(img);
                }}
            }}
        }});
    </script>
</body>
</html>
        """.replace('{assembly_plot_img}', '').replace('{bin_quality_plot_img}', '').replace('{taxonomy_plot_img}', '').replace('{novelty_plot_img}', '').replace('{abundance_plot_img}', '')

def main():
    """Main function"""
    args = parse_arguments()
    setup_logging(args.log_level)
    
    logging.info("Starting Phase 1 report generation")
    
    # Initialize report generator
    generator = Phase1ReportGenerator(
        args.sample_name,
        args.results_dir,
        args.output
    )
    
    # Load all data
    generator.load_all_data()
    
    # Calculate summary statistics
    generator.calculate_summary_statistics()
    
    # Generate plots
    generator.generate_plots()
    
    # Generate HTML report
    report_file = generator.generate_html_report(args.template)
    
    print(f"\n📊 Phase 1 Report Generated Successfully!")
    print("=" * 50)
    print(f"Sample: {args.sample_name}")
    print(f"Report: {report_file}")
    print(f"Total contigs: {generator.summary_stats.get('total_contigs', 0):,}")
    print(f"Total MAGs: {generator.summary_stats.get('total_bins', 0)}")
    print(f"High-quality MAGs: {generator.summary_stats.get('high_quality_bins', 0)}")
    
    logging.info("Phase 1 report generation completed successfully")

if __name__ == "__main__":
    main()


