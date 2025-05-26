#!/usr/bin/env python3
"""
🧪 MAG (Metagenome-Assembled Genome) Abundance Calculator
Part of the Metagenomics Pipeline - Phase 1

Calculates abundance for MAGs (bins) by aggregating contig-level abundance
data and applying quality-based weighting and normalization.
"""

import argparse
import os
import sys
import logging
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import json
from collections import defaultdict
import glob

def setup_logging(log_level="INFO"):
    """Setup logging configuration"""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Calculate MAG abundance from contig abundance and binning results"
    )
    parser.add_argument(
        "--contig-abundance", "-c",
        required=True,
        help="Path to contig abundance TSV file"
    )
    parser.add_argument(
        "--bins-dir", "-b",
        required=True,
        help="Directory containing MAG/bin FASTA files"
    )
    parser.add_argument(
        "--bin-quality",
        help="Path to bin quality assessment file (CheckM or similar)"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "--sample-name",
        default="sample",
        help="Sample name for output files"
    )
    parser.add_argument(
        "--min-completeness",
        type=float,
        default=50.0,
        help="Minimum completeness for MAG inclusion (default: 50%)"
    )
    parser.add_argument(
        "--max-contamination",
        type=float,
        default=10.0,
        help="Maximum contamination for MAG inclusion (default: 10%)"
    )
    parser.add_argument(
        "--quality-weighting",
        action="store_true",
        help="Apply quality-based weighting to abundance calculations"
    )
    parser.add_argument(
        "--normalization",
        choices=["relative", "total_coverage", "genome_size", "all"],
        default="relative",
        help="Normalization method for MAG abundance"
    )
    parser.add_argument(
        "--min-contig-length",
        type=int,
        default=1000,
        help="Minimum contig length to include in MAG abundance"
    )
    parser.add_argument(
        "--taxonomy-file",
        help="Path to taxonomy classification file for MAGs"
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

class MAGAbundanceCalculator:
    """Main class for calculating MAG abundance"""
    
    def __init__(self, min_completeness=50.0, max_contamination=10.0,
                 quality_weighting=False, min_contig_length=1000):
        self.min_completeness = min_completeness
        self.max_contamination = max_contamination
        self.quality_weighting = quality_weighting
        self.min_contig_length = min_contig_length
        
        self.contig_abundance = None
        self.mag_composition = {}
        self.mag_quality = {}
        self.mag_taxonomy = {}
        self.mag_abundance = {}
        self.statistics = {}
    
    def load_contig_abundance(self, abundance_file):
        """Load contig abundance data"""
        logging.info(f"Loading contig abundance from {abundance_file}")
        
        try:
            self.contig_abundance = pd.read_csv(abundance_file, sep='\t')
            logging.info(f"Loaded abundance data for {len(self.contig_abundance)} contigs")
            
            # Ensure required columns are present
            required_cols = ['contig_id']
            missing_cols = [col for col in required_cols if col not in self.contig_abundance.columns]
            
            if missing_cols:
                logging.error(f"Missing required columns: {missing_cols}")
                sys.exit(1)
            
            # Identify abundance columns
            self.abundance_columns = [col for col in self.contig_abundance.columns 
                                    if any(metric in col.lower() for metric in 
                                          ['coverage', 'abundance', 'tpm', 'rpkm'])]
            
            logging.info(f"Found abundance columns: {self.abundance_columns}")
            
        except Exception as e:
            logging.error(f"Error loading contig abundance: {e}")
            sys.exit(1)
    
    def load_mag_composition(self, bins_dir):
        """Load MAG composition from bin FASTA files"""
        logging.info(f"Loading MAG composition from {bins_dir}")
        
        bins_path = Path(bins_dir)
        if not bins_path.exists():
            logging.error(f"Bins directory not found: {bins_dir}")
            sys.exit(1)
        
        # Find all FASTA files in bins directory
        fasta_patterns = ['*.fasta', '*.fa', '*.fas', '*.fna']
        fasta_files = []
        
        for pattern in fasta_patterns:
            fasta_files.extend(glob.glob(str(bins_path / pattern)))
        
        if not fasta_files:
            logging.error(f"No FASTA files found in {bins_dir}")
            sys.exit(1)
        
        logging.info(f"Found {len(fasta_files)} bin files")
        
        # Load contig composition for each MAG
        for fasta_file in fasta_files:
            mag_id = Path(fasta_file).stem
            contigs = []
            total_length = 0
            
            try:
                with open(fasta_file, 'r') as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        if len(record.seq) >= self.min_contig_length:
                            contigs.append(record.id)
                            total_length += len(record.seq)
                
                self.mag_composition[mag_id] = {
                    'contigs': contigs,
                    'num_contigs': len(contigs),
                    'total_length': total_length,
                    'fasta_file': fasta_file
                }
                
            except Exception as e:
                logging.warning(f"Error reading {fasta_file}: {e}")
                continue
        
        logging.info(f"Loaded composition for {len(self.mag_composition)} MAGs")
    
    def load_mag_quality(self, quality_file):
        """Load MAG quality assessment data"""
        if not quality_file or not os.path.exists(quality_file):
            logging.warning("MAG quality file not provided or not found")
            # Create default quality data
            for mag_id in self.mag_composition:
                self.mag_quality[mag_id] = {
                    'completeness': 70.0,  # Default values
                    'contamination': 5.0,
                    'strain_heterogeneity': 0.0,
                    'quality_score': 65.0
                }
            return
        
        logging.info(f"Loading MAG quality from {quality_file}")
        
        try:
            # Try to detect file format and load accordingly
            if quality_file.endswith('.tsv'):
                quality_df = pd.read_csv(quality_file, sep='\t')
            else:
                quality_df = pd.read_csv(quality_file)
            
            # Standardize column names (handle different CheckM output formats)
            column_mapping = {
                'Bin Id': 'mag_id',
                'bin_id': 'mag_id',
                'MAG': 'mag_id',
                'Completeness': 'completeness',
                'completeness': 'completeness',
                'Contamination': 'contamination',
                'contamination': 'contamination',
                'Strain heterogeneity': 'strain_heterogeneity',
                'strain_heterogeneity': 'strain_heterogeneity'
            }
            
            for old_name, new_name in column_mapping.items():
                if old_name in quality_df.columns:
                    quality_df = quality_df.rename(columns={old_name: new_name})
            
            # Load quality data for each MAG
            for _, row in quality_df.iterrows():
                mag_id = str(row['mag_id'])
                
                # Remove common suffixes to match bin file names
                for suffix in ['.fa', '.fasta', '.fna']:
                    if mag_id.endswith(suffix):
                        mag_id = mag_id[:-len(suffix)]
                
                completeness = float(row.get('completeness', 0))
                contamination = float(row.get('contamination', 0))
                strain_het = float(row.get('strain_heterogeneity', 0))
                
                # Calculate quality score (similar to MIMAG standards)
                quality_score = completeness - 5 * contamination
                
                self.mag_quality[mag_id] = {
                    'completeness': completeness,
                    'contamination': contamination,
                    'strain_heterogeneity': strain_het,
                    'quality_score': quality_score
                }
            
            logging.info(f"Loaded quality data for {len(self.mag_quality)} MAGs")
            
        except Exception as e:
            logging.error(f"Error loading MAG quality: {e}")
    
    def load_mag_taxonomy(self, taxonomy_file):
        """Load MAG taxonomy data"""
        if not taxonomy_file or not os.path.exists(taxonomy_file):
            logging.warning("MAG taxonomy file not provided")
            return
        
        logging.info(f"Loading MAG taxonomy from {taxonomy_file}")
        
        try:
            if taxonomy_file.endswith('.tsv'):
                tax_df = pd.read_csv(taxonomy_file, sep='\t')
            else:
                tax_df = pd.read_csv(taxonomy_file)
            
            # Standardize column names
            column_mapping = {
                'user_genome': 'mag_id',
                'mag_id': 'mag_id',
                'bin_id': 'mag_id',
                'classification': 'taxonomy',
                'taxonomy': 'taxonomy',
                'lineage': 'taxonomy'
            }
            
            for old_name, new_name in column_mapping.items():
                if old_name in tax_df.columns:
                    tax_df = tax_df.rename(columns={old_name: new_name})
            
            # Load taxonomy for each MAG
            for _, row in tax_df.iterrows():
                mag_id = str(row['mag_id'])
                
                # Remove common suffixes
                for suffix in ['.fa', '.fasta', '.fna']:
                    if mag_id.endswith(suffix):
                        mag_id = mag_id[:-len(suffix)]
                
                taxonomy = row.get('taxonomy', 'Unclassified')
                
                # Parse taxonomy string
                tax_levels = self.parse_taxonomy_string(taxonomy)
                
                self.mag_taxonomy[mag_id] = {
                    'full_taxonomy': taxonomy,
                    **tax_levels
                }
            
            logging.info(f"Loaded taxonomy for {len(self.mag_taxonomy)} MAGs")
            
        except Exception as e:
            logging.error(f"Error loading MAG taxonomy: {e}")
    
    def parse_taxonomy_string(self, taxonomy):
        """Parse taxonomy string into levels"""
        levels = {
            'domain': 'Unclassified',
            'phylum': 'Unclassified',
            'class': 'Unclassified',
            'order': 'Unclassified',
            'family': 'Unclassified',
            'genus': 'Unclassified',
            'species': 'Unclassified'
        }
        
        if pd.isna(taxonomy) or taxonomy == '':
            return levels
        
        # Handle GTDB-style taxonomy
        if ';' in taxonomy:
            parts = taxonomy.split(';')
            level_prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
            level_names = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
            
            for part in parts:
                part = part.strip()
                for prefix, level_name in zip(level_prefixes, level_names):
                    if part.startswith(prefix):
                        value = part[3:].strip()
                        if value and value != '':
                            levels[level_name] = value
                        break
        
        return levels
    
    def filter_mags_by_quality(self):
        """Filter MAGs based on quality thresholds"""
        logging.info("Filtering MAGs by quality thresholds")
        
        high_quality_mags = []
        medium_quality_mags = []
        low_quality_mags = []
        
        for mag_id in self.mag_composition:
            if mag_id not in self.mag_quality:
                logging.warning(f"No quality data for MAG {mag_id}, skipping")
                continue
            
            quality = self.mag_quality[mag_id]
            completeness = quality['completeness']
            contamination = quality['contamination']
            
            # MIMAG quality standards
            if completeness >= 90 and contamination < 5:
                high_quality_mags.append(mag_id)
            elif completeness >= self.min_completeness and contamination <= self.max_contamination:
                medium_quality_mags.append(mag_id)
            else:
                low_quality_mags.append(mag_id)
        
        logging.info(f"Quality filtering results:")
        logging.info(f"  High quality MAGs: {len(high_quality_mags)}")
        logging.info(f"  Medium quality MAGs: {len(medium_quality_mags)}")
        logging.info(f"  Low quality MAGs: {len(low_quality_mags)}")
        
        # Store quality categories
        self.quality_categories = {
            'high': high_quality_mags,
            'medium': medium_quality_mags,
            'low': low_quality_mags
        }
        
        # Return MAGs that pass quality filters
        return high_quality_mags + medium_quality_mags
    
    def calculate_mag_abundance(self, normalization_methods):
        """Calculate abundance for each MAG"""
        logging.info("Calculating MAG abundance")
        
        if isinstance(normalization_methods, str):
            if normalization_methods == "all":
                methods = ["relative", "total_coverage", "genome_size"]
            else:
                methods = [normalization_methods]
        else:
            methods = normalization_methods
        
        # Filter MAGs by quality
        valid_mags = self.filter_mags_by_quality()
        
        for mag_id in valid_mags:
            mag_data = self.mag_composition[mag_id]
            mag_contigs = mag_data['contigs']
            
            # Get abundance data for MAG contigs
            mag_contig_data = self.contig_abundance[
                self.contig_abundance['contig_id'].isin(mag_contigs)
            ]
            
            if mag_contig_data.empty:
                logging.warning(f"No abundance data found for MAG {mag_id} contigs")
                continue
            
            # Calculate abundance metrics
            abundance_metrics = {}
            
            # Raw abundance (sum of contig abundances)
            for col in self.abundance_columns:
                if col in mag_contig_data.columns:
                    abundance_metrics[f'total_{col}'] = mag_contig_data[col].sum()
                    abundance_metrics[f'mean_{col}'] = mag_contig_data[col].mean()
                    abundance_metrics[f'median_{col}'] = mag_contig_data[col].median()
            
            # Calculate coverage-weighted abundance
            if 'mean_coverage' in mag_contig_data.columns and 'length' in mag_contig_data.columns:
                # Weight by contig length
                total_length = mag_contig_data['length'].sum()
                if total_length > 0:
                    weighted_coverage = (mag_contig_data['mean_coverage'] * 
                                       mag_contig_data['length']).sum() / total_length
                    abundance_metrics['weighted_mean_coverage'] = weighted_coverage
            
            # Apply quality weighting if enabled
            if self.quality_weighting and mag_id in self.mag_quality:
                quality = self.mag_quality[mag_id]
                quality_weight = (quality['completeness'] / 100) * (1 - quality['contamination'] / 100)
                
                for metric in abundance_metrics:
                    if metric.startswith('total_') or metric.startswith('mean_'):
                        abundance_metrics[f'{metric}_quality_weighted'] = (
                            abundance_metrics[metric] * quality_weight
                        )
            
            # Store MAG-level information
            abundance_metrics.update({
                'mag_id': mag_id,
                'num_contigs': len(mag_contigs),
                'total_length': mag_data['total_length'],
                'contigs_with_abundance': len(mag_contig_data),
                'abundance_coverage': len(mag_contig_data) / len(mag_contigs) * 100
            })
            
            # Add quality information
            if mag_id in self.mag_quality:
                abundance_metrics.update(self.mag_quality[mag_id])
            
            # Add taxonomy information
            if mag_id in self.mag_taxonomy:
                abundance_metrics.update(self.mag_taxonomy[mag_id])
            
            self.mag_abundance[mag_id] = abundance_metrics
        
        # Apply normalization methods
        self.apply_normalization(methods)
        
        # Calculate summary statistics
        self.calculate_summary_statistics()
        
        logging.info(f"Calculated abundance for {len(self.mag_abundance)} MAGs")
    
    def apply_normalization(self, methods):
        """Apply normalization methods to MAG abundance"""
        logging.info(f"Applying normalization methods: {methods}")
        
        if not self.mag_abundance:
            return
        
        # Get primary abundance metric (use first available coverage column)
        primary_metric = None
        for col in self.abundance_columns:
            total_col = f'total_{col}'
            if total_col in list(self.mag_abundance.values())[0]:
                primary_metric = total_col
                break
        
        if not primary_metric:
            logging.warning("No suitable abundance metric found for normalization")
            return
        
        # Extract values for normalization
        mag_abundances = {mag_id: data[primary_metric] 
                         for mag_id, data in self.mag_abundance.items()}
        
        total_abundance = sum(mag_abundances.values())
        
        for mag_id in self.mag_abundance:
            data = self.mag_abundance[mag_id]
            abundance = mag_abundances[mag_id]
            
            if "relative" in methods:
                # Relative abundance (percentage)
                if total_abundance > 0:
                    data['relative_abundance'] = (abundance / total_abundance) * 100
                else:
                    data['relative_abundance'] = 0
            
            if "total_coverage" in methods:
                # Normalize by total coverage
                if total_abundance > 0:
                    data['normalized_by_total_coverage'] = abundance / total_abundance
                else:
                    data['normalized_by_total_coverage'] = 0
            
            if "genome_size" in methods:
                # Abundance per Mb (normalize by genome size)
                genome_size_mb = data['total_length'] / 1e6
                if genome_size_mb > 0:
                    data['abundance_per_mb'] = abundance / genome_size_mb
                else:
                    data['abundance_per_mb'] = 0
    
    def calculate_summary_statistics(self):
        """Calculate summary statistics for MAG abundance"""
        if not self.mag_abundance:
            self.statistics = {}
            return
        
        # Basic statistics
        abundance_values = []
        quality_scores = []
        genome_sizes = []
        
        for data in self.mag_abundance.values():
            if 'relative_abundance' in data:
                abundance_values.append(data['relative_abundance'])
            if 'quality_score' in data:
                quality_scores.append(data['quality_score'])
            if 'total_length' in data:
                genome_sizes.append(data['total_length'])
        
        self.statistics = {
            'total_mags': len(self.mag_abundance),
            'high_quality_mags': len(self.quality_categories.get('high', [])),
            'medium_quality_mags': len(self.quality_categories.get('medium', [])),
            'low_quality_mags': len(self.quality_categories.get('low', []))
        }
        
        if abundance_values:
            self.statistics.update({
                'mean_relative_abundance': np.mean(abundance_values),
                'median_relative_abundance': np.median(abundance_values),
                'std_relative_abundance': np.std(abundance_values),
                'max_relative_abundance': np.max(abundance_values),
                'min_relative_abundance': np.min(abundance_values)
            })
        
        if quality_scores:
            self.statistics.update({
                'mean_quality_score': np.mean(quality_scores),
                'median_quality_score': np.median(quality_scores)
            })
        
        if genome_sizes:
            self.statistics.update({
                'mean_genome_size': np.mean(genome_sizes),
                'median_genome_size': np.median(genome_sizes),
                'total_genomic_content': np.sum(genome_sizes)
            })
        
        # Taxonomic diversity
        if self.mag_taxonomy:
            taxa_levels = ['phylum', 'class', 'order', 'family', 'genus']
            for level in taxa_levels:
                unique_taxa = set()
                for mag_id in self.mag_abundance:
                    if mag_id in self.mag_taxonomy:
                        taxon = self.mag_taxonomy[mag_id].get(level, 'Unclassified')
                        if taxon != 'Unclassified':
                            unique_taxa.add(taxon)
                
                self.statistics[f'unique_{level}_count'] = len(unique_taxa)
    
    def get_abundance_dataframe(self):
        """Convert MAG abundance data to pandas DataFrame"""
        if not self.mag_abundance:
            return pd.DataFrame()
        
        df = pd.DataFrame.from_dict(self.mag_abundance, orient='index')
        df.reset_index(drop=True, inplace=True)
        
        # Sort by relative abundance if available
        if 'relative_abundance' in df.columns:
            df = df.sort_values('relative_abundance', ascending=False)
        elif 'weighted_mean_coverage' in df.columns:
            df = df.sort_values('weighted_mean_coverage', ascending=False)
        
        return df

def create_visualization_plots(calculator, output_dir, sample_name):
    """Create visualization plots for MAG abundance"""
    logging.info("Creating MAG abundance visualization plots")
    
    # Create plots directory
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    df = calculator.get_abundance_dataframe()
    
    if df.empty:
        logging.warning("No MAG data available for plotting")
        return
    
    # 1. MAG abundance overview
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Relative abundance distribution
    if 'relative_abundance' in df.columns:
        ax = axes[0, 0]
        values = df['relative_abundance']
        values = values[values > 0]
        
        if len(values) > 0:
            ax.hist(values, bins=30, alpha=0.7, edgecolor='black')
            ax.set_xlabel('Relative Abundance (%)')
            ax.set_ylabel('Frequency')
            ax.set_title(f'{sample_name}: MAG Relative Abundance Distribution')
            ax.set_yscale('log')
            ax.grid(True, alpha=0.3)
    
    # Genome size vs abundance
    if 'total_length' in df.columns and 'relative_abundance' in df.columns:
        ax = axes[0, 1]
        scatter = ax.scatter(df['total_length'] / 1e6, df['relative_abundance'],
                           c=df.get('quality_score', np.ones(len(df))),
                           cmap='viridis', alpha=0.7, s=50)
        ax.set_xlabel('Genome Size (Mb)')
        ax.set_ylabel('Relative Abundance (%)')
        ax.set_title(f'{sample_name}: Genome Size vs Abundance')
        ax.set_yscale('log')
        plt.colorbar(scatter, ax=ax, label='Quality Score')
        ax.grid(True, alpha=0.3)
    
    # Quality distribution
    if 'completeness' in df.columns and 'contamination' in df.columns:
        ax = axes[1, 0]
        scatter = ax.scatter(df['completeness'], df['contamination'],
                           c=df.get('relative_abundance', np.ones(len(df))),
                           cmap='plasma', alpha=0.7, s=50)
        ax.set_xlabel('Completeness (%)')
        ax.set_ylabel('Contamination (%)')
        ax.set_title(f'{sample_name}: MAG Quality Distribution')
        
        # Add quality thresholds
        ax.axhline(y=calculator.max_contamination, color='red', linestyle='--', alpha=0.5)
        ax.axvline(x=calculator.min_completeness, color='red', linestyle='--', alpha=0.5)
        
        plt.colorbar(scatter, ax=ax, label='Relative Abundance (%)')
        ax.grid(True, alpha=0.3)
    
    # MAG count by quality category
    ax = axes[1, 1]
    if hasattr(calculator, 'quality_categories'):
        categories = ['High', 'Medium', 'Low']
        counts = [
            len(calculator.quality_categories.get('high', [])),
            len(calculator.quality_categories.get('medium', [])),
            len(calculator.quality_categories.get('low', []))
        ]
        colors = ['green', 'orange', 'red']
        
        bars = ax.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
        ax.set_ylabel('Number of MAGs')
        ax.set_title(f'{sample_name}: MAG Quality Categories')
        ax.grid(True, alpha=0.3)
        
        # Add counts on bars
        for bar, count in zip(bars, counts):
            if count > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts) * 0.01,
                       str(count), ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(plots_dir / f"{sample_name}_mag_abundance_overview.png",
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Top abundant MAGs
    if 'relative_abundance' in df.columns:
        top_mags = df.head(20)
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Top MAGs by abundance
        y_pos = np.arange(len(top_mags))
        ax1.barh(y_pos, top_mags['relative_abundance'], alpha=0.7)
        
        # Create labels (MAG ID or genus if available)
        labels = []
        for _, row in top_mags.iterrows():
            if 'genus' in row and row['genus'] != 'Unclassified':
                labels.append(f"{row['mag_id'][:10]}... ({row['genus']})")
            else:
                labels.append(f"{row['mag_id'][:15]}...")
        
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(labels)
        ax1.set_xlabel('Relative Abundance (%)')
        ax1.set_title(f'{sample_name}: Top 20 MAGs by Abundance')
        ax1.grid(True, alpha=0.3)
        
        # Quality vs abundance for top MAGs
        if 'quality_score' in top_mags.columns:
            colors = ['green' if q >= 70 else 'orange' if q >= 50 else 'red' 
                     for q in top_mags['quality_score']]
            
            ax2.barh(y_pos, top_mags['quality_score'], color=colors, alpha=0.7)
            ax2.set_yticks(y_pos)
            ax2.set_yticklabels(labels)
            ax2.set_xlabel('Quality Score')
            ax2.set_title(f'{sample_name}: Quality Scores for Top MAGs')
            ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plots_dir / f"{sample_name}_top_abundant_mags.png",
                    dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Taxonomic composition (if taxonomy available)
    if calculator.mag_taxonomy:
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        
        # Phylum-level composition
        phylum_abundance = defaultdict(float)
        for _, row in df.iterrows():
            mag_id = row['mag_id']
            abundance = row.get('relative_abundance', 0)
            
            if mag_id in calculator.mag_taxonomy:
                phylum = calculator.mag_taxonomy[mag_id].get('phylum', 'Unclassified')
                phylum_abundance[phylum] += abundance
        
        # Plot top phyla
        sorted_phyla = sorted(phylum_abundance.items(), key=lambda x: x[1], reverse=True)
        top_phyla = sorted_phyla[:10]
        
        if top_phyla:
            phyla, abundances = zip(*top_phyla)
            
            ax = axes[0]
            wedges, texts, autotexts = ax.pie(abundances, labels=phyla, autopct='%1.1f%%')
            ax.set_title(f'{sample_name}: Taxonomic Composition (Phylum)')
        
        # Genus-level for top MAGs
        if len(df) > 0:
            top_10_mags = df.head(10)
            genera = []
            
            for _, row in top_10_mags.iterrows():
                mag_id = row['mag_id']
                if mag_id in calculator.mag_taxonomy:
                    genus = calculator.mag_taxonomy[mag_id].get('genus', 'Unclassified')
                    genera.append(genus[:20])  # Truncate long names
                else:
                    genera.append('Unclassified')
            
            ax = axes[1]
            bars = ax.bar(range(len(genera)), top_10_mags['relative_abundance'])
            ax.set_xticks(range(len(genera)))
            ax.set_xticklabels(genera, rotation=45, ha='right')
            ax.set_ylabel('Relative Abundance (%)')
            ax.set_title(f'{sample_name}: Top 10 MAGs by Genus')
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plots_dir / f"{sample_name}_mag_taxonomy.png",
                    dpi=300, bbox_inches='tight')
        plt.close()
    
    logging.info(f"Plots saved to {plots_dir}")

def save_results(calculator, output_dir, sample_name):
    """Save MAG abundance results"""
    logging.info("Saving MAG abundance results")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save main MAG abundance table
    df = calculator.get_abundance_dataframe()
    abundance_file = output_path / f"{sample_name}_mag_abundance.tsv"
    df.to_csv(abundance_file, sep='\t', index=False)
    
    # Save high-quality MAGs only
    if 'quality_score' in df.columns:
        high_quality = df[df['quality_score'] >= 70]
        if not high_quality.empty:
            hq_file = output_path / f"{sample_name}_high_quality_mag_abundance.tsv"
            high_quality.to_csv(hq_file, sep='\t', index=False)
    
    # Save statistics
    stats_file = output_path / f"{sample_name}_mag_abundance_stats.json"
    with open(stats_file, 'w') as f:
        json.dump(calculator.statistics, f, indent=2)
    
    # Save abundance matrix for cross-sample analysis
    if 'relative_abundance' in df.columns:
        abundance_matrix = df[['mag_id', 'relative_abundance']].copy()
        abundance_matrix.columns = ['mag_id', f'{sample_name}_relative_abundance']
        
        matrix_file = output_path / f"{sample_name}_mag_abundance_matrix.tsv"
        abundance_matrix.to_csv(matrix_file, sep='\t', index=False)
    
    # Save taxonomy summary if available
    if calculator.mag_taxonomy:
        tax_summary = []
        for mag_id, tax_data in calculator.mag_taxonomy.items():
            if mag_id in calculator.mag_abundance:
                abundance_data = calculator.mag_abundance[mag_id]
                summary_row = {
                    'mag_id': mag_id,
                    'relative_abundance': abundance_data.get('relative_abundance', 0),
                    **tax_data
                }
                tax_summary.append(summary_row)
        
        if tax_summary:
            tax_df = pd.DataFrame(tax_summary)
            tax_file = output_path / f"{sample_name}_mag_taxonomy_abundance.tsv"
            tax_df.to_csv(tax_file, sep='\t', index=False)
    
    logging.info(f"Results saved to {output_path}")
    
    # Print summary
    stats = calculator.statistics
    print(f"\n🧪 MAG Abundance Summary for {sample_name}")
    print("=" * 50)
    print(f"Total MAGs analyzed: {stats.get('total_mags', 0):,}")
    print(f"High quality MAGs: {stats.get('high_quality_mags', 0):,}")
    print(f"Medium quality MAGs: {stats.get('medium_quality_mags', 0):,}")
    print(f"Low quality MAGs: {stats.get('low_quality_mags', 0):,}")
    
    if 'mean_relative_abundance' in stats:
        print(f"Mean relative abundance: {stats['mean_relative_abundance']:.2f}%")
        print(f"Max relative abundance: {stats['max_relative_abundance']:.2f}%")
    
    if 'mean_genome_size' in stats:
        print(f"Mean genome size: {stats['mean_genome_size']/1e6:.2f} Mb")
        print(f"Total genomic content: {stats['total_genomic_content']/1e6:.1f} Mb")

def main():
    """Main function"""
    args = parse_arguments()
    setup_logging(args.log_level)
    
    logging.info("Starting MAG abundance calculation")
    
    # Validate inputs
    if not os.path.exists(args.contig_abundance):
        logging.error(f"Contig abundance file not found: {args.contig_abundance}")
        sys.exit(1)
    
    if not os.path.exists(args.bins_dir):
        logging.error(f"Bins directory not found: {args.bins_dir}")
        sys.exit(1)
    
    # Initialize calculator
    calculator = MAGAbundanceCalculator(
        min_completeness=args.min_completeness,
        max_contamination=args.max_contamination,
        quality_weighting=args.quality_weighting,
        min_contig_length=args.min_contig_length
    )
    
    # Load data
    calculator.load_contig_abundance(args.contig_abundance)
    calculator.load_mag_composition(args.bins_dir)
    calculator.load_mag_quality(args.bin_quality)
    calculator.load_mag_taxonomy(args.taxonomy_file)
    
    # Calculate abundance
    calculator.calculate_mag_abundance(args.normalization)
    
    # Save results
    save_results(calculator, args.output, args.sample_name)
    
    # Generate plots if requested
    if args.plots:
        create_visualization_plots(calculator, args.output, args.sample_name)
    
    logging.info("MAG abundance calculation completed successfully")

if __name__ == "__main__":
    main()

