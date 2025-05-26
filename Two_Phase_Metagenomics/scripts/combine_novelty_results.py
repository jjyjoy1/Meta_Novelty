#!/usr/bin/env python3
"""
🔗 Combined Novelty Results Integration
Part of the Metagenomics Pipeline - Phase 1

Combines homology-based and ML-based novelty detection results using
weighted scoring with confidence estimation and method agreement analysis.
"""

import argparse
import os
import sys
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
from scipy import stats
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import cohen_kappa_score
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
        description="Combine homology and ML novelty detection results"
    )
    parser.add_argument(
        "--homology-results", "-h",
        required=True,
        help="Path to homology novelty results TSV file"
    )
    parser.add_argument(
        "--ml-results", "-m",
        required=True,
        help="Path to ML novelty results TSV file"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for combined results"
    )
    parser.add_argument(
        "--sample-name",
        default="sample",
        help="Sample name for output files"
    )
    parser.add_argument(
        "--homology-weight",
        type=float,
        default=0.4,
        help="Weight for homology-based scores (default: 0.4)"
    )
    parser.add_argument(
        "--ml-weight",
        type=float,
        default=0.6,
        help="Weight for ML-based scores (default: 0.6)"
    )
    parser.add_argument(
        "--confidence-boost",
        type=float,
        default=0.2,
        help="Confidence boost factor for method agreement (default: 0.2)"
    )
    parser.add_argument(
        "--agreement-threshold",
        type=float,
        default=0.1,
        help="Threshold for method agreement (default: 0.1)"
    )
    parser.add_argument(
        "--high-novelty-threshold",
        type=float,
        default=75.0,
        help="Threshold for high novelty classification (default: 75)"
    )
    parser.add_argument(
        "--medium-novelty-threshold",
        type=float,
        default=50.0,
        help="Threshold for medium novelty classification (default: 50)"
    )
    parser.add_argument(
        "--normalize-scores",
        action="store_true",
        help="Normalize input scores before combining"
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

class NoveltyIntegrator:
    """Main class for integrating novelty detection results"""
    
    def __init__(self, homology_weight=0.4, ml_weight=0.6, 
                 confidence_boost=0.2, agreement_threshold=0.1):
        self.homology_weight = homology_weight
        self.ml_weight = ml_weight
        self.confidence_boost = confidence_boost
        self.agreement_threshold = agreement_threshold
        
        # Ensure weights sum to 1
        total_weight = homology_weight + ml_weight
        if total_weight != 1.0:
            logging.warning(f"Weights don't sum to 1.0 ({total_weight}), normalizing...")
            self.homology_weight = homology_weight / total_weight
            self.ml_weight = ml_weight / total_weight
        
        self.homology_data = None
        self.ml_data = None
        self.combined_data = None
        self.statistics = {}
    
    def load_data(self, homology_file, ml_file):
        """Load homology and ML results"""
        logging.info("Loading novelty detection results")
        
        # Load homology results
        try:
            if homology_file.endswith('.tsv'):
                self.homology_data = pd.read_csv(homology_file, sep='\t')
            else:
                self.homology_data = pd.read_csv(homology_file)
            
            logging.info(f"Loaded homology results: {len(self.homology_data)} sequences")
            
        except Exception as e:
            logging.error(f"Error loading homology results: {e}")
            sys.exit(1)
        
        # Load ML results
        try:
            if ml_file.endswith('.tsv'):
                self.ml_data = pd.read_csv(ml_file, sep='\t')
            else:
                self.ml_data = pd.read_csv(ml_file)
            
            logging.info(f"Loaded ML results: {len(self.ml_data)} sequences")
            
        except Exception as e:
            logging.error(f"Error loading ML results: {e}")
            sys.exit(1)
    
    def prepare_data(self, normalize_scores=False):
        """Prepare and align data from both methods"""
        logging.info("Preparing and aligning data")
        
        # Handle different column names and data structures
        homology_df = self.prepare_homology_data()
        ml_df = self.prepare_ml_data()
        
        # Merge datasets on sequence ID
        merged_df = self.merge_datasets(homology_df, ml_df)
        
        if merged_df.empty:
            logging.error("No overlapping sequences found between datasets")
            sys.exit(1)
        
        logging.info(f"Found {len(merged_df)} overlapping sequences")
        
        # Normalize scores if requested
        if normalize_scores:
            merged_df = self.normalize_scores(merged_df)
        
        self.combined_data = merged_df
        return merged_df
    
    def prepare_homology_data(self):
        """Prepare homology data with standardized columns"""
        df = self.homology_data.copy()
        
        # Standardize column names
        column_mapping = {
            'sequence_id': 'sequence_id',
            'contig_id': 'sequence_id',
            'homology_novelty_score': 'homology_score',
            'novelty_score': 'homology_score',
            'novelty_level': 'homology_level',
            'level': 'homology_level',
            'total_hits': 'total_hits',
            'best_identity': 'best_identity',
            'sequence_length': 'sequence_length',
            'length': 'sequence_length'
        }
        
        # Rename columns if they exist
        for old_name, new_name in column_mapping.items():
            if old_name in df.columns:
                df = df.rename(columns={old_name: new_name})
        
        # Ensure required columns exist
        required_cols = ['sequence_id', 'homology_score']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            logging.error(f"Missing required columns in homology data: {missing_cols}")
            sys.exit(1)
        
        return df
    
    def prepare_ml_data(self):
        """Prepare ML data with standardized columns"""
        df = self.ml_data.copy()
        
        # Handle windowed sequences by aggregating to contig level
        if 'original_contig_id' in df.columns:
            logging.info("Aggregating windowed ML results to contig level")
            df = self.aggregate_windowed_ml_results(df)
        
        # Standardize column names
        column_mapping = {
            'sequence_id': 'sequence_id',
            'contig_id': 'sequence_id',
            'original_contig_id': 'sequence_id',
            'ml_novelty_score': 'ml_score',
            'novelty_score': 'ml_score',
            'mean_novelty_score': 'ml_score',
            'ml_prediction': 'ml_prediction',
            'prediction': 'ml_prediction',
            'sequence_length': 'sequence_length',
            'length': 'sequence_length'
        }
        
        # Rename columns if they exist
        for old_name, new_name in column_mapping.items():
            if old_name in df.columns:
                df = df.rename(columns={old_name: new_name})
        
        # Ensure required columns exist
        required_cols = ['sequence_id', 'ml_score']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            logging.error(f"Missing required columns in ML data: {missing_cols}")
            sys.exit(1)
        
        return df
    
    def aggregate_windowed_ml_results(self, df):
        """Aggregate windowed ML results to contig level"""
        # Group by original contig ID and aggregate
        agg_funcs = {
            'ml_novelty_score': 'mean',
            'anomaly_score': 'mean',
            'sequence_length': 'sum'
        }
        
        # Add other numeric columns for aggregation
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            if col not in agg_funcs and col != 'window_index':
                agg_funcs[col] = 'mean'
        
        aggregated = df.groupby('original_contig_id').agg(agg_funcs).reset_index()
        aggregated = aggregated.rename(columns={'original_contig_id': 'sequence_id'})
        
        # Calculate additional metrics
        window_stats = df.groupby('original_contig_id').agg({
            'ml_novelty_score': ['std', 'max', 'min', 'count'],
            'ml_prediction': lambda x: (x == 'Novel').sum()
        }).reset_index()
        
        # Flatten column names
        window_stats.columns = ['sequence_id', 'ml_score_std', 'ml_score_max', 
                               'ml_score_min', 'total_windows', 'novel_windows']
        
        # Merge with aggregated data
        aggregated = aggregated.merge(window_stats, on='sequence_id', how='left')
        
        # Calculate novelty percentage for windowed sequences
        aggregated['window_novelty_percentage'] = (
            aggregated['novel_windows'] / aggregated['total_windows'] * 100
        )
        
        return aggregated
    
    def merge_datasets(self, homology_df, ml_df):
        """Merge homology and ML datasets"""
        # Merge on sequence ID
        merged = homology_df.merge(ml_df, on='sequence_id', how='inner', 
                                  suffixes=('_homology', '_ml'))
        
        # Handle sequence length conflicts
        if 'sequence_length_homology' in merged.columns and 'sequence_length_ml' in merged.columns:
            # Use homology length as primary (more reliable for full contigs)
            merged['sequence_length'] = merged['sequence_length_homology']
            merged = merged.drop(['sequence_length_homology', 'sequence_length_ml'], axis=1)
        elif 'sequence_length_homology' in merged.columns:
            merged['sequence_length'] = merged['sequence_length_homology']
            merged = merged.drop(['sequence_length_homology'], axis=1)
        elif 'sequence_length_ml' in merged.columns:
            merged['sequence_length'] = merged['sequence_length_ml']
            merged = merged.drop(['sequence_length_ml'], axis=1)
        
        return merged
    
    def normalize_scores(self, df):
        """Normalize novelty scores to 0-100 range"""
        logging.info("Normalizing novelty scores")
        
        scaler = MinMaxScaler(feature_range=(0, 100))
        
        # Normalize homology scores
        if 'homology_score' in df.columns:
            df['homology_score_normalized'] = scaler.fit_transform(
                df[['homology_score']]
            ).flatten()
            df['homology_score'] = df['homology_score_normalized']
            df = df.drop(['homology_score_normalized'], axis=1)
        
        # Normalize ML scores
        if 'ml_score' in df.columns:
            df['ml_score_normalized'] = scaler.fit_transform(
                df[['ml_score']]
            ).flatten()
            df['ml_score'] = df['ml_score_normalized']
            df = df.drop(['ml_score_normalized'], axis=1)
        
        return df
    
    def calculate_combined_scores(self, high_threshold=75.0, medium_threshold=50.0):
        """Calculate combined novelty scores with confidence estimation"""
        logging.info("Calculating combined novelty scores")
        
        df = self.combined_data.copy()
        
        # Calculate method agreement
        df['score_difference'] = abs(df['homology_score'] - df['ml_score'])
        df['method_agreement'] = df['score_difference'] <= (self.agreement_threshold * 100)
        
        # Calculate confidence boost based on agreement
        df['confidence_multiplier'] = np.where(
            df['method_agreement'], 
            1 + self.confidence_boost, 
            1.0
        )
        
        # Calculate weighted average
        df['combined_score_base'] = (
            self.homology_weight * df['homology_score'] + 
            self.ml_weight * df['ml_score']
        )
        
        # Apply confidence boost
        df['combined_novelty_score'] = df['combined_score_base'] * df['confidence_multiplier']
        
        # Cap at 100
        df['combined_novelty_score'] = np.minimum(df['combined_novelty_score'], 100.0)
        
        # Classify combined novelty levels
        df['combined_novelty_level'] = pd.cut(
            df['combined_novelty_score'],
            bins=[-np.inf, medium_threshold, high_threshold, np.inf],
            labels=['Low', 'Medium', 'High']
        )
        
        # Calculate confidence score
        df['confidence_score'] = self.calculate_confidence_scores(df)
        
        # Flag high-confidence novel sequences
        df['high_confidence_novel'] = (
            (df['combined_novelty_score'] >= high_threshold) & 
            (df['confidence_score'] >= 0.7)
        )
        
        self.combined_data = df
        
        # Calculate statistics
        self.calculate_statistics()
        
        return df
    
    def calculate_confidence_scores(self, df):
        """Calculate confidence scores based on multiple factors"""
        # Factors affecting confidence:
        # 1. Method agreement (lower difference = higher confidence)
        # 2. Score magnitude (extreme scores = higher confidence)
        # 3. Data quality indicators
        
        confidence_scores = np.zeros(len(df))
        
        # Agreement-based confidence (0-0.4 points)
        max_diff = 100.0
        agreement_confidence = 0.4 * (1 - df['score_difference'] / max_diff)
        confidence_scores += agreement_confidence
        
        # Magnitude-based confidence (0-0.3 points)
        # High or low scores are more confident
        score_extremity = np.maximum(
            df['combined_score_base'] / 100,  # High scores
            (100 - df['combined_score_base']) / 100  # Low scores (inverted)
        )
        magnitude_confidence = 0.3 * score_extremity
        confidence_scores += magnitude_confidence
        
        # Data quality confidence (0-0.3 points)
        quality_confidence = np.zeros(len(df))
        
        # Homology data quality
        if 'total_hits' in df.columns:
            # More hits = more confident for low novelty
            # Fewer hits = more confident for high novelty
            hit_confidence = np.where(
                df['combined_score_base'] > 50,  # High novelty
                0.15 * (1 - np.minimum(df['total_hits'] / 10, 1)),  # Fewer hits better
                0.15 * np.minimum(df['total_hits'] / 10, 1)  # More hits better
            )
            quality_confidence += hit_confidence
        
        # ML data quality
        if 'total_windows' in df.columns:
            # More windows analyzed = higher confidence
            window_confidence = 0.15 * np.minimum(df['total_windows'] / 5, 1)
            quality_confidence += window_confidence
        elif 'sequence_length' in df.columns:
            # Longer sequences = higher confidence
            length_confidence = 0.15 * np.minimum(df['sequence_length'] / 2000, 1)
            quality_confidence += length_confidence
        
        confidence_scores += quality_confidence
        
        # Normalize to 0-1 range
        confidence_scores = np.clip(confidence_scores, 0, 1)
        
        return confidence_scores
    
    def calculate_statistics(self):
        """Calculate comprehensive statistics"""
        df = self.combined_data
        
        # Basic statistics
        self.statistics = {
            'total_sequences': len(df),
            'homology_mean': df['homology_score'].mean(),
            'homology_std': df['homology_score'].std(),
            'ml_mean': df['ml_score'].mean(),
            'ml_std': df['ml_score'].std(),
            'combined_mean': df['combined_novelty_score'].mean(),
            'combined_std': df['combined_novelty_score'].std(),
        }
        
        # Method agreement statistics
        self.statistics.update({
            'method_agreement_rate': df['method_agreement'].mean() * 100,
            'mean_score_difference': df['score_difference'].mean(),
            'median_score_difference': df['score_difference'].median(),
            'correlation_coefficient': stats.pearsonr(df['homology_score'], df['ml_score'])[0]
        })
        
        # Novelty level distributions
        combined_levels = df['combined_novelty_level'].value_counts()
        for level in ['Low', 'Medium', 'High']:
            self.statistics[f'combined_{level.lower()}_count'] = combined_levels.get(level, 0)
            self.statistics[f'combined_{level.lower()}_percentage'] = (
                combined_levels.get(level, 0) / len(df) * 100
            )
        
        # High-confidence novel sequences
        self.statistics['high_confidence_novel_count'] = df['high_confidence_novel'].sum()
        self.statistics['high_confidence_novel_percentage'] = (
            df['high_confidence_novel'].sum() / len(df) * 100
        )
        
        # Confidence statistics
        self.statistics.update({
            'mean_confidence': df['confidence_score'].mean(),
            'median_confidence': df['confidence_score'].median(),
            'high_confidence_sequences': (df['confidence_score'] >= 0.7).sum()
        })
        
        # Method-specific agreement
        if 'homology_level' in df.columns and 'ml_prediction' in df.columns:
            # Convert to binary classifications for kappa calculation
            homology_binary = df['homology_level'].isin(['High', 'Medium'])
            ml_binary = df['ml_prediction'] == 'Novel'
            
            if len(set(homology_binary)) > 1 and len(set(ml_binary)) > 1:
                kappa = cohen_kappa_score(homology_binary, ml_binary)
                self.statistics['cohen_kappa'] = kappa
        
        logging.info("Combined novelty analysis statistics:")
        logging.info(f"  Total sequences: {self.statistics['total_sequences']}")
        logging.info(f"  Method agreement rate: {self.statistics['method_agreement_rate']:.1f}%")
        logging.info(f"  Mean combined score: {self.statistics['combined_mean']:.1f}")
        logging.info(f"  High-confidence novel: {self.statistics['high_confidence_novel_count']}")

def create_visualization_plots(integrator, output_dir, sample_name):
    """Create comprehensive visualization plots"""
    logging.info("Creating visualization plots")
    
    # Create plots directory
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    df = integrator.combined_data
    
    # 1. Method comparison and integration overview
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Score correlation
    ax = axes[0, 0]
    scatter = ax.scatter(df['homology_score'], df['ml_score'], 
                        c=df['combined_novelty_score'], cmap='viridis', alpha=0.6)
    ax.plot([0, 100], [0, 100], 'r--', alpha=0.5, label='Perfect Agreement')
    ax.set_xlabel('Homology Novelty Score')
    ax.set_ylabel('ML Novelty Score')
    ax.set_title(f'{sample_name}: Method Correlation')
    ax.legend()
    plt.colorbar(scatter, ax=ax, label='Combined Score')
    
    # Score distributions
    ax = axes[0, 1]
    ax.hist([df['homology_score'], df['ml_score'], df['combined_novelty_score']], 
           bins=30, alpha=0.6, label=['Homology', 'ML', 'Combined'])
    ax.set_xlabel('Novelty Score')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{sample_name}: Score Distributions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Method agreement
    ax = axes[0, 2]
    agreement_data = df['method_agreement'].value_counts()
    colors = ['lightcoral', 'lightgreen']
    wedges, texts, autotexts = ax.pie(agreement_data.values, 
                                     labels=['Disagreement', 'Agreement'],
                                     colors=colors, autopct='%1.1f%%')
    ax.set_title(f'{sample_name}: Method Agreement')
    
    # Score difference distribution
    ax = axes[1, 0]
    ax.hist(df['score_difference'], bins=30, alpha=0.7, edgecolor='black')
    ax.axvline(df['score_difference'].mean(), color='red', linestyle='--',
              label=f'Mean: {df["score_difference"].mean():.1f}')
    ax.set_xlabel('Score Difference |Homology - ML|')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{sample_name}: Score Agreement Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Combined novelty levels
    ax = axes[1, 1]
    level_counts = df['combined_novelty_level'].value_counts()
    colors = {'High': 'red', 'Medium': 'orange', 'Low': 'green'}
    plot_colors = [colors.get(level, 'gray') for level in level_counts.index]
    
    bars = ax.bar(level_counts.index, level_counts.values, 
                 color=plot_colors, alpha=0.7, edgecolor='black')
    ax.set_ylabel('Number of Sequences')
    ax.set_title(f'{sample_name}: Combined Novelty Levels')
    ax.grid(True, alpha=0.3)
    
    # Add counts on bars
    for bar, count in zip(bars, level_counts.values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(level_counts) * 0.01,
               str(count), ha='center', va='bottom')
    
    # Confidence vs combined score
    ax = axes[1, 2]
    scatter = ax.scatter(df['combined_novelty_score'], df['confidence_score'],
                        c=df['method_agreement'].astype(int), cmap='RdYlGn', alpha=0.6)
    ax.set_xlabel('Combined Novelty Score')
    ax.set_ylabel('Confidence Score')
    ax.set_title(f'{sample_name}: Confidence vs Novelty')
    plt.colorbar(scatter, ax=ax, label='Agreement (0=No, 1=Yes)')
    
    plt.tight_layout()
    plt.savefig(plots_dir / f"{sample_name}_combined_novelty_overview.png",
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Detailed analysis plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # High-confidence novel sequences
    ax = axes[0, 0]
    novel_mask = df['high_confidence_novel']
    ax.scatter(df.loc[~novel_mask, 'sequence_length'], 
              df.loc[~novel_mask, 'combined_novelty_score'],
              alpha=0.4, s=20, label='Other', color='gray')
    ax.scatter(df.loc[novel_mask, 'sequence_length'], 
              df.loc[novel_mask, 'combined_novelty_score'],
              alpha=0.8, s=30, label='High-Confidence Novel', color='red')
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Combined Novelty Score')
    ax.set_title(f'{sample_name}: High-Confidence Novel Sequences')
    ax.set_xscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Method weights visualization
    ax = axes[0, 1]
    weights = [integrator.homology_weight, integrator.ml_weight]
    labels = ['Homology', 'ML']
    colors = ['skyblue', 'lightgreen']
    
    wedges, texts, autotexts = ax.pie(weights, labels=labels, colors=colors,
                                     autopct='%1.1f%%', startangle=90)
    ax.set_title(f'{sample_name}: Method Weights')
    
    # Confidence distribution
    ax = axes[1, 0]
    ax.hist(df['confidence_score'], bins=30, alpha=0.7, edgecolor='black')
    ax.axvline(df['confidence_score'].mean(), color='red', linestyle='--',
              label=f'Mean: {df["confidence_score"].mean():.2f}')
    ax.axvline(0.7, color='orange', linestyle='--', label='High Confidence Threshold')
    ax.set_xlabel('Confidence Score')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{sample_name}: Confidence Score Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Score improvement analysis
    ax = axes[1, 1]
    max_individual = np.maximum(df['homology_score'], df['ml_score'])
    improvement = df['combined_novelty_score'] - max_individual
    
    ax.hist(improvement, bins=30, alpha=0.7, edgecolor='black')
    ax.axvline(0, color='red', linestyle='-', alpha=0.5, label='No Improvement')
    ax.axvline(improvement.mean(), color='orange', linestyle='--',
              label=f'Mean: {improvement.mean():.1f}')
    ax.set_xlabel('Score Improvement (Combined - Best Individual)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{sample_name}: Integration Benefit')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(plots_dir / f"{sample_name}_combined_novelty_detailed.png",
                dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Plots saved to {plots_dir}")

def save_results(integrator, output_dir, sample_name):
    """Save integrated results"""
    logging.info("Saving integrated results")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save main combined results
    main_results = integrator.combined_data.copy()
    
    # Select and order columns for output
    output_columns = [
        'sequence_id', 'sequence_length',
        'homology_score', 'ml_score', 'combined_novelty_score',
        'combined_novelty_level', 'confidence_score',
        'method_agreement', 'high_confidence_novel'
    ]
    
    # Add optional columns if they exist
    optional_columns = [
        'homology_level', 'ml_prediction', 'total_hits', 'best_identity',
        'total_windows', 'window_novelty_percentage', 'score_difference'
    ]
    
    for col in optional_columns:
        if col in main_results.columns:
            output_columns.append(col)
    
    # Filter to existing columns
    existing_columns = [col for col in output_columns if col in main_results.columns]
    results_output = main_results[existing_columns]
    
    # Sort by combined novelty score (descending)
    results_output = results_output.sort_values('combined_novelty_score', ascending=False)
    
    # Save main results
    results_file = output_path / f"{sample_name}_combined_novelty.tsv"
    results_output.to_csv(results_file, sep='\t', index=False)
    
    # Save high-confidence novel sequences
    high_conf_novel = results_output[results_output['high_confidence_novel'] == True]
    if not high_conf_novel.empty:
        novel_file = output_path / f"{sample_name}_high_confidence_novel.tsv"
        high_conf_novel.to_csv(novel_file, sep='\t', index=False)
    
    # Save statistics
    stats_file = output_path / f"{sample_name}_combined_novelty_stats.json"
    with open(stats_file, 'w') as f:
        json.dump(integrator.statistics, f, indent=2)
    
    # Save summary report
    summary_file = output_path / f"{sample_name}_novelty_integration_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"🔗 Combined Novelty Detection Summary - {sample_name}\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Integration Parameters:\n")
        f.write(f"  Homology weight: {integrator.homology_weight:.2f}\n")
        f.write(f"  ML weight: {integrator.ml_weight:.2f}\n")
        f.write(f"  Confidence boost: {integrator.confidence_boost:.2f}\n")
        f.write(f"  Agreement threshold: {integrator.agreement_threshold:.2f}\n\n")
        
        f.write("Dataset Statistics:\n")
        f.write(f"  Total sequences: {integrator.statistics['total_sequences']:,}\n")
        f.write(f"  Method agreement rate: {integrator.statistics['method_agreement_rate']:.1f}%\n")
        f.write(f"  Score correlation: {integrator.statistics.get('correlation_coefficient', 0):.3f}\n")
        f.write(f"  Mean score difference: {integrator.statistics['mean_score_difference']:.1f}\n\n")
        
        f.write("Novelty Level Distribution:\n")
        f.write(f"  High novelty: {integrator.statistics['combined_high_count']:,} "
               f"({integrator.statistics['combined_high_percentage']:.1f}%)\n")
        f.write(f"  Medium novelty: {integrator.statistics['combined_medium_count']:,} "
               f"({integrator.statistics['combined_medium_percentage']:.1f}%)\n")
        f.write(f"  Low novelty: {integrator.statistics['combined_low_count']:,} "
               f"({integrator.statistics['combined_low_percentage']:.1f}%)\n\n")
        
        f.write("Confidence Analysis:\n")
        f.write(f"  Mean confidence: {integrator.statistics['mean_confidence']:.3f}\n")
        f.write(f"  High-confidence sequences: {integrator.statistics['high_confidence_sequences']:,}\n")
        f.write(f"  High-confidence novel: {integrator.statistics['high_confidence_novel_count']:,} "
               f"({integrator.statistics['high_confidence_novel_percentage']:.1f}%)\n")
        
        if 'cohen_kappa' in integrator.statistics:
            f.write(f"  Cohen's kappa: {integrator.statistics['cohen_kappa']:.3f}\n")
    
    logging.info(f"Results saved to {output_path}")
    
    # Print summary to console
    stats = integrator.statistics
    print(f"\n🔗 Combined Novelty Detection Summary for {sample_name}")
    print("=" * 60)
    print(f"Total sequences analyzed: {stats['total_sequences']:,}")
    print(f"Method agreement rate: {stats['method_agreement_rate']:.1f}%")
    print(f"Score correlation: {stats.get('correlation_coefficient', 0):.3f}")
    print(f"Mean combined novelty score: {stats['combined_mean']:.1f}")
    print(f"High-confidence novel sequences: {stats['high_confidence_novel_count']:,} "
          f"({stats['high_confidence_novel_percentage']:.1f}%)")
    print(f"High novelty sequences: {stats['combined_high_count']:,}")
    print(f"Medium novelty sequences: {stats['combined_medium_count']:,}")
    print(f"Low novelty sequences: {stats['combined_low_count']:,}")

def main():
    """Main function"""
    args = parse_arguments()
    setup_logging(args.log_level)
    
    logging.info("Starting combined novelty analysis")
    
    # Validate inputs
    for file_path in [args.homology_results, args.ml_results]:
        if not os.path.exists(file_path):
            logging.error(f"Input file not found: {file_path}")
            sys.exit(1)
    
    # Initialize integrator
    integrator = NoveltyIntegrator(
        homology_weight=args.homology_weight,
        ml_weight=args.ml_weight,
        confidence_boost=args.confidence_boost,
        agreement_threshold=args.agreement_threshold
    )
    
    # Load and process data
    integrator.load_data(args.homology_results, args.ml_results)
    integrator.prepare_data(normalize_scores=args.normalize_scores)
    
    # Calculate combined scores
    combined_results = integrator.calculate_combined_scores(
        high_threshold=args.high_novelty_threshold,
        medium_threshold=args.medium_novelty_threshold
    )
    
    # Save results
    save_results(integrator, args.output, args.sample_name)
    
    # Generate plots if requested
    if args.plots:
        create_visualization_plots(integrator, args.output, args.sample_name)
    
    logging.info("Combined novelty analysis completed successfully")

if __name__ == "__main__":
    main()


