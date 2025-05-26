#!/usr/bin/env python3
"""
Pan-Novelty Analysis Across Samples
Comprehensive analysis of novelty patterns across multiple samples
"""

import sys
import os
import pandas as pd
import numpy as np
import json
import logging
import argparse
from pathlib import Path
from collections import defaultdict

from scipy import stats
from scipy.stats import mannwhitneyu, ttest_ind, wilcoxon
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

def setup_logging(log_file):
    """Setup logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def load_novelty_data(novelty_files):
    """Load novelty data from all samples"""
    logging.info("Loading novelty data from all samples")
    
    all_novelty_data = []
    
    for novelty_file in novelty_files:
        if os.path.exists(novelty_file):
            try:
                # Extract sample name from file path
                sample_name = extract_sample_name_from_path(novelty_file)
                
                df = pd.read_csv(novelty_file, sep='\t')
                df['sample_id'] = sample_name
                all_novelty_data.append(df)
                
                logging.info(f"Loaded novelty data for {sample_name}: {len(df)} features")
                
            except Exception as e:
                logging.error(f"Error loading {novelty_file}: {str(e)}")
        else:
            logging.warning(f"Novelty file not found: {novelty_file}")
    
    if not all_novelty_data:
        raise ValueError("No novelty data could be loaded")
    
    # Combine all novelty data
    combined_df = pd.concat(all_novelty_data, ignore_index=True)
    logging.info(f"Combined novelty data: {len(combined_df)} total entries from {len(all_novelty_data)} samples")
    
    return combined_df

def extract_sample_name_from_path(file_path):
    """Extract sample name from file path"""
    # Assume structure: .../sample_name/novelty/combined_novelty_results.tsv
    path_parts = Path(file_path).parts
    for i, part in enumerate(path_parts):
        if part == 'novelty' and i > 0:
            return path_parts[i-1]
    
    # Fallback: use parent directory name
    return Path(file_path).parent.parent.name

def analyze_pan_novelty_distribution(novelty_df, thresholds):
    """Analyze the distribution of novelty across all samples"""
    logging.info("Analyzing pan-novelty distribution")
    
    # Overall novelty statistics
    novelty_stats = {
        'total_features': len(novelty_df),
        'unique_features': len(novelty_df['contig_id'].unique()),
        'total_samples': len(novelty_df['sample_id'].unique()),
        'mean_novelty_score': novelty_df['combined_novelty_score'].mean(),
        'median_novelty_score': novelty_df['combined_novelty_score'].median(),
        'std_novelty_score': novelty_df['combined_novelty_score'].std(),
        'min_novelty_score': novelty_df['combined_novelty_score'].min(),
        'max_novelty_score': novelty_df['combined_novelty_score'].max()
    }
    
    # Novelty categories
    low_threshold = thresholds['low']
    medium_threshold = thresholds['medium']
    high_threshold = thresholds['high']
    
    novelty_categories = pd.cut(
        novelty_df['combined_novelty_score'],
        bins=[0, low_threshold, medium_threshold, high_threshold, 1.0],
        labels=['Low', 'Medium', 'High', 'Very High'],
        include_lowest=True
    )
    
    category_counts = novelty_categories.value_counts()
    novelty_stats['category_distribution'] = category_counts.to_dict()
    novelty_stats['category_percentages'] = (category_counts / len(novelty_df) * 100).to_dict()
    
    return novelty_stats

def identify_core_and_rare_novelty(novelty_df, abundance_df, core_threshold, rare_threshold):
    """Identify core and rare novel features"""
    logging.info("Identifying core and rare novelty features")
    
    # Calculate prevalence of each feature across samples
    feature_sample_counts = novelty_df.groupby('contig_id')['sample_id'].nunique()
    total_samples = len(novelty_df['sample_id'].unique())
    feature_prevalence = feature_sample_counts / total_samples
    
    # Get mean novelty score for each feature
    feature_novelty = novelty_df.groupby('contig_id')['combined_novelty_score'].agg(['mean', 'std', 'max']).fillna(0)
    
    # Combine with abundance data if available
    if abundance_df is not None and not abundance_df.empty:
        # Calculate mean abundance across samples for each feature
        feature_abundance = abundance_df.mean(axis=1)
        feature_analysis = pd.concat([feature_novelty, feature_prevalence.rename('prevalence')], axis=1)
        feature_analysis = feature_analysis.join(feature_abundance.rename('mean_abundance'), how='left')
    else:
        feature_analysis = pd.concat([feature_novelty, feature_prevalence.rename('prevalence')], axis=1)
        feature_analysis['mean_abundance'] = 0
    
    feature_analysis = feature_analysis.fillna(0)
    
    # Classify features
    core_novelty = feature_analysis[
        (feature_analysis['prevalence'] >= core_threshold) & 
        (feature_analysis['mean'] >= 0.5)  # Mean novelty score above 0.5
    ].copy()
    
    rare_novelty = feature_analysis[
        (feature_analysis['prevalence'] <= rare_threshold) & 
        (feature_analysis['mean'] >= 0.7)  # Higher novelty threshold for rare features
    ].copy()
    
    # Sort by novelty score
    core_novelty = core_novelty.sort_values('mean', ascending=False)
    rare_novelty = rare_novelty.sort_values('mean', ascending=False)
    
    logging.info(f"Identified {len(core_novelty)} core novel features")
    logging.info(f"Identified {len(rare_novelty)} rare novel features")
    
    return core_novelty, rare_novelty, feature_analysis

def detect_novelty_emergence_patterns(novelty_df, metadata, emergence_config):
    """Detect patterns of novelty emergence across samples"""
    logging.info("Detecting novelty emergence patterns")
    
    # Merge with metadata to get temporal information
    if 'collection_date' in metadata.columns:
        novelty_with_time = novelty_df.merge(
            metadata[['sample_id', 'collection_date']], 
            on='sample_id', 
            how='left'
        )
        novelty_with_time['collection_date'] = pd.to_datetime(novelty_with_time['collection_date'])
        
        emergence_patterns = detect_temporal_emergence(novelty_with_time, emergence_config)
    else:
        emergence_patterns = {}
        logging.warning("No temporal data available for emergence detection")
    
    # Detect sudden appearance patterns
    sudden_emergence = detect_sudden_emergence(novelty_df, emergence_config)
    emergence_patterns['sudden_emergence'] = sudden_emergence
    
    # Detect co-emergence patterns
    co_emergence = detect_co_emergence_patterns(novelty_df, emergence_config)
    emergence_patterns['co_emergence'] = co_emergence
    
    return emergence_patterns

def detect_temporal_emergence(novelty_with_time, emergence_config):
    """Detect temporal patterns in novelty emergence"""
    temporal_patterns = {}
    
    # Group by time periods and analyze novelty trends
    novelty_with_time = novelty_with_time.sort_values('collection_date')
    
    # Create time bins (e.g., weekly, monthly)
    time_bins = pd.cut(
        novelty_with_time['collection_date'],
        bins=10,  # Divide into 10 time periods
        labels=False
    )
    novelty_with_time['time_bin'] = time_bins
    
    # Analyze novelty trends over time
    time_trends = novelty_with_time.groupby('time_bin').agg({
        'combined_novelty_score': ['mean', 'std', 'count'],
        'contig_id': 'nunique'
    }).round(4)
    
    # Flatten column names
    time_trends.columns = ['_'.join(col).strip() for col in time_trends.columns]
    
    # Detect significant trends
    if len(time_trends) > 3:
        # Test for correlation between time and novelty
        time_corr, time_p = stats.pearsonr(
            time_trends.index, 
            time_trends['combined_novelty_score_mean']
        )
        
        temporal_patterns['time_correlation'] = {
            'correlation': time_corr,
            'p_value': time_p,
            'significant': time_p < 0.05,
            'trend_direction': 'increasing' if time_corr > 0 else 'decreasing'
        }
    
    temporal_patterns['time_trends'] = time_trends
    
    return temporal_patterns

def detect_sudden_emergence(novelty_df, emergence_config):
    """Detect features that suddenly appear with high novelty"""
    min_fold_change = emergence_config.get('min_fold_change', 2.0)
    min_samples_present = emergence_config.get('min_samples_present', 2)
    
    sudden_emergence = []
    
    # Group by feature and analyze appearance patterns
    for feature_id, feature_data in novelty_df.groupby('contig_id'):
        if len(feature_data) >= min_samples_present:
            # Check if feature appears suddenly with high novelty
            high_novelty_samples = feature_data[feature_data['combined_novelty_score'] > 0.7]
            
            if len(high_novelty_samples) >= min_samples_present:
                sudden_emergence.append({
                    'contig_id': feature_id,
                    'samples_present': len(feature_data),
                    'high_novelty_samples': len(high_novelty_samples),
                    'max_novelty_score': feature_data['combined_novelty_score'].max(),
                    'mean_novelty_score': feature_data['combined_novelty_score'].mean(),
                    'sample_ids': feature_data['sample_id'].tolist()
                })
    
    sudden_emergence_df = pd.DataFrame(sudden_emergence)
    if not sudden_emergence_df.empty:
        sudden_emergence_df = sudden_emergence_df.sort_values('max_novelty_score', ascending=False)
    
    return sudden_emergence_df

def detect_co_emergence_patterns(novelty_df, emergence_config):
    """Detect features that co-emerge together"""
    # Find features that appear together in samples
    sample_features = novelty_df.groupby('sample_id')['contig_id'].apply(set).to_dict()
    
    co_emergence_patterns = []
    
    # Calculate Jaccard similarity between feature sets
    feature_ids = novelty_df['contig_id'].unique()
    
    for i, feature1 in enumerate(feature_ids):
        for j, feature2 in enumerate(feature_ids[i+1:], i+1):
            # Find samples where both features appear
            feature1_samples = set(novelty_df[novelty_df['contig_id'] == feature1]['sample_id'])
            feature2_samples = set(novelty_df[novelty_df['contig_id'] == feature2]['sample_id'])
            
            intersection = feature1_samples.intersection(feature2_samples)
            union = feature1_samples.union(feature2_samples)
            
            if len(union) > 0:
                jaccard_similarity = len(intersection) / len(union)
                
                if jaccard_similarity > 0.5 and len(intersection) >= 2:  # Threshold for co-emergence
                    co_emergence_patterns.append({
                        'feature1': feature1,
                        'feature2': feature2,
                        'jaccard_similarity': jaccard_similarity,
                        'shared_samples': len(intersection),
                        'total_samples': len(union),
                        'shared_sample_ids': list(intersection)
                    })
    
    co_emergence_df = pd.DataFrame(co_emergence_patterns)
    if not co_emergence_df.empty:
        co_emergence_df = co_emergence_df.sort_values('jaccard_similarity', ascending=False)
    
    return co_emergence_df

def perform_differential_novelty_analysis(novelty_df, metadata):
    """Perform differential novelty analysis between groups"""
    logging.info("Performing differential novelty analysis")
    
    differential_results = {}
    
    # Test for differences between categorical variables
    categorical_cols = metadata.select_dtypes(include=['object', 'category']).columns
    
    for col in categorical_cols:
        if col in ['sample_id']:
            continue
            
        logging.info(f"Testing differential novelty for {col}")
        
        # Merge novelty data with metadata
        novelty_with_groups = novelty_df.merge(
            metadata[['sample_id', col]], 
            on='sample_id', 
            how='left'
        )
        
        # Get unique groups
        groups = novelty_with_groups[col].dropna().unique()
        
        if len(groups) >= 2:
            differential_results[col] = perform_group_comparison(novelty_with_groups, col, groups)
    
    return differential_results

def perform_group_comparison(novelty_with_groups, group_col, groups):
    """Compare novelty between groups"""
    results = {}
    
    # Overall group comparison
    group_stats = []
    group_data = []
    
    for group in groups:
        group_novelty = novelty_with_groups[novelty_with_groups[group_col] == group]['combined_novelty_score']
        
        group_stats.append({
            'group': group,
            'count': len(group_novelty),
            'mean_novelty': group_novelty.mean(),
            'median_novelty': group_novelty.median(),
            'std_novelty': group_novelty.std()
        })
        
        group_data.append(group_novelty.values)
    
    results['group_statistics'] = pd.DataFrame(group_stats)
    
    # Statistical tests
    if len(groups) == 2:
        # Two-group comparison
        stat, p_value = mannwhitneyu(group_data[0], group_data[1], alternative='two-sided')
        results['statistical_test'] = {
            'test': 'mann_whitney_u',
            'statistic': stat,
            'p_value': p_value,
            'significant': p_value < 0.05
        }
    elif len(groups) > 2:
        # Multi-group comparison
        stat, p_value = stats.kruskal(*group_data)
        results['statistical_test'] = {
            'test': 'kruskal_wallis',
            'statistic': stat,
            'p_value': p_value,
            'significant': p_value < 0.05
        }
    
    return results

def cluster_samples_by_novelty(novelty_df, n_clusters=None):
    """Cluster samples based on their novelty profiles"""
    logging.info("Clustering samples by novelty profiles")
    
    # Create sample-feature novelty matrix
    novelty_matrix = novelty_df.pivot_table(
        index='sample_id',
        columns='contig_id',
        values='combined_novelty_score',
        fill_value=0
    )
    
    # Standardize features
    scaler = StandardScaler()
    novelty_scaled = scaler.fit_transform(novelty_matrix)
    
    # Determine optimal number of clusters if not specified
    if n_clusters is None:
        n_clusters = min(8, max(2, len(novelty_matrix) // 3))
    
    # Perform K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(novelty_scaled)
    
    # Create cluster results
    cluster_results = pd.DataFrame({
        'sample_id': novelty_matrix.index,
        'cluster': cluster_labels
    })
    
    # Calculate cluster characteristics
    cluster_characteristics = []
    for cluster_id in range(n_clusters):
        cluster_samples = novelty_matrix.index[cluster_labels == cluster_id]
        cluster_novelty = novelty_df[novelty_df['sample_id'].isin(cluster_samples)]
        
        cluster_characteristics.append({
            'cluster_id': cluster_id,
            'sample_count': len(cluster_samples),
            'mean_novelty': cluster_novelty['combined_novelty_score'].mean(),
            'std_novelty': cluster_novelty['combined_novelty_score'].std(),
            'unique_features': len(cluster_novelty['contig_id'].unique()),
            'high_novelty_features': (cluster_novelty['combined_novelty_score'] > 0.8).sum()
        })
    
    cluster_char_df = pd.DataFrame(cluster_characteristics)
    
    return cluster_results, cluster_char_df

def create_novelty_visualizations(novelty_df, core_novelty, rare_novelty, 
                                emergence_patterns, output_dir):
    """Create comprehensive novelty visualizations"""
    logging.info("Creating novelty visualizations")
    
    # 1. Novelty distribution plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Overall novelty score distribution
    axes[0, 0].hist(novelty_df['combined_novelty_score'], bins=50, alpha=0.7, edgecolor='black')
    axes[0, 0].set_title('Distribution of Novelty Scores')
    axes[0, 0].set_xlabel('Combined Novelty Score')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].axvline(novelty_df['combined_novelty_score'].mean(), color='red', linestyle='--', label='Mean')
    axes[0, 0].legend()
    
    # Novelty by sample
    sample_novelty = novelty_df.groupby('sample_id')['combined_novelty_score'].mean()
    axes[0, 1].bar(range(len(sample_novelty)), sample_novelty.values)
    axes[0, 1].set_title('Mean Novelty Score by Sample')
    axes[0, 1].set_xlabel('Sample Index')
    axes[0, 1].set_ylabel('Mean Novelty Score')
    axes[0, 1].tick_params(axis='x', rotation=45)
    
    # Core vs rare novelty comparison
    if not core_novelty.empty and not rare_novelty.empty:
        axes[1, 0].scatter(core_novelty['prevalence'], core_novelty['mean'], 
                          alpha=0.6, label='Core Novelty', s=50)
        axes[1, 0].scatter(rare_novelty['prevalence'], rare_novelty['mean'], 
                          alpha=0.6, label='Rare Novelty', s=50, color='red')
        axes[1, 0].set_xlabel('Prevalence Across Samples')
        axes[1, 0].set_ylabel('Mean Novelty Score')
        axes[1, 0].set_title('Core vs Rare Novel Features')
        axes[1, 0].legend()
    
    # Novelty score vs feature count
    feature_counts = novelty_df['contig_id'].value_counts()
    feature_novelty_mean = novelty_df.groupby('contig_id')['combined_novelty_score'].mean()
    
    axes[1, 1].scatter(feature_counts.values, 
                      [feature_novelty_mean[fid] for fid in feature_counts.index],
                      alpha=0.6)
    axes[1, 1].set_xlabel('Feature Frequency Across Samples')
    axes[1, 1].set_ylabel('Mean Novelty Score')
    axes[1, 1].set_title('Novelty vs Feature Frequency')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'novelty_distribution_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Emergence patterns visualization
    if 'sudden_emergence' in emergence_patterns and not emergence_patterns['sudden_emergence'].empty:
        sudden_df = emergence_patterns['sudden_emergence']
        
        plt.figure(figsize=(12, 8))
        
        plt.subplot(2, 1, 1)
        plt.scatter(sudden_df['samples_present'], sudden_df['max_novelty_score'], alpha=0.6)
        plt.xlabel('Number of Samples Present')
        plt.ylabel('Maximum Novelty Score')
        plt.title('Sudden Emergence Pattern Analysis')
        
        plt.subplot(2, 1, 2)
        plt.hist(sudden_df['max_novelty_score'], bins=20, alpha=0.7, edgecolor='black')
        plt.xlabel('Maximum Novelty Score')
        plt.ylabel('Frequency')
        plt.title('Distribution of Maximum Novelty for Emergent Features')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'emergence_patterns.png'), dpi=300, bbox_inches='tight')
        plt.close()

def main():
    parser = argparse.ArgumentParser(description='Pan-novelty analysis across samples')
    parser.add_argument('--novelty-files', nargs='+', required=True,
                        help='List of novelty result files from single-sample analysis')
    parser.add_argument('--abundance-matrix', required=True,
                        help='Abundance matrix file')
    parser.add_argument('--metadata', required=True,
                        help='Sample metadata file')
    parser.add_argument('--output-pan-novelty', required=True,
                        help='Output pan-novelty results file')
    parser.add_argument('--output-emergence-patterns', required=True,
                        help='Output emergence patterns file')
    parser.add_argument('--output-core-novelty', required=True,
                        help='Output core novelty features file')
    parser.add_argument('--output-rare-novelty', required=True,
                        help='Output rare novelty features file')
    parser.add_argument('--output-summary', required=True,
                        help='Output summary JSON file')
    parser.add_argument('--log-file', required=True,
                        help='Log file path')
    
    # Analysis parameters
    parser.add_argument('--min-samples', type=int, default=3,
                        help='Minimum number of samples for analysis')
    parser.add_argument('--low-threshold', type=float, default=0.3,
                        help='Low novelty threshold')
    parser.add_argument('--medium-threshold', type=float, default=0.6,
                        help='Medium novelty threshold')
    parser.add_argument('--high-threshold', type=float, default=0.8,
                        help='High novelty threshold')
    parser.add_argument('--core-threshold', type=float, default=0.8,
                        help='Core novelty prevalence threshold')
    parser.add_argument('--rare-threshold', type=float, default=0.1,
                        help='Rare novelty prevalence threshold')
    parser.add_argument('--min-fold-change', type=float, default=2.0,
                        help='Minimum fold change for emergence detection')
    
    args = parser.parse_args()
    
    # Setup logging
    os.makedirs(os.path.dirname(args.log_file), exist_ok=True)
    setup_logging(args.log_file)
    
    try:
        logging.info("Starting pan-novelty analysis")
        
        # Load novelty data
        novelty_df = load_novelty_data(args.novelty_files)
        
        # Load abundance data
        abundance_df = None
        if os.path.exists(args.abundance_matrix):
            abundance_df = pd.read_csv(args.abundance_matrix, sep='\t', index_col=0)
            logging.info(f"Loaded abundance matrix: {abundance_df.shape}")
        
        # Load metadata
        metadata = pd.read_csv(args.metadata, sep='\t')
        
        # Check minimum samples requirement
        unique_samples = len(novelty_df['sample_id'].unique())
        if unique_samples < args.min_samples:
            logging.warning(f"Only {unique_samples} samples available, minimum required: {args.min_samples}")
        
        # Analysis configuration
        thresholds = {
            'low': args.low_threshold,
            'medium': args.medium_threshold,
            'high': args.high_threshold
        }
        
        emergence_config = {
            'min_fold_change': args.min_fold_change,
            'min_samples_present': 2,
            'statistical_test': 'mannwhitneyu',
            'p_value_threshold': 0.05
        }
        
        # Analyze pan-novelty distribution
        novelty_stats = analyze_pan_novelty_distribution(novelty_df, thresholds)
        
        # Identify core and rare novelty
        core_novelty, rare_novelty, feature_analysis = identify_core_and_rare_novelty(
            novelty_df, abundance_df, args.core_threshold, args.rare_threshold
        )
        
        # Detect emergence patterns
        emergence_patterns = detect_novelty_emergence_patterns(
            novelty_df, metadata, emergence_config
        )
        
        # Perform differential analysis
        differential_results = perform_differential_novelty_analysis(novelty_df, metadata)
        
        # Cluster samples by novelty
        cluster_results, cluster_characteristics = cluster_samples_by_novelty(novelty_df)
        
        # Create output directories
        for output_file in [args.output_pan_novelty, args.output_emergence_patterns,
                           args.output_core_novelty, args.output_rare_novelty, args.output_summary]:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Save results
        logging.info("Saving pan-novelty analysis results")
        
        # Pan-novelty results
        feature_analysis.to_csv(args.output_pan_novelty, sep='\t')
        
        # Emergence patterns
        if 'sudden_emergence' in emergence_patterns and not emergence_patterns['sudden_emergence'].empty:
            emergence_patterns['sudden_emergence'].to_csv(args.output_emergence_patterns, sep='\t', index=False)
        else:
            # Create empty file with headers
            pd.DataFrame(columns=['contig_id', 'samples_present', 'high_novelty_samples', 
                                'max_novelty_score', 'mean_novelty_score']).to_csv(
                args.output_emergence_patterns, sep='\t', index=False)
        
        # Core and rare novelty
        core_novelty.to_csv(args.output_core_novelty, sep='\t')
        rare_novelty.to_csv(args.output_rare_novelty, sep='\t')
        
        # Summary
        summary = {
            'novelty_statistics': novelty_stats,
            'emergence_patterns': {k: v.to_dict() if hasattr(v, 'to_dict') else v 
                                 for k, v in emergence_patterns.items()},
            'differential_analysis': differential_results,
            'cluster_analysis': {
                'cluster_results': cluster_results.to_dict(),
                'cluster_characteristics': cluster_characteristics.to_dict()
            },
            'feature_counts': {
                'total_features': len(feature_analysis),
                'core_novel_features': len(core_novelty),
                'rare_novel_features': len(rare_novelty)
            }
        }
        
        with open(args.output_summary, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        
        # Create visualizations
        output_dir = os.path.dirname(args.output_pan_novelty)
        create_novelty_visualizations(novelty_df, core_novelty, rare_novelty, 
                                    emergence_patterns, output_dir)
        
        logging.info("Pan-novelty analysis completed successfully")
        
        # Print summary
        print("\n" + "="*80)
        print("PAN-NOVELTY ANALYSIS SUMMARY")
        print("="*80)
        print(f"Total samples analyzed: {unique_samples}")
        print(f"Total novel features: {novelty_stats['unique_features']}")
        print(f"Mean novelty score: {novelty_stats['mean_novelty_score']:.4f}")
        print(f"Core novel features: {len(core_novelty)}")
        print(f"Rare novel features: {len(rare_novelty)}")
        
        print(f"\nNovelty distribution:")
        for category, count in novelty_stats['category_distribution'].items():
            percentage = novelty_stats['category_percentages'][category]
            print(f"  {category}: {count} ({percentage:.1f}%)")
        
        if 'sudden_emergence' in emergence_patterns and not emergence_patterns['sudden_emergence'].empty:
            print(f"\nEmergent features detected: {len(emergence_patterns['sudden_emergence'])}")
            
        print("="*80)
        
    except Exception as e:
        logging.error(f"Error in pan-novelty analysis: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()

