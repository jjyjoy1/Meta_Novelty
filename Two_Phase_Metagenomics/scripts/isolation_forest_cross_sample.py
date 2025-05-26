#!/usr/bin/env python3
"""
Isolation Forest Cross-Sample Analysis
Multi-feature anomaly detection for identifying unusual samples and features
"""

import sys
import os
import pandas as pd
import numpy as np
import json
import logging
import argparse
from pathlib import Path

from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns

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

def load_and_integrate_features(abundance_matrix_file, vae_embeddings_file, 
                               novelty_data_file, metadata_file, feature_config):
    """Load and integrate multiple feature types for anomaly detection"""
    logging.info("Loading and integrating feature data")
    
    features = {}
    feature_names = []
    
    # Load abundance data
    if feature_config.get('abundance_based', True):
        logging.info("Loading abundance features")
        abundance_df = pd.read_csv(abundance_matrix_file, sep='\t', index_col=0)
        
        # Calculate diversity metrics
        diversity_features = calculate_diversity_features(abundance_df)
        features['diversity'] = diversity_features
        feature_names.extend([f"diversity_{col}" for col in diversity_features.columns])
        
        # Calculate abundance summary statistics
        abundance_stats = calculate_abundance_statistics(abundance_df)
        features['abundance_stats'] = abundance_stats
        feature_names.extend([f"abundance_{col}" for col in abundance_stats.columns])
    
    # Load VAE embeddings
    if feature_config.get('vae_embeddings', True) and os.path.exists(vae_embeddings_file):
        logging.info("Loading VAE embeddings")
        vae_embeddings = np.load(vae_embeddings_file)
        
        # Assume samples are in same order as abundance matrix
        sample_ids = abundance_df.columns if 'abundance_df' in locals() else None
        if sample_ids is not None and len(vae_embeddings) == len(sample_ids):
            vae_df = pd.DataFrame(
                vae_embeddings, 
                index=sample_ids,
                columns=[f"vae_dim_{i}" for i in range(vae_embeddings.shape[1])]
            )
            features['vae'] = vae_df
            feature_names.extend(vae_df.columns.tolist())
        else:
            logging.warning("VAE embeddings dimension mismatch with abundance data")
    
    # Load novelty scores
    if feature_config.get('novelty_scores', True) and os.path.exists(novelty_data_file):
        logging.info("Loading novelty features")
        novelty_features = calculate_novelty_features(novelty_data_file, abundance_df.columns)
        if novelty_features is not None:
            features['novelty'] = novelty_features
            feature_names.extend([f"novelty_{col}" for col in novelty_features.columns])
    
    # Load temporal features if metadata contains temporal information
    if feature_config.get('temporal_features', True):
        logging.info("Calculating temporal features")
        temporal_features = calculate_temporal_features(metadata_file, abundance_df.columns)
        if temporal_features is not None:
            features['temporal'] = temporal_features
            feature_names.extend([f"temporal_{col}" for col in temporal_features.columns])
    
    # Combine all features
    combined_features = []
    common_samples = None
    
    for feature_type, feature_df in features.items():
        if common_samples is None:
            common_samples = set(feature_df.index)
        else:
            common_samples = common_samples.intersection(set(feature_df.index))
    
    common_samples = sorted(list(common_samples))
    logging.info(f"Common samples across all feature types: {len(common_samples)}")
    
    for feature_type, feature_df in features.items():
        feature_subset = feature_df.loc[common_samples]
        combined_features.append(feature_subset)
    
    if combined_features:
        final_features = pd.concat(combined_features, axis=1)
        # Handle any remaining NaN values
        final_features = final_features.fillna(final_features.mean())
    else:
        raise ValueError("No features could be loaded")
    
    logging.info(f"Final feature matrix shape: {final_features.shape}")
    return final_features, feature_names

def calculate_diversity_features(abundance_df):
    """Calculate diversity metrics for each sample"""
    diversity_metrics = {}
    
    for sample in abundance_df.columns:
        abundances = abundance_df[sample]
        non_zero_abundances = abundances[abundances > 0]
        
        # Basic diversity metrics
        observed_otus = len(non_zero_abundances)
        shannon = -np.sum(non_zero_abundances * np.log(non_zero_abundances + 1e-10))
        simpson = 1 - np.sum(non_zero_abundances ** 2)
        
        # Evenness metrics
        pielou = shannon / np.log(observed_otus) if observed_otus > 1 else 0
        
        # Dominance metrics
        berger_parker = non_zero_abundances.max() if len(non_zero_abundances) > 0 else 0
        
        # Rarity metrics
        singletons = np.sum(abundances == abundances[abundances > 0].min()) if len(non_zero_abundances) > 0 else 0
        rare_taxa_fraction = np.sum(abundances <= np.percentile(non_zero_abundances, 10)) / len(abundances) if len(non_zero_abundances) > 0 else 0
        
        diversity_metrics[sample] = {
            'observed_otus': observed_otus,
            'shannon': shannon,
            'simpson': simpson,
            'pielou_evenness': pielou,
            'berger_parker': berger_parker,
            'singletons': singletons,
            'rare_taxa_fraction': rare_taxa_fraction
        }
    
    return pd.DataFrame(diversity_metrics).T

def calculate_abundance_statistics(abundance_df):
    """Calculate abundance distribution statistics for each sample"""
    stats = {}
    
    for sample in abundance_df.columns:
        abundances = abundance_df[sample]
        non_zero_abundances = abundances[abundances > 0]
        
        if len(non_zero_abundances) > 0:
            stats[sample] = {
                'total_abundance': abundances.sum(),
                'mean_abundance': non_zero_abundances.mean(),
                'median_abundance': non_zero_abundances.median(),
                'std_abundance': non_zero_abundances.std(),
                'cv_abundance': non_zero_abundances.std() / non_zero_abundances.mean() if non_zero_abundances.mean() > 0 else 0,
                'abundance_range': non_zero_abundances.max() - non_zero_abundances.min(),
                'abundance_iqr': np.percentile(non_zero_abundances, 75) - np.percentile(non_zero_abundances, 25),
                'skewness': calculate_skewness(non_zero_abundances),
                'kurtosis': calculate_kurtosis(non_zero_abundances)
            }
        else:
            stats[sample] = {col: 0 for col in ['total_abundance', 'mean_abundance', 'median_abundance', 
                                              'std_abundance', 'cv_abundance', 'abundance_range', 
                                              'abundance_iqr', 'skewness', 'kurtosis']}
    
    return pd.DataFrame(stats).T

def calculate_skewness(data):
    """Calculate skewness of data"""
    if len(data) < 3:
        return 0
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    if std == 0:
        return 0
    return np.mean(((data - mean) / std) ** 3)

def calculate_kurtosis(data):
    """Calculate kurtosis of data"""
    if len(data) < 4:
        return 0
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    if std == 0:
        return 0
    return np.mean(((data - mean) / std) ** 4) - 3

def calculate_novelty_features(novelty_file, sample_ids):
    """Calculate novelty-based features for samples"""
    try:
        novelty_df = pd.read_csv(novelty_file, sep='\t')
        
        if 'sample_id' not in novelty_df.columns:
            logging.warning("No sample_id column in novelty data")
            return None
        
        novelty_features = {}
        
        for sample_id in sample_ids:
            sample_novelty = novelty_df[novelty_df['sample_id'] == sample_id]
            
            if len(sample_novelty) > 0:
                novelty_features[sample_id] = {
                    'mean_novelty_score': sample_novelty['combined_novelty_score'].mean(),
                    'max_novelty_score': sample_novelty['combined_novelty_score'].max(),
                    'std_novelty_score': sample_novelty['combined_novelty_score'].std(),
                    'high_novelty_count': (sample_novelty['combined_novelty_score'] > 0.8).sum(),
                    'medium_novelty_count': ((sample_novelty['combined_novelty_score'] > 0.6) & 
                                           (sample_novelty['combined_novelty_score'] <= 0.8)).sum(),
                    'novel_feature_fraction': (sample_novelty['combined_novelty_score'] > 0.5).sum() / len(sample_novelty)
                }
            else:
                novelty_features[sample_id] = {col: 0 for col in ['mean_novelty_score', 'max_novelty_score', 
                                                                'std_novelty_score', 'high_novelty_count', 
                                                                'medium_novelty_count', 'novel_feature_fraction']}
        
        return pd.DataFrame(novelty_features).T
        
    except Exception as e:
        logging.warning(f"Could not load novelty features: {str(e)}")
        return None

def calculate_temporal_features(metadata_file, sample_ids):
    """Calculate temporal features if temporal data is available"""
    try:
        metadata = pd.read_csv(metadata_file, sep='\t')
        
        if 'collection_date' not in metadata.columns:
            logging.info("No temporal data available")
            return None
        
        # Convert collection_date to datetime
        metadata['collection_date'] = pd.to_datetime(metadata['collection_date'])
        
        temporal_features = {}
        
        # Calculate time-based features
        min_date = metadata['collection_date'].min()
        
        for sample_id in sample_ids:
            sample_metadata = metadata[metadata['sample_id'] == sample_id]
            
            if len(sample_metadata) > 0:
                collection_date = sample_metadata['collection_date'].iloc[0]
                days_from_start = (collection_date - min_date).days
                
                # Day of week (0-6)
                day_of_week = collection_date.weekday()
                
                # Month of year (1-12)
                month_of_year = collection_date.month
                
                # Season (encoded as 0-3)
                season = (collection_date.month - 1) // 3
                
                temporal_features[sample_id] = {
                    'days_from_start': days_from_start,
                    'day_of_week': day_of_week,
                    'month_of_year': month_of_year,
                    'season': season,
                    'is_weekend': 1 if day_of_week >= 5 else 0
                }
            else:
                temporal_features[sample_id] = {col: 0 for col in ['days_from_start', 'day_of_week', 
                                                                 'month_of_year', 'season', 'is_weekend']}
        
        return pd.DataFrame(temporal_features).T
        
    except Exception as e:
        logging.warning(f"Could not calculate temporal features: {str(e)}")
        return None

def run_isolation_forest_analysis(features, config):
    """Run Isolation Forest anomaly detection"""
    logging.info("Running Isolation Forest analysis")
    
    # Scale features
    scaler_type = config.get('scaling_method', 'robust')
    if scaler_type == 'standard':
        scaler = StandardScaler()
    elif scaler_type == 'robust':
        scaler = RobustScaler()
    elif scaler_type == 'minmax':
        scaler = MinMaxScaler()
    else:
        scaler = StandardScaler()
    
    features_scaled = scaler.fit_transform(features)
    
    # Initialize Isolation Forest
    isolation_forest = IsolationForest(
        n_estimators=config.get('n_estimators', 200),
        contamination=config.get('contamination', 0.1),
        random_state=config.get('random_seed', 42),
        n_jobs=-1
    )
    
    # Fit and predict
    anomaly_scores = isolation_forest.fit_predict(features_scaled)
    decision_scores = isolation_forest.decision_function(features_scaled)
    
    # Convert to anomaly probability (0-1 scale)
    anomaly_probabilities = 1 - (decision_scores - decision_scores.min()) / (decision_scores.max() - decision_scores.min())
    
    results = {
        'sample_ids': features.index.tolist(),
        'anomaly_labels': anomaly_scores,  # -1 for anomalies, 1 for normal
        'decision_scores': decision_scores,
        'anomaly_probabilities': anomaly_probabilities,
        'is_anomaly': anomaly_scores == -1,
        'scaler': scaler,
        'isolation_forest': isolation_forest
    }
    
    logging.info(f"Identified {np.sum(anomaly_scores == -1)} anomalous samples out of {len(anomaly_scores)}")
    
    return results

def calculate_feature_importance(isolation_forest, features, n_estimators_sample=50):
    """Calculate feature importance based on isolation path lengths"""
    logging.info("Calculating feature importance")
    
    feature_importance = np.zeros(features.shape[1])
    
    # Sample a subset of estimators for efficiency
    n_estimators_total = len(isolation_forest.estimators_)
    estimator_indices = np.random.choice(
        n_estimators_total, 
        size=min(n_estimators_sample, n_estimators_total), 
        replace=False
    )
    
    for idx in estimator_indices:
        estimator = isolation_forest.estimators_[idx]
        
        # Get the features used in this estimator
        features_used = isolation_forest.estimators_features_[idx]
        
        # Calculate feature importance based on tree structure
        if hasattr(estimator, 'tree_'):
            tree = estimator.tree_
            
            # Simple importance: count how often each feature is used for splitting
            for i in range(tree.node_count):
                if tree.children_left[i] != tree.children_right[i]:  # Internal node
                    feature_idx = tree.feature[i]
                    if feature_idx >= 0:  # Valid feature
                        original_feature_idx = features_used[feature_idx]
                        feature_importance[original_feature_idx] += 1
    
    # Normalize importance scores
    if feature_importance.sum() > 0:
        feature_importance = feature_importance / feature_importance.sum()
    
    # Create feature importance dataframe
    importance_df = pd.DataFrame({
        'feature': features.columns,
        'importance': feature_importance
    }).sort_values('importance', ascending=False)
    
    return importance_df

def perform_outlier_analysis(results, features, metadata_file):
    """Perform detailed analysis of detected outliers"""
    logging.info("Performing detailed outlier analysis")
    
    outlier_indices = np.where(results['is_anomaly'])[0]
    outlier_samples = [results['sample_ids'][i] for i in outlier_indices]
    
    outlier_analysis = {
        'outlier_samples': outlier_samples,
        'outlier_count': len(outlier_samples),
        'outlier_percentage': len(outlier_samples) / len(results['sample_ids']) * 100
    }
    
    # Load metadata for outlier characterization
    try:
        metadata = pd.read_csv(metadata_file, sep='\t')
        
        if len(outlier_samples) > 0:
            outlier_metadata = metadata[metadata['sample_id'].isin(outlier_samples)]
            
            # Analyze outlier characteristics
            outlier_characteristics = {}
            
            for column in metadata.columns:
                if column != 'sample_id':
                    if metadata[column].dtype in ['object', 'category']:
                        # Categorical variable
                        outlier_counts = outlier_metadata[column].value_counts()
                        total_counts = metadata[column].value_counts()
                        
                        outlier_characteristics[column] = {
                            'type': 'categorical',
                            'outlier_distribution': outlier_counts.to_dict(),
                            'total_distribution': total_counts.to_dict()
                        }
                    else:
                        # Numerical variable
                        if not outlier_metadata[column].isna().all():
                            outlier_characteristics[column] = {
                                'type': 'numerical',
                                'outlier_mean': outlier_metadata[column].mean(),
                                'outlier_std': outlier_metadata[column].std(),
                                'total_mean': metadata[column].mean(),
                                'total_std': metadata[column].std()
                            }
            
            outlier_analysis['characteristics'] = outlier_characteristics
        
    except Exception as e:
        logging.warning(f"Could not analyze outlier characteristics: {str(e)}")
    
    return outlier_analysis

def create_visualization_plots(results, features, feature_importance, output_dir):
    """Create visualization plots for anomaly detection results"""
    logging.info("Creating visualization plots")
    
    # 1. Anomaly score distribution
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Decision scores histogram
    axes[0, 0].hist(results['decision_scores'], bins=30, alpha=0.7, edgecolor='black')
    axes[0, 0].set_title('Distribution of Decision Scores')
    axes[0, 0].set_xlabel('Decision Score')
    axes[0, 0].set_ylabel('Frequency')
    
    # Anomaly probability distribution
    axes[0, 1].hist(results['anomaly_probabilities'], bins=30, alpha=0.7, edgecolor='black', color='orange')
    axes[0, 1].set_title('Distribution of Anomaly Probabilities')
    axes[0, 1].set_xlabel('Anomaly Probability')
    axes[0, 1].set_ylabel('Frequency')
    
    # Feature importance (top 15)
    top_features = feature_importance.head(15)
    axes[1, 0].barh(range(len(top_features)), top_features['importance'])
    axes[1, 0].set_yticks(range(len(top_features)))
    axes[1, 0].set_yticklabels(top_features['feature'])
    axes[1, 0].set_title('Top 15 Feature Importance')
    axes[1, 0].set_xlabel('Importance Score')
    
    # Anomaly vs normal samples scatter plot (PCA)
    try:
        pca = PCA(n_components=2)
        features_pca = pca.fit_transform(features)
        
        normal_mask = ~results['is_anomaly']
        anomaly_mask = results['is_anomaly']
        
        axes[1, 1].scatter(features_pca[normal_mask, 0], features_pca[normal_mask, 1], 
                          alpha=0.6, label='Normal', s=50)
        axes[1, 1].scatter(features_pca[anomaly_mask, 0], features_pca[anomaly_mask, 1], 
                          alpha=0.8, label='Anomaly', s=100, color='red', marker='^')
        axes[1, 1].set_title(f'PCA Visualization (Var Explained: {pca.explained_variance_ratio_.sum():.3f})')
        axes[1, 1].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.3f})')
        axes[1, 1].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.3f})')
        axes[1, 1].legend()
        
    except Exception as e:
        logging.warning(f"Could not create PCA plot: {str(e)}")
        axes[1, 1].text(0.5, 0.5, 'PCA plot unavailable', ha='center', va='center', transform=axes[1, 1].transAxes)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'isolation_forest_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Feature correlation heatmap for top features
    try:
        top_feature_names = feature_importance.head(20)['feature'].tolist()
        top_feature_data = features[top_feature_names]
        
        correlation_matrix = top_feature_data.corr()
        
        plt.figure(figsize=(12, 10))
        sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0, 
                   square=True, fmt='.2f', cbar_kws={'shrink': 0.8})
        plt.title('Feature Correlation Heatmap (Top 20 Important Features)')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'feature_correlation_heatmap.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        logging.warning(f"Could not create correlation heatmap: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Isolation Forest cross-sample analysis')
    parser.add_argument('--abundance-matrix', required=True,
                        help='Abundance matrix file')
    parser.add_argument('--vae-embeddings', required=True,
                        help='VAE embeddings file (npy format)')
    parser.add_argument('--novelty-data', required=True,
                        help='Pan-novelty analysis results file')
    parser.add_argument('--metadata', required=True,
                        help='Sample metadata file')
    parser.add_argument('--output-anomaly-scores', required=True,
                        help='Output anomaly scores file')
    parser.add_argument('--output-sample-anomalies', required=True,
                        help='Output sample anomaly analysis file')
    parser.add_argument('--output-feature-importance', required=True,
                        help='Output feature importance file')
    parser.add_argument('--output-outlier-analysis', required=True,
                        help='Output outlier analysis JSON file')
    parser.add_argument('--log-file', required=True,
                        help='Log file path')
    
    # Isolation Forest parameters
    parser.add_argument('--n-estimators', type=int, default=200,
                        help='Number of estimators in Isolation Forest')
    parser.add_argument('--contamination', type=float, default=0.1,
                        help='Expected proportion of anomalies')
    parser.add_argument('--scaling-method', choices=['standard', 'robust', 'minmax'], default='robust',
                        help='Feature scaling method')
    parser.add_argument('--random-seed', type=int, default=42,
                        help='Random seed')
    
    # Feature selection
    parser.add_argument('--use-abundance', action='store_true', default=True,
                        help='Include abundance-based features')
    parser.add_argument('--use-vae-embeddings', action='store_true', default=True,
                        help='Include VAE embeddings')
    parser.add_argument('--use-novelty-scores', action='store_true', default=True,
                        help='Include novelty scores')
    parser.add_argument('--use-temporal-features', action='store_true', default=True,
                        help='Include temporal features')
    
    args = parser.parse_args()
    
    # Setup logging
    os.makedirs(os.path.dirname(args.log_file), exist_ok=True)
    setup_logging(args.log_file)
    
    try:
        logging.info("Starting Isolation Forest cross-sample analysis")
        
        # Feature configuration
        feature_config = {
            'abundance_based': args.use_abundance,
            'vae_embeddings': args.use_vae_embeddings,
            'novelty_scores': args.use_novelty_scores,
            'temporal_features': args.use_temporal_features
        }
        
        # Load and integrate features
        features, feature_names = load_and_integrate_features(
            args.abundance_matrix,
            args.vae_embeddings,
            args.novelty_data,
            args.metadata,
            feature_config
        )
        
        # Isolation Forest configuration
        if_config = {
            'n_estimators': args.n_estimators,
            'contamination': args.contamination,
            'scaling_method': args.scaling_method,
            'random_seed': args.random_seed
        }
        
        # Run Isolation Forest analysis
        results = run_isolation_forest_analysis(features, if_config)
        
        # Calculate feature importance
        feature_importance = calculate_feature_importance(results['isolation_forest'], features)
        
        # Perform outlier analysis
        outlier_analysis = perform_outlier_analysis(results, features, args.metadata)
        
        # Create output directories
        for output_file in [args.output_anomaly_scores, args.output_sample_anomalies, 
                           args.output_feature_importance, args.output_outlier_analysis]:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Save anomaly scores
        anomaly_df = pd.DataFrame({
            'sample_id': results['sample_ids'],
            'anomaly_label': results['anomaly_labels'],
            'decision_score': results['decision_scores'],
            'anomaly_probability': results['anomaly_probabilities'],
            'is_anomaly': results['is_anomaly']
        })
        anomaly_df.to_csv(args.output_anomaly_scores, sep='\t', index=False)
        
        # Save sample anomaly analysis
        sample_anomaly_df = anomaly_df[anomaly_df['is_anomaly']].copy()
        sample_anomaly_df.to_csv(args.output_sample_anomalies, sep='\t', index=False)
        
        # Save feature importance
        feature_importance.to_csv(args.output_feature_importance, sep='\t', index=False)
        
        # Save outlier analysis
        with open(args.output_outlier_analysis, 'w') as f:
            json.dump(outlier_analysis, f, indent=2, default=str)
        
        # Create visualizations
        output_dir = os.path.dirname(args.output_anomaly_scores)
        create_visualization_plots(results, features, feature_importance, output_dir)
        
        logging.info("Isolation Forest analysis completed successfully")
        
        # Print summary
        print("\n" + "="*80)
        print("ISOLATION FOREST ANALYSIS SUMMARY")
        print("="*80)
        print(f"Total samples analyzed: {len(results['sample_ids'])}")
        print(f"Features used: {features.shape[1]}")
        print(f"Anomalous samples detected: {outlier_analysis['outlier_count']}")
        print(f"Anomaly percentage: {outlier_analysis['outlier_percentage']:.2f}%")
        print(f"Contamination threshold: {args.contamination}")
        print(f"Top 3 important features:")
        for i, (_, row) in enumerate(feature_importance.head(3).iterrows()):
            print(f"  {i+1}. {row['feature']}: {row['importance']:.4f}")
        
        if outlier_analysis['outlier_count'] > 0:
            print(f"\nAnomalous samples:")
            for sample in outlier_analysis['outlier_samples'][:10]:  # Show first 10
                print(f"  - {sample}")
            if outlier_analysis['outlier_count'] > 10:
                print(f"  ... and {outlier_analysis['outlier_count'] - 10} more")
        
        print("="*80)
        
    except Exception as e:
        logging.error(f"Error in Isolation Forest analysis: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()

