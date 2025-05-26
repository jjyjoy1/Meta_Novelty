#!/usr/bin/env python3
"""
Temporal Analysis of Microbiome Data
Comprehensive time-series analysis for longitudinal microbiome studies
"""

import sys
import os
import pandas as pd
import numpy as np
import json
import logging
import argparse
from pathlib import Path
from datetime import datetime, timedelta

from scipy import stats
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

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

def load_and_prepare_temporal_data(abundance_matrix, novelty_data, metadata, config):
    """Load and prepare data for temporal analysis"""
    logging.info("Loading and preparing temporal data")
    
    # Load abundance data
    abundance_df = pd.read_csv(abundance_matrix, sep='\t', index_col=0)
    
    # Load novelty data if available
    novelty_df = None
    if os.path.exists(novelty_data):
        novelty_df = pd.read_csv(novelty_data, sep='\t', index_col=0)
        logging.info(f"Loaded novelty data: {novelty_df.shape}")
    
    # Load metadata
    metadata_df = pd.read_csv(metadata, sep='\t')
    
    # Check if temporal data is available
    if 'collection_date' not in metadata_df.columns:
        raise ValueError("No 'collection_date' column found in metadata")
    
    # Convert collection_date to datetime
    metadata_df['collection_date'] = pd.to_datetime(metadata_df['collection_date'])
    
    # Filter samples with valid dates
    valid_dates = metadata_df['collection_date'].notna()
    metadata_df = metadata_df[valid_dates]
    
    # Filter abundance data to include only samples with temporal information
    common_samples = set(abundance_df.columns).intersection(set(metadata_df['sample_id']))
    abundance_filtered = abundance_df[sorted(list(common_samples))]
    
    # Sort by collection date
    metadata_sorted = metadata_df.sort_values('collection_date')
    sample_order = metadata_sorted['sample_id'].tolist()
    abundance_sorted = abundance_filtered[sample_order]
    
    logging.info(f"Temporal dataset: {len(common_samples)} samples with dates")
    logging.info(f"Date range: {metadata_sorted['collection_date'].min()} to {metadata_sorted['collection_date'].max()}")
    
    return abundance_sorted, novelty_df, metadata_sorted

def calculate_time_intervals(metadata_df, time_unit='days'):
    """Calculate time intervals from the first sample"""
    logging.info(f"Calculating time intervals in {time_unit}")
    
    first_date = metadata_df['collection_date'].min()
    
    if time_unit == 'days':
        time_intervals = (metadata_df['collection_date'] - first_date).dt.days
    elif time_unit == 'weeks':
        time_intervals = (metadata_df['collection_date'] - first_date).dt.days / 7
    elif time_unit == 'months':
        # Approximate months as 30.44 days
        time_intervals = (metadata_df['collection_date'] - first_date).dt.days / 30.44
    else:
        raise ValueError(f"Unsupported time unit: {time_unit}")
    
    metadata_df = metadata_df.copy()
    metadata_df['time_interval'] = time_intervals
    
    return metadata_df

def detect_temporal_trends(abundance_data, metadata_df, config):
    """Detect temporal trends in microbial abundance"""
    logging.info("Detecting temporal trends")
    
    trend_method = config.get('method', 'mann_kendall')
    significance_threshold = config.get('significance_threshold', 0.05)
    
    trend_results = []
    time_intervals = metadata_df['time_interval'].values
    
    for feature in abundance_data.index:
        abundances = abundance_data.loc[feature].values
        
        # Remove zeros for trend analysis
        non_zero_mask = abundances > 0
        if np.sum(non_zero_mask) < 3:  # Need at least 3 points
            continue
        
        filtered_times = time_intervals[non_zero_mask]
        filtered_abundances = abundances[non_zero_mask]
        
        if trend_method == 'mann_kendall':
            trend_result = mann_kendall_test(filtered_times, filtered_abundances)
        elif trend_method == 'linear_regression':
            trend_result = linear_regression_test(filtered_times, filtered_abundances)
        else:
            continue
        
        trend_result['feature_id'] = feature
        trend_result['n_timepoints'] = len(filtered_abundances)
        trend_result['mean_abundance'] = np.mean(filtered_abundances)
        trend_result['abundance_range'] = np.max(filtered_abundances) - np.min(filtered_abundances)
        
        trend_results.append(trend_result)
    
    trend_df = pd.DataFrame(trend_results)
    
    if not trend_df.empty:
        # Apply multiple testing correction
        from statsmodels.stats.multitest import multipletests
        
        p_values = trend_df['p_value'].values
        rejected, p_adjusted, _, _ = multipletests(p_values, method='fdr_bh')
        
        trend_df['p_adjusted'] = p_adjusted
        trend_df['significant'] = p_adjusted < significance_threshold
        
        # Sort by significance
        trend_df = trend_df.sort_values('p_adjusted')
        
        logging.info(f"Detected {np.sum(trend_df['significant'])} significant temporal trends")
    
    return trend_df

def mann_kendall_test(time_values, data_values):
    """Perform Mann-Kendall test for trend detection"""
    n = len(data_values)
    
    # Calculate Mann-Kendall statistic
    s = 0
    for i in range(n-1):
        for j in range(i+1, n):
            if data_values[j] > data_values[i]:
                s += 1
            elif data_values[j] < data_values[i]:
                s -= 1
    
    # Calculate variance
    var_s = n * (n - 1) * (2 * n + 5) / 18
    
    # Calculate Z statistic
    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0
    
    # Calculate p-value (two-tailed)
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))
    
    # Determine trend direction
    if s > 0:
        trend = 'increasing'
    elif s < 0:
        trend = 'decreasing'
    else:
        trend = 'no_trend'
    
    return {
        'statistic': s,
        'z_score': z,
        'p_value': p_value,
        'trend_direction': trend,
        'test_method': 'mann_kendall'
    }

def linear_regression_test(time_values, data_values):
    """Perform linear regression test for trend detection"""
    from scipy.stats import linregress
    
    slope, intercept, r_value, p_value, std_err = linregress(time_values, data_values)
    
    # Determine trend direction
    if slope > 0:
        trend = 'increasing'
    elif slope < 0:
        trend = 'decreasing'
    else:
        trend = 'no_trend'
    
    return {
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_value**2,
        'p_value': p_value,
        'std_error': std_err,
        'trend_direction': trend,
        'test_method': 'linear_regression'
    }

def detect_changepoints(abundance_data, metadata_df, config):
    """Detect changepoints in time series data"""
    logging.info("Detecting changepoints in temporal data")
    
    method = config.get('method', 'pelt')
    penalty = config.get('penalty', 'BIC')
    min_size = config.get('min_size', 3)
    
    changepoint_results = []
    
    try:
        import ruptures as rpt
        
        for feature in abundance_data.index:
            abundances = abundance_data.loc[feature].values
            
            # Skip features with too many zeros
            if np.sum(abundances > 0) < min_size * 2:
                continue
            
            # Apply smoothing to reduce noise
            if len(abundances) > 5:
                smoothed = savgol_filter(abundances, window_length=min(5, len(abundances)//2*2+1), polyorder=1)
            else:
                smoothed = abundances
            
            # Detect changepoints
            if method == 'pelt':
                algo = rpt.Pelt(model="rbf").fit(smoothed.reshape(-1, 1))
            elif method == 'binseg':
                algo = rpt.Binseg(model="rbf").fit(smoothed.reshape(-1, 1))
            else:
                continue
            
            # Get changepoints
            changepoints = algo.predict(pen=penalty)
            
            # Remove the last point (always equal to the length)
            changepoints = changepoints[:-1]
            
            if len(changepoints) > 0:
                # Convert to actual time points
                changepoint_times = metadata_df['time_interval'].iloc[changepoints].tolist()
                changepoint_dates = metadata_df['collection_date'].iloc[changepoints].tolist()
                
                changepoint_results.append({
                    'feature_id': feature,
                    'n_changepoints': len(changepoints),
                    'changepoint_indices': changepoints,
                    'changepoint_times': changepoint_times,
                    'changepoint_dates': [str(date) for date in changepoint_dates],
                    'detection_method': method
                })
    
    except ImportError:
        logging.warning("ruptures package not available for changepoint detection")
        # Simple changepoint detection based on variance
        changepoint_results = simple_changepoint_detection(abundance_data, metadata_df, min_size)
    
    changepoint_df = pd.DataFrame(changepoint_results)
    
    if not changepoint_df.empty:
        logging.info(f"Detected changepoints in {len(changepoint_df)} features")
    
    return changepoint_df

def simple_changepoint_detection(abundance_data, metadata_df, min_size):
    """Simple changepoint detection based on variance changes"""
    logging.info("Using simple variance-based changepoint detection")
    
    changepoint_results = []
    
    for feature in abundance_data.index:
        abundances = abundance_data.loc[feature].values
        
        if len(abundances) < min_size * 2:
            continue
        
        # Calculate rolling variance
        window_size = max(3, len(abundances) // 4)
        rolling_var = pd.Series(abundances).rolling(window=window_size).var()
        
        # Find points where variance changes significantly
        var_changes = np.abs(np.diff(rolling_var.fillna(0)))
        
        # Find peaks in variance changes
        threshold = np.percentile(var_changes, 90)
        changepoints = np.where(var_changes > threshold)[0]
        
        if len(changepoints) > 0:
            # Filter changepoints that are too close
            filtered_changepoints = []
            for cp in changepoints:
                if not filtered_changepoints or cp - filtered_changepoints[-1] >= min_size:
                    filtered_changepoints.append(cp)
            
            if filtered_changepoints:
                changepoint_times = metadata_df['time_interval'].iloc[filtered_changepoints].tolist()
                changepoint_dates = metadata_df['collection_date'].iloc[filtered_changepoints].tolist()
                
                changepoint_results.append({
                    'feature_id': feature,
                    'n_changepoints': len(filtered_changepoints),
                    'changepoint_indices': filtered_changepoints,
                    'changepoint_times': changepoint_times,
                    'changepoint_dates': [str(date) for date in changepoint_dates],
                    'detection_method': 'variance_based'
                })
    
    return changepoint_results

def perform_seasonal_decomposition(abundance_data, metadata_df, config):
    """Perform seasonal decomposition of time series"""
    logging.info("Performing seasonal decomposition")
    
    if not config.get('enabled', True):
        return pd.DataFrame()
    
    model = config.get('model', 'additive')
    period = config.get('period', None)
    
    seasonal_results = []
    
    try:
        from statsmodels.tsa.seasonal import seasonal_decompose
        
        # Determine period if not specified
        if period is None:
            # Estimate period based on data frequency
            time_diffs = np.diff(metadata_df['time_interval'])
            median_interval = np.median(time_diffs)
            
            if median_interval <= 1:  # Daily data
                period = 7  # Weekly seasonality
            elif median_interval <= 7:  # Weekly data
                period = 4  # Monthly seasonality
            else:
                period = 12  # Yearly seasonality for monthly data
        
        for feature in abundance_data.index:
            abundances = abundance_data.loc[feature].values
            
            # Skip features with insufficient data or too many zeros
            if len(abundances) < period * 2 or np.sum(abundances > 0) < len(abundances) * 0.3:
                continue
            
            try:
                # Perform seasonal decomposition
                decomposition = seasonal_decompose(
                    abundances, 
                    model=model, 
                    period=period,
                    extrapolate_trend='freq'
                )
                
                # Calculate seasonal strength
                seasonal_var = np.var(decomposition.seasonal)
                residual_var = np.var(decomposition.resid[~np.isnan(decomposition.resid)])
                seasonal_strength = seasonal_var / (seasonal_var + residual_var) if (seasonal_var + residual_var) > 0 else 0
                
                # Calculate trend strength
                trend_var = np.var(decomposition.trend[~np.isnan(decomposition.trend)])
                trend_strength = trend_var / (trend_var + residual_var) if (trend_var + residual_var) > 0 else 0
                
                seasonal_results.append({
                    'feature_id': feature,
                    'seasonal_strength': seasonal_strength,
                    'trend_strength': trend_strength,
                    'period_used': period,
                    'model_type': model,
                    'mean_abundance': np.mean(abundances[abundances > 0]) if np.sum(abundances > 0) > 0 else 0
                })
                
            except Exception as e:
                logging.debug(f"Could not decompose feature {feature}: {str(e)}")
                continue
    
    except ImportError:
        logging.warning("statsmodels not available for seasonal decomposition")
    
    seasonal_df = pd.DataFrame(seasonal_results)
    
    if not seasonal_df.empty:
        logging.info(f"Performed seasonal decomposition for {len(seasonal_df)} features")
    
    return seasonal_df

def analyze_community_stability(abundance_data, metadata_df):
    """Analyze community stability over time"""
    logging.info("Analyzing community stability")
    
    stability_metrics = []
    
    # Calculate diversity metrics over time
    for i, sample in enumerate(abundance_data.columns):
        abundances = abundance_data[sample].values
        non_zero_abundances = abundances[abundances > 0]
        
        # Basic diversity metrics
        observed_otus = len(non_zero_abundances)
        shannon = -np.sum(non_zero_abundances * np.log(non_zero_abundances + 1e-10))
        simpson = 1 - np.sum(non_zero_abundances ** 2)
        
        # Evenness
        pielou = shannon / np.log(observed_otus) if observed_otus > 1 else 0
        
        # Community dissimilarity from first sample
        if i == 0:
            reference_community = abundances
            bray_curtis_distance = 0
        else:
            bray_curtis_distance = calculate_bray_curtis_distance(reference_community, abundances)
        
        # Community dissimilarity from previous sample
        if i == 0:
            temporal_distance = 0
        else:
            prev_abundances = abundance_data.iloc[:, i-1].values
            temporal_distance = calculate_bray_curtis_distance(prev_abundances, abundances)
        
        stability_metrics.append({
            'sample_id': sample,
            'time_point': i,
            'time_interval': metadata_df['time_interval'].iloc[i],
            'observed_otus': observed_otus,
            'shannon_diversity': shannon,
            'simpson_diversity': simpson,
            'pielou_evenness': pielou,
            'bray_curtis_from_baseline': bray_curtis_distance,
            'bray_curtis_from_previous': temporal_distance
        })
    
    stability_df = pd.DataFrame(stability_metrics)
    
    # Calculate overall stability metrics
    overall_stability = {
        'mean_shannon': stability_df['shannon_diversity'].mean(),
        'cv_shannon': stability_df['shannon_diversity'].std() / stability_df['shannon_diversity'].mean(),
        'mean_bray_curtis_baseline': stability_df['bray_curtis_from_baseline'].mean(),
        'mean_bray_curtis_temporal': stability_df['bray_curtis_from_previous'].mean(),
        'community_stability_index': 1 - stability_df['bray_curtis_from_previous'].mean()
    }
    
    return stability_df, overall_stability

def calculate_bray_curtis_distance(community1, community2):
    """Calculate Bray-Curtis dissimilarity between two communities"""
    # Normalize to relative abundances
    rel_abund1 = community1 / np.sum(community1) if np.sum(community1) > 0 else community1
    rel_abund2 = community2 / np.sum(community2) if np.sum(community2) > 0 else community2
    
    # Calculate Bray-Curtis dissimilarity
    numerator = np.sum(np.abs(rel_abund1 - rel_abund2))
    denominator = np.sum(rel_abund1 + rel_abund2)
    
    if denominator > 0:
        return numerator / denominator
    else:
        return 0

def create_temporal_visualizations(abundance_data, metadata_df, trend_results, 
                                 changepoint_results, stability_data, output_file):
    """Create comprehensive temporal visualizations"""
    logging.info("Creating temporal visualizations")
    
    fig = plt.figure(figsize=(20, 16))
    
    # 1. Community stability over time
    ax1 = plt.subplot(3, 3, 1)
    plt.plot(stability_data['time_interval'], stability_data['shannon_diversity'], 'o-', alpha=0.7)
    plt.xlabel('Time')
    plt.ylabel('Shannon Diversity')
    plt.title('Shannon Diversity Over Time')
    plt.grid(True, alpha=0.3)
    
    # 2. Community dissimilarity from baseline
    ax2 = plt.subplot(3, 3, 2)
    plt.plot(stability_data['time_interval'], stability_data['bray_curtis_from_baseline'], 'o-', alpha=0.7, color='red')
    plt.xlabel('Time')
    plt.ylabel('Bray-Curtis Distance')
    plt.title('Community Drift from Baseline')
    plt.grid(True, alpha=0.3)
    
    # 3. Trend direction distribution
    ax3 = plt.subplot(3, 3, 3)
    if not trend_results.empty:
        trend_counts = trend_results['trend_direction'].value_counts()
        plt.bar(trend_counts.index, trend_counts.values)
        plt.xlabel('Trend Direction')
        plt.ylabel('Number of Features')
        plt.title('Distribution of Temporal Trends')
        plt.xticks(rotation=45)
    
    # 4. Top trending features
    ax4 = plt.subplot(3, 3, 4)
    if not trend_results.empty:
        significant_trends = trend_results[trend_results['significant']].head(5)
        if len(significant_trends) > 0:
            feature_to_plot = significant_trends.iloc[0]['feature_id']
            abundances = abundance_data.loc[feature_to_plot].values
            plt.plot(metadata_df['time_interval'], abundances, 'o-', alpha=0.7)
            plt.xlabel('Time')
            plt.ylabel('Abundance')
            plt.title(f'Top Trending Feature: {feature_to_plot[:20]}...')
            plt.grid(True, alpha=0.3)
    
    # 5. Changepoint frequency over time
    ax5 = plt.subplot(3, 3, 5)
    if not changepoint_results.empty:
        all_changepoints = []
        for _, row in changepoint_results.iterrows():
            all_changepoints.extend(row['changepoint_times'])
        
        if all_changepoints:
            plt.hist(all_changepoints, bins=20, alpha=0.7, edgecolor='black')
            plt.xlabel('Time')
            plt.ylabel('Number of Changepoints')
            plt.title('Changepoint Frequency Distribution')
    
    # 6. Feature abundance heatmap (top variable features)
    ax6 = plt.subplot(3, 3, 6)
    # Select top 20 most variable features
    feature_variance = abundance_data.var(axis=1).sort_values(ascending=False)
    top_features = feature_variance.head(20).index
    
    heatmap_data = abundance_data.loc[top_features]
    # Normalize each feature
    heatmap_normalized = heatmap_data.div(heatmap_data.max(axis=1), axis=0)
    
    sns.heatmap(heatmap_normalized, cmap='viridis', cbar=True, 
               xticklabels=False, yticklabels=[f[:15] for f in top_features])
    plt.title('Top Variable Features Over Time')
    plt.xlabel('Time Points')
    
    # 7. Temporal autocorrelation
    ax7 = plt.subplot(3, 3, 7)
    # Calculate autocorrelation for Shannon diversity
    shannon_values = stability_data['shannon_diversity'].values
    if len(shannon_values) > 3:
        autocorr = [np.corrcoef(shannon_values[:-lag], shannon_values[lag:])[0,1] 
                   for lag in range(1, min(10, len(shannon_values)//2))]
        plt.plot(range(1, len(autocorr)+1), autocorr, 'o-')
        plt.xlabel('Lag')
        plt.ylabel('Autocorrelation')
        plt.title('Shannon Diversity Autocorrelation')
        plt.grid(True, alpha=0.3)
    
    # 8. Community stability metrics
    ax8 = plt.subplot(3, 3, 8)
    stability_metrics = ['shannon_diversity', 'simpson_diversity', 'pielou_evenness']
    for metric in stability_metrics:
        normalized_values = (stability_data[metric] - stability_data[metric].min()) / (stability_data[metric].max() - stability_data[metric].min())
        plt.plot(stability_data['time_interval'], normalized_values, 'o-', alpha=0.7, label=metric)
    plt.xlabel('Time')
    plt.ylabel('Normalized Value')
    plt.title('Normalized Diversity Metrics')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 9. Feature emergence/disappearance
    ax9 = plt.subplot(3, 3, 9)
    # Count features present at each time point
    presence_counts = (abundance_data > 0).sum(axis=0)
    plt.plot(metadata_df['time_interval'], presence_counts, 'o-', alpha=0.7)
    plt.xlabel('Time')
    plt.ylabel('Number of Present Features')
    plt.title('Feature Richness Over Time')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Temporal analysis of microbiome data')
    parser.add_argument('--abundance-matrix', required=True,
                        help='Abundance matrix file')
    parser.add_argument('--novelty-data', required=True,
                        help='Pan-novelty analysis results file')
    parser.add_argument('--metadata', required=True,
                        help='Sample metadata file')
    parser.add_argument('--output-temporal-trends', required=True,
                        help='Output temporal trends file')
    parser.add_argument('--output-changepoint-detection', required=True,
                        help='Output changepoint detection file')
    parser.add_argument('--output-seasonal-decomposition', required=True,
                        help='Output seasonal decomposition file')
    parser.add_argument('--output-trend-significance', required=True,
                        help='Output trend significance tests file')
    parser.add_argument('--output-temporal-plots', required=True,
                        help='Output temporal visualization HTML file')
    parser.add_argument('--log-file', required=True,
                        help='Log file path')
    
    # Analysis parameters
    parser.add_argument('--time-unit', choices=['days', 'weeks', 'months'], default='days',
                        help='Time unit for analysis')
    parser.add_argument('--trend-method', choices=['mann_kendall', 'linear_regression'], default='mann_kendall',
                        help='Method for trend detection')
    parser.add_argument('--changepoint-method', choices=['pelt', 'binseg', 'variance'], default='pelt',
                        help='Method for changepoint detection')
    parser.add_argument('--seasonal-model', choices=['additive', 'multiplicative'], default='additive',
                        help='Seasonal decomposition model')
    parser.add_argument('--enable-seasonal', action='store_true', default=True,
                        help='Enable seasonal decomposition')
    
    args = parser.parse_args()
    
    # Setup logging
    os.makedirs(os.path.dirname(args.log_file), exist_ok=True)
    setup_logging(args.log_file)
    
    try:
        logging.info("Starting temporal analysis")
        
        # Check if temporal analysis is enabled
        if not args.enable_seasonal and args.changepoint_method and args.trend_method:
            logging.info("Temporal analysis is enabled")
        
        # Load and prepare data
        abundance_data, novelty_data, metadata_df = load_and_prepare_temporal_data(
            args.abundance_matrix, args.novelty_data, args.metadata, {}
        )
        
        # Calculate time intervals
        metadata_with_time = calculate_time_intervals(metadata_df, args.time_unit)
        
        # Trend detection configuration
        trend_config = {
            'method': args.trend_method,
            'significance_threshold': 0.05
        }
        
        # Detect temporal trends
        trend_results = detect_temporal_trends(abundance_data, metadata_with_time, trend_config)
        
        # Changepoint detection configuration
        changepoint_config = {
            'method': args.changepoint_method,
            'penalty': 'BIC',
            'min_size': 3
        }
        
        # Detect changepoints
        changepoint_results = detect_changepoints(abundance_data, metadata_with_time, changepoint_config)
        
        # Seasonal decomposition configuration
        seasonal_config = {
            'enabled': args.enable_seasonal,
            'model': args.seasonal_model,
            'period': None  # Auto-detect
        }
        
        # Perform seasonal decomposition
        seasonal_results = perform_seasonal_decomposition(abundance_data, metadata_with_time, seasonal_config)
        
        # Analyze community stability
        stability_data, overall_stability = analyze_community_stability(abundance_data, metadata_with_time)
        
        # Create output directories
        for output_file in [args.output_temporal_trends, args.output_changepoint_detection,
                           args.output_seasonal_decomposition, args.output_trend_significance,
                           args.output_temporal_plots]:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Save results
        logging.info("Saving temporal analysis results")
        
        # Temporal trends
        if not trend_results.empty:
            trend_results.to_csv(args.output_temporal_trends, sep='\t', index=False)
        else:
            pd.DataFrame().to_csv(args.output_temporal_trends, sep='\t', index=False)
        
        # Changepoint detection
        if not changepoint_results.empty:
            changepoint_results.to_csv(args.output_changepoint_detection, sep='\t', index=False)
        else:
            pd.DataFrame().to_csv(args.output_changepoint_detection, sep='\t', index=False)
        
        # Seasonal decomposition
        if not seasonal_results.empty:
            seasonal_results.to_csv(args.output_seasonal_decomposition, sep='\t', index=False)
        else:
            pd.DataFrame().to_csv(args.output_seasonal_decomposition, sep='\t', index=False)
        
        # Trend significance tests (same as temporal trends for now)
        if not trend_results.empty:
            significant_trends = trend_results[trend_results['significant']]
            significant_trends.to_csv(args.output_trend_significance, sep='\t', index=False)
        else:
            pd.DataFrame().to_csv(args.output_trend_significance, sep='\t', index=False)
        
        # Create visualizations
        plot_file = args.output_temporal_plots.replace('.html', '.png')
        create_temporal_visualizations(
            abundance_data, metadata_with_time, trend_results,
            changepoint_results, stability_data, plot_file
        )
        
        # Create HTML report with embedded plot
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Temporal Analysis Results</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                h1, h2 {{ color: #333; }}
                .summary {{ background-color: #f5f5f5; padding: 20px; border-radius: 5px; }}
                .metric {{ margin: 10px 0; }}
            </style>
        </head>
        <body>
            <h1>Temporal Analysis Results</h1>
            
            <div class="summary">
                <h2>Analysis Summary</h2>
                <div class="metric"><strong>Time range:</strong> {metadata_df['collection_date'].min()} to {metadata_df['collection_date'].max()}</div>
                <div class="metric"><strong>Number of time points:</strong> {len(metadata_with_time)}</div>
                <div class="metric"><strong>Significant trends detected:</strong> {np.sum(trend_results['significant']) if not trend_results.empty else 0}</div>
                <div class="metric"><strong>Features with changepoints:</strong> {len(changepoint_results) if not changepoint_results.empty else 0}</div>
                <div class="metric"><strong>Community stability index:</strong> {overall_stability['community_stability_index']:.3f}</div>
            </div>
            
            <h2>Temporal Visualization</h2>
            <img src="{os.path.basename(plot_file)}" alt="Temporal Analysis Plots" style="max-width: 100%;">
            
        </body>
        </html>
        """
        
        with open(args.output_temporal_plots, 'w') as f:
            f.write(html_content)
        
        logging.info("Temporal analysis completed successfully")
        
        # Print summary
        print("\n" + "="*80)
        print("TEMPORAL ANALYSIS SUMMARY")
        print("="*80)
        print(f"Time range: {metadata_df['collection_date'].min()} to {metadata_df['collection_date'].max()}")
        print(f"Number of time points: {len(metadata_with_time)}")
        print(f"Time unit: {args.time_unit}")
        
        if not trend_results.empty:
            print(f"Significant trends detected: {np.sum(trend_results['significant'])}")
            print(f"Increasing trends: {np.sum((trend_results['trend_direction'] == 'increasing') & trend_results['significant'])}")
            print(f"Decreasing trends: {np.sum((trend_results['trend_direction'] == 'decreasing') & trend_results['significant'])}")
        
        if not changepoint_results.empty:
            print(f"Features with changepoints: {len(changepoint_results)}")
            avg_changepoints = changepoint_results['n_changepoints'].mean()
            print(f"Average changepoints per feature: {avg_changepoints:.2f}")
        
        print(f"Community stability index: {overall_stability['community_stability_index']:.3f}")
        print(f"Mean Shannon diversity: {overall_stability['mean_shannon']:.3f}")
        print(f"Shannon diversity CV: {overall_stability['cv_shannon']:.3f}")
        
        print("="*80)
        
    except Exception as e:
        logging.error(f"Error in temporal analysis: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()

