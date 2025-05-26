#!/usr/bin/env python3
"""
Comprehensive Report Generation
Generate comprehensive HTML and PDF reports from cross-sample analysis results
"""

import sys
import os
import pandas as pd
import numpy as np
import json
import logging
import argparse
from pathlib import Path
from datetime import datetime
import base64

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
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

def load_analysis_results(file_paths):
    """Load all analysis results"""
    logging.info("Loading analysis results")
    
    results = {}
    
    # Define expected files and their types
    file_types = {
        'abundance_matrix': 'csv',
        'mag_abundance_matrix': 'csv',
        'dereplicated_mags': 'csv',
        'vae_embeddings': 'npy',
        'vae_metrics': 'json',
        'isolation_forest_results': 'csv',
        'pan_novelty_results': 'csv',
        'novelty_emergence': 'csv',
        'temporal_trends': 'csv',
        'pathogen_risk': 'csv',
        'high_risk_alerts': 'csv',
        'metadata': 'csv'
    }
    
    for key, file_path in file_paths.items():
        if key in file_types and file_path and os.path.exists(file_path):
            try:
                file_type = file_types[key]
                
                if file_type == 'csv':
                    if key == 'metadata':
                        results[key] = pd.read_csv(file_path, sep='\t')
                    else:
                        results[key] = pd.read_csv(file_path, sep='\t', index_col=0)
                elif file_type == 'npy':
                    results[key] = np.load(file_path)
                elif file_type == 'json':
                    with open(file_path, 'r') as f:
                        results[key] = json.load(f)
                
                logging.info(f"Loaded {key}: {file_path}")
                
            except Exception as e:
                logging.warning(f"Could not load {key} from {file_path}: {str(e)}")
                results[key] = None
        else:
            logging.warning(f"File not found for {key}: {file_path}")
            results[key] = None
    
    return results

def generate_executive_summary(results):
    """Generate executive summary statistics"""
    logging.info("Generating executive summary")
    
    summary = {
        'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'dataset_overview': {},
        'key_findings': {},
        'risk_assessment': {},
        'recommendations': []
    }
    
    # Dataset overview
    if results['abundance_matrix'] is not None:
        abundance_df = results['abundance_matrix']
        summary['dataset_overview'] = {
            'total_samples': abundance_df.shape[1],
            'total_features': abundance_df.shape[0],
            'data_sparsity': (abundance_df == 0).sum().sum() / abundance_df.size,
            'mean_features_per_sample': (abundance_df > 0).sum(axis=0).mean()
        }
    
    # MAG analysis
    if results['dereplicated_mags'] is not None:
        mag_df = results['dereplicated_mags']
        summary['dataset_overview']['total_mags'] = len(mag_df)
        summary['dataset_overview']['high_quality_mags'] = len(mag_df[
            (mag_df['completeness'] >= 90) & (mag_df['contamination'] <= 5)
        ])
    
    # ML analysis
    if results['vae_metrics'] is not None:
        vae_metrics = results['vae_metrics']
        summary['key_findings']['ml_analysis'] = {
            'vae_reconstruction_correlation': vae_metrics.get('mean_reconstruction_correlation', 0),
            'community_embedding_quality': 'Good' if vae_metrics.get('mean_reconstruction_correlation', 0) > 0.8 else 'Moderate'
        }
    
    # Novelty analysis
    if results['pan_novelty_results'] is not None:
        novelty_df = results['pan_novelty_results']
        high_novelty = novelty_df[novelty_df['mean'] > 0.8] if 'mean' in novelty_df.columns else pd.DataFrame()
        summary['key_findings']['novelty_analysis'] = {
            'novel_features_detected': len(novelty_df),
            'high_novelty_features': len(high_novelty),
            'novelty_discovery_rate': len(high_novelty) / len(novelty_df) if len(novelty_df) > 0 else 0
        }
    
    # Temporal analysis
    if results['temporal_trends'] is not None:
        temporal_df = results['temporal_trends']
        if 'significant' in temporal_df.columns:
            significant_trends = temporal_df[temporal_df['significant']] if len(temporal_df) > 0 else pd.DataFrame()
            summary['key_findings']['temporal_analysis'] = {
                'significant_temporal_trends': len(significant_trends),
                'increasing_trends': len(significant_trends[significant_trends['trend_direction'] == 'increasing']) if len(significant_trends) > 0 else 0,
                'decreasing_trends': len(significant_trends[significant_trends['trend_direction'] == 'decreasing']) if len(significant_trends) > 0 else 0
            }
    
    # Risk assessment
    if results['pathogen_risk'] is not None:
        pathogen_df = results['pathogen_risk']
        high_risk_count = len(pathogen_df[pathogen_df['risk_level'] == 'high']) if 'risk_level' in pathogen_df.columns else 0
        summary['risk_assessment'] = {
            'total_pathogen_detections': len(pathogen_df),
            'high_risk_detections': high_risk_count,
            'risk_level': 'HIGH' if high_risk_count > 0 else 'MODERATE' if len(pathogen_df) > 0 else 'LOW'
        }
    
    # Generate recommendations
    summary['recommendations'] = generate_recommendations(results, summary)
    
    return summary

def generate_recommendations(results, summary):
    """Generate actionable recommendations based on analysis results"""
    recommendations = []
    
    # Risk-based recommendations
    risk_level = summary.get('risk_assessment', {}).get('risk_level', 'LOW')
    if risk_level == 'HIGH':
        recommendations.append({
            'priority': 'CRITICAL',
            'category': 'Pathogen Detection',
            'recommendation': 'Immediate investigation of high-risk pathogen detections required',
            'action_items': [
                'Review high-risk alerts in detail',
                'Implement enhanced biosafety protocols',
                'Consider additional confirmatory testing',
                'Monitor affected samples closely'
            ]
        })
    
    # Novelty-based recommendations
    novelty_rate = summary.get('key_findings', {}).get('novelty_analysis', {}).get('novelty_discovery_rate', 0)
    if novelty_rate > 0.1:
        recommendations.append({
            'priority': 'HIGH',
            'category': 'Novel Discovery',
            'recommendation': 'Significant novel genetic material detected requiring further investigation',
            'action_items': [
                'Sequence novel contigs for detailed characterization',
                'Perform functional annotation of novel features',
                'Investigate potential biotechnological applications',
                'Consider intellectual property implications'
            ]
        })
    
    # Temporal analysis recommendations
    temporal_findings = summary.get('key_findings', {}).get('temporal_analysis', {})
    if temporal_findings:
        significant_trends = temporal_findings.get('significant_temporal_trends', 0)
        if significant_trends > 10:
            recommendations.append({
                'priority': 'MEDIUM',
                'category': 'Temporal Dynamics',
                'recommendation': 'Significant temporal changes detected in microbial community',
                'action_items': [
                    'Investigate environmental factors driving changes',
                    'Implement more frequent monitoring',
                    'Assess stability of ecosystem services',
                    'Consider intervention strategies if needed'
                ]
            })
    
    # Data quality recommendations
    sparsity = summary.get('dataset_overview', {}).get('data_sparsity', 0)
    if sparsity > 0.8:
        recommendations.append({
            'priority': 'MEDIUM',
            'category': 'Data Quality',
            'recommendation': 'High data sparsity detected - consider optimizing sampling strategy',
            'action_items': [
                'Evaluate sequencing depth adequacy',
                'Consider targeted enrichment protocols',
                'Optimize DNA extraction methods',
                'Increase sampling frequency if resources allow'
            ]
        })
    
    # General recommendations
    recommendations.append({
        'priority': 'LOW',
        'category': 'Continuous Monitoring',
        'recommendation': 'Establish regular monitoring protocol for ongoing surveillance',
        'action_items': [
            'Schedule regular re-analysis of samples',
            'Implement automated alert systems',
            'Maintain updated pathogen databases',
            'Train personnel on interpretation of results'
        ]
    })
    
    return recommendations

def create_summary_visualizations(results, output_dir):
    """Create summary visualizations for the report"""
    logging.info("Creating summary visualizations")
    
    visualization_files = {}
    
    # 1. Dataset overview
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Sample composition
    if results['abundance_matrix'] is not None:
        abundance_df = results['abundance_matrix']
        
        # Features per sample
        features_per_sample = (abundance_df > 0).sum(axis=0)
        axes[0, 0].hist(features_per_sample, bins=20, alpha=0.7, edgecolor='black')
        axes[0, 0].set_title('Features per Sample Distribution')
        axes[0, 0].set_xlabel('Number of Features')
        axes[0, 0].set_ylabel('Number of Samples')
        
        # Abundance distribution
        all_abundances = abundance_df.values.flatten()
        non_zero_abundances = all_abundances[all_abundances > 0]
        axes[0, 1].hist(np.log10(non_zero_abundances + 1), bins=30, alpha=0.7, edgecolor='black')
        axes[0, 1].set_title('Abundance Distribution (log10)')
        axes[0, 1].set_xlabel('log10(Abundance + 1)')
        axes[0, 1].set_ylabel('Frequency')
    
    # MAG quality if available
    if results['dereplicated_mags'] is not None:
        mag_df = results['dereplicated_mags']
        
        axes[1, 0].scatter(mag_df['completeness'], mag_df['contamination'], alpha=0.6)
        axes[1, 0].set_xlabel('Completeness (%)')
        axes[1, 0].set_ylabel('Contamination (%)')
        axes[1, 0].set_title('MAG Quality Assessment')
        axes[1, 0].axhline(y=5, color='red', linestyle='--', label='5% contamination')
        axes[1, 0].axvline(x=90, color='red', linestyle='--', label='90% completeness')
        axes[1, 0].legend()
    
    # Risk assessment
    if results['pathogen_risk'] is not None:
        pathogen_df = results['pathogen_risk']
        if 'risk_level' in pathogen_df.columns:
            risk_counts = pathogen_df['risk_level'].value_counts()
            colors = {'high': 'red', 'medium': 'orange', 'low': 'yellow'}
            
            bars = axes[1, 1].bar(risk_counts.index, risk_counts.values)
            for bar, risk_level in zip(bars, risk_counts.index):
                bar.set_color(colors.get(risk_level, 'gray'))
            
            axes[1, 1].set_title('Pathogen Risk Assessment')
            axes[1, 1].set_xlabel('Risk Level')
            axes[1, 1].set_ylabel('Number of Detections')
    
    plt.tight_layout()
    overview_file = os.path.join(output_dir, 'dataset_overview.png')
    plt.savefig(overview_file, dpi=300, bbox_inches='tight')
    plt.close()
    visualization_files['dataset_overview'] = overview_file
    
    # 2. Analysis results summary
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Novelty analysis
    if results['pan_novelty_results'] is not None:
        novelty_df = results['pan_novelty_results']
        if 'mean' in novelty_df.columns:
            axes[0, 0].hist(novelty_df['mean'], bins=20, alpha=0.7, edgecolor='black')
            axes[0, 0].set_title('Novelty Score Distribution')
            axes[0, 0].set_xlabel('Mean Novelty Score')
            axes[0, 0].set_ylabel('Number of Features')
            axes[0, 0].axvline(0.8, color='red', linestyle='--', label='High novelty')
            axes[0, 0].legend()
    
    # Temporal trends
    if results['temporal_trends'] is not None:
        temporal_df = results['temporal_trends']
        if 'trend_direction' in temporal_df.columns:
            trend_counts = temporal_df['trend_direction'].value_counts()
            axes[0, 1].pie(trend_counts.values, labels=trend_counts.index, autopct='%1.1f%%')
            axes[0, 1].set_title('Temporal Trend Directions')
    
    # ML analysis results
    if results['isolation_forest_results'] is not None:
        if_df = results['isolation_forest_results']
        if 'anomaly_probability' in if_df.columns:
            axes[1, 0].hist(if_df['anomaly_probability'], bins=20, alpha=0.7, edgecolor='black')
            axes[1, 0].set_title('Anomaly Probability Distribution')
            axes[1, 0].set_xlabel('Anomaly Probability')
            axes[1, 0].set_ylabel('Number of Samples')
    
    # Empty plot for future use
    axes[1, 1].text(0.5, 0.5, 'Additional Analysis\nResults', ha='center', va='center', 
                   transform=axes[1, 1].transAxes, fontsize=12)
    axes[1, 1].set_title('Future Extensions')
    
    plt.tight_layout()
    analysis_file = os.path.join(output_dir, 'analysis_summary.png')
    plt.savefig(analysis_file, dpi=300, bbox_inches='tight')
    plt.close()
    visualization_files['analysis_summary'] = analysis_file
    
    return visualization_files

def encode_image_to_base64(image_path):
    """Encode image to base64 for embedding in HTML"""
    if os.path.exists(image_path):
        with open(image_path, 'rb') as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    return None

def generate_html_report(results, summary, visualization_files, config, output_file):
    """Generate comprehensive HTML report"""
    logging.info("Generating HTML report")
    
    # Encode images to base64 for embedding
    embedded_images = {}
    for key, img_path in visualization_files.items():
        encoded_img = encode_image_to_base64(img_path)
        if encoded_img:
            embedded_images[key] = encoded_img
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Comprehensive Metagenomics Analysis Report</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                line-height: 1.6;
                margin: 0;
                padding: 20px;
                background-color: #f5f5f5;
            }}
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background: white;
                padding: 30px;
                border-radius: 10px;
                box-shadow: 0 0 20px rgba(0,0,0,0.1);
            }}
            .header {{
                text-align: center;
                border-bottom: 3px solid #2c3e50;
                padding-bottom: 20px;
                margin-bottom: 30px;
            }}
            .header h1 {{
                color: #2c3e50;
                margin-bottom: 10px;
                font-size: 2.5em;
            }}
            .header .subtitle {{
                color: #7f8c8d;
                font-size: 1.2em;
            }}
            .executive-summary {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 25px;
                border-radius: 10px;
                margin-bottom: 30px;
            }}
            .section {{
                margin-bottom: 40px;
                padding: 20px;
                border-left: 4px solid #3498db;
                background: #f8f9fa;
                border-radius: 5px;
            }}
            .section h2 {{
                color: #2c3e50;
                margin-top: 0;
                font-size: 1.8em;
            }}
            .section h3 {{
                color: #34495e;
                margin-top: 25px;
                font-size: 1.4em;
            }}
            .metrics-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                gap: 20px;
                margin: 20px 0;
            }}
            .metric-card {{
                background: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                text-align: center;
            }}
            .metric-value {{
                font-size: 2em;
                font-weight: bold;
                color: #3498db;
            }}
            .metric-label {{
                color: #7f8c8d;
                margin-top: 5px;
            }}
            .risk-high {{ color: #e74c3c; }}
            .risk-medium {{ color: #f39c12; }}
            .risk-low {{ color: #27ae60; }}
            .alert-box {{
                padding: 15px;
                border-radius: 5px;
                margin: 15px 0;
            }}
            .alert-critical {{
                background: #ffebee;
                border-left: 4px solid #e74c3c;
                color: #c62828;
            }}
            .alert-high {{
                background: #fff3e0;
                border-left: 4px solid #f39c12;
                color: #e65100;
            }}
            .alert-medium {{
                background: #e8f5e8;
                border-left: 4px solid #4caf50;
                color: #2e7d32;
            }}
            .recommendations {{
                background: #f0f8ff;
                padding: 20px;
                border-radius: 10px;
                margin: 20px 0;
            }}
            .recommendation-item {{
                margin: 15px 0;
                padding: 15px;
                background: white;
                border-left: 4px solid #3498db;
                border-radius: 5px;
            }}
            .action-items {{
                margin-top: 10px;
            }}
            .action-items li {{
                margin: 5px 0;
            }}
            .visualization {{
                text-align: center;
                margin: 30px 0;
            }}
            .visualization img {{
                max-width: 100%;
                height: auto;
                border: 1px solid #ddd;
                border-radius: 8px;
                box-shadow: 0 4px 8px rgba(0,0,0,0.1);
            }}
            .data-table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
            }}
            .data-table th,
            .data-table td {{
                padding: 12px;
                text-align: left;
                border-bottom: 1px solid #ddd;
            }}
            .data-table th {{
                background-color: #f2f2f2;
                font-weight: bold;
            }}
            .footer {{
                text-align: center;
                margin-top: 50px;
                padding-top: 20px;
                border-top: 1px solid #ddd;
                color: #7f8c8d;
            }}
            .toc {{
                background: #f8f9fa;
                padding: 20px;
                border-radius: 8px;
                margin: 20px 0;
            }}
            .toc ul {{
                list-style-type: none;
                padding-left: 0;
            }}
            .toc li {{
                margin: 8px 0;
            }}
            .toc a {{
                text-decoration: none;
                color: #3498db;
            }}
            .toc a:hover {{
                text-decoration: underline;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>🧬 Comprehensive Metagenomics Analysis Report</h1>
                <p class="subtitle">Assembly-First Pipeline with ML-Enhanced Novelty Detection</p>
                <p><strong>Generated:</strong> {summary['analysis_date']}</p>
            </div>
    """
    
    # Table of Contents
    html_content += """
            <div class="toc">
                <h2>📋 Table of Contents</h2>
                <ul>
                    <li><a href="#executive-summary">Executive Summary</a></li>
                    <li><a href="#dataset-overview">Dataset Overview</a></li>
                    <li><a href="#key-findings">Key Findings</a></li>
                    <li><a href="#risk-assessment">Risk Assessment</a></li>
                    <li><a href="#detailed-results">Detailed Analysis Results</a></li>
                    <li><a href="#recommendations">Recommendations</a></li>
                    <li><a href="#methodology">Methodology</a></li>
                </ul>
            </div>
    """
    
    # Executive Summary
    html_content += f"""
            <div id="executive-summary" class="executive-summary">
                <h2>📊 Executive Summary</h2>
                <div class="metrics-grid">
                    <div class="metric-card">
                        <div class="metric-value">{summary.get('dataset_overview', {}).get('total_samples', 'N/A')}</div>
                        <div class="metric-label">Total Samples</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{summary.get('dataset_overview', {}).get('total_features', 'N/A'):,}</div>
                        <div class="metric-label">Total Features</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{summary.get('key_findings', {}).get('novelty_analysis', {}).get('high_novelty_features', 'N/A')}</div>
                        <div class="metric-label">High Novelty Features</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value risk-{summary.get('risk_assessment', {}).get('risk_level', 'low').lower()}">{summary.get('risk_assessment', {}).get('risk_level', 'LOW')}</div>
                        <div class="metric-label">Risk Level</div>
                    </div>
                </div>
            </div>
    """
    
    # Dataset Overview Section
    dataset_overview = summary.get('dataset_overview', {})
    html_content += f"""
            <div id="dataset-overview" class="section">
                <h2>📈 Dataset Overview</h2>
                <div class="metrics-grid">
                    <div class="metric-card">
                        <div class="metric-value">{dataset_overview.get('total_samples', 'N/A')}</div>
                        <div class="metric-label">Samples Analyzed</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{dataset_overview.get('total_features', 'N/A'):,}</div>
                        <div class="metric-label">Genomic Features</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{dataset_overview.get('total_mags', 'N/A')}</div>
                        <div class="metric-label">Metagenome-Assembled Genomes</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{dataset_overview.get('data_sparsity', 0):.1%}</div>
                        <div class="metric-label">Data Sparsity</div>
                    </div>
                </div>
    """
    
    # Add dataset overview visualization
    if 'dataset_overview' in embedded_images:
        html_content += f"""
                <div class="visualization">
                    <h3>Dataset Composition Analysis</h3>
                    <img src="data:image/png;base64,{embedded_images['dataset_overview']}" alt="Dataset Overview">
                </div>
        """
    
    html_content += "</div>"
    
    # Key Findings Section
    key_findings = summary.get('key_findings', {})
    html_content += """
            <div id="key-findings" class="section">
                <h2>🔍 Key Findings</h2>
    """
    
    # ML Analysis Findings
    ml_findings = key_findings.get('ml_analysis', {})
    if ml_findings:
        html_content += f"""
                <h3>🤖 Machine Learning Analysis</h3>
                <p><strong>Community Embedding Quality:</strong> {ml_findings.get('community_embedding_quality', 'N/A')}</p>
                <p><strong>VAE Reconstruction Correlation:</strong> {ml_findings.get('vae_reconstruction_correlation', 'N/A'):.3f}</p>
        """
    
    # Novelty Analysis Findings
    novelty_findings = key_findings.get('novelty_analysis', {})
    if novelty_findings:
        discovery_rate = novelty_findings.get('novelty_discovery_rate', 0)
        html_content += f"""
                <h3>🧬 Novelty Discovery</h3>
                <p><strong>Novel Features Detected:</strong> {novelty_findings.get('novel_features_detected', 'N/A'):,}</p>
                <p><strong>High Novelty Features:</strong> {novelty_findings.get('high_novelty_features', 'N/A'):,}</p>
                <p><strong>Discovery Rate:</strong> {discovery_rate:.1%}</p>
                
                {"<div class='alert-box alert-high'><strong>Significant Novel Content:</strong> High rate of novel genetic material detected requiring further investigation.</div>" if discovery_rate > 0.1 else ""}
        """
    
    # Temporal Analysis Findings
    temporal_findings = key_findings.get('temporal_analysis', {})
    if temporal_findings:
        html_content += f"""
                <h3>⏰ Temporal Dynamics</h3>
                <p><strong>Significant Trends:</strong> {temporal_findings.get('significant_temporal_trends', 'N/A')}</p>
                <p><strong>Increasing Trends:</strong> {temporal_findings.get('increasing_trends', 'N/A')}</p>
                <p><strong>Decreasing Trends:</strong> {temporal_findings.get('decreasing_trends', 'N/A')}</p>
        """
    
    # Add analysis summary visualization
    if 'analysis_summary' in embedded_images:
        html_content += f"""
                <div class="visualization">
                    <h3>Analysis Results Summary</h3>
                    <img src="data:image/png;base64,{embedded_images['analysis_summary']}" alt="Analysis Summary">
                </div>
        """
    
    html_content += "</div>"
    
    # Risk Assessment Section
    risk_assessment = summary.get('risk_assessment', {})
    html_content += f"""
            <div id="risk-assessment" class="section">
                <h2>🚨 Risk Assessment</h2>
                
                <div class="alert-box alert-{'critical' if risk_assessment.get('risk_level') == 'HIGH' else 'high' if risk_assessment.get('risk_level') == 'MODERATE' else 'medium'}">
                    <strong>Overall Risk Level: {risk_assessment.get('risk_level', 'LOW')}</strong>
                </div>
                
                <div class="metrics-grid">
                    <div class="metric-card">
                        <div class="metric-value">{risk_assessment.get('total_pathogen_detections', 'N/A')}</div>
                        <div class="metric-label">Total Pathogen Detections</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value risk-high">{risk_assessment.get('high_risk_detections', 'N/A')}</div>
                        <div class="metric-label">High-Risk Detections</div>
                    </div>
                </div>
                
                {"<p><strong>⚠️ Action Required:</strong> High-risk pathogen detections require immediate attention and investigation.</p>" if risk_assessment.get('high_risk_detections', 0) > 0 else ""}
            </div>
    """
    
    # Detailed Results Section
    html_content += """
            <div id="detailed-results" class="section">
                <h2>📋 Detailed Analysis Results</h2>
    """
    
    # Add detailed results tables and summaries here
    if results['pathogen_risk'] is not None and not results['pathogen_risk'].empty:
        high_risk_pathogens = results['pathogen_risk'][results['pathogen_risk']['risk_level'] == 'high'].head(10)
        if not high_risk_pathogens.empty:
            html_content += """
                <h3>🦠 High-Risk Pathogen Detections</h3>
                <table class="data-table">
                    <thead>
                        <tr>
                            <th>Sample ID</th>
                            <th>Contig ID</th>
                            <th>Risk Score</th>
                            <th>Pathogen Type</th>
                            <th>Description</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for _, pathogen in high_risk_pathogens.iterrows():
                html_content += f"""
                        <tr>
                            <td>{pathogen.get('sample_id', 'N/A')}</td>
                            <td>{str(pathogen.get('contig_id', 'N/A'))[:20]}...</td>
                            <td class="risk-high">{pathogen.get('final_risk_score', 0):.3f}</td>
                            <td>{pathogen.get('pathogen_type', 'N/A')}</td>
                            <td>{str(pathogen.get('best_hit_description', 'N/A'))[:50]}...</td>
                        </tr>
                """
            
            html_content += """
                    </tbody>
                </table>
            """
    
    html_content += "</div>"
    
    # Recommendations Section
    recommendations = summary.get('recommendations', [])
    html_content += """
            <div id="recommendations" class="section">
                <h2>💡 Recommendations</h2>
                <div class="recommendations">
    """
    
    for rec in recommendations:
        priority_class = rec['priority'].lower()
        html_content += f"""
                    <div class="recommendation-item">
                        <h4><span class="risk-{priority_class}">[{rec['priority']}]</span> {rec['category']}</h4>
                        <p><strong>{rec['recommendation']}</strong></p>
                        <div class="action-items">
                            <strong>Action Items:</strong>
                            <ul>
        """
        
        for action in rec.get('action_items', []):
            html_content += f"<li>{action}</li>"
        
        html_content += """
                            </ul>
                        </div>
                    </div>
        """
    
    html_content += """
                </div>
            </div>
    """
    
    # Methodology Section
    html_content += """
            <div id="methodology" class="section">
                <h2>🔬 Methodology</h2>
                <p>This analysis was performed using a comprehensive two-phase metagenomics pipeline:</p>
                
                <h3>Phase 1: Single-Sample Processing</h3>
                <ul>
                    <li><strong>Quality Control:</strong> FastQC and fastp for read quality assessment and filtering</li>
                    <li><strong>Assembly:</strong> MetaSPAdes for de novo metagenomic assembly</li>
                    <li><strong>Binning:</strong> Multi-algorithm approach for genome binning</li>
                    <li><strong>Taxonomy:</strong> GTDB-Tk for taxonomic classification</li>
                    <li><strong>Novelty Detection:</strong> Dual approach combining homology (BLAST) and ML (DNABert + Isolation Forest)</li>
                </ul>
                
                <h3>Phase 2: Cross-Sample Analysis</h3>
                <ul>
                    <li><strong>MAG Dereplication:</strong> FastANI-based clustering for representative genome selection</li>
                    <li><strong>Community Embedding:</strong> Variational Autoencoder for low-dimensional community representation</li>
                    <li><strong>Anomaly Detection:</strong> Isolation Forest for sample and feature anomaly identification</li>
                    <li><strong>Temporal Analysis:</strong> Time-series analysis for trend detection and changepoint identification</li>
                    <li><strong>Pathogen Screening:</strong> Multi-database BLAST with ML-enhanced risk scoring</li>
                </ul>
                
                <h3>Key Innovations</h3>
                <ul>
                    <li><strong>Assembly-First Approach:</strong> Prioritizes de novo assembly for comprehensive novelty detection</li>
                    <li><strong>Dual Novelty Detection:</strong> Combines traditional homology with modern ML pattern recognition</li>
                    <li><strong>Integrated ML/DL:</strong> VAEs and Isolation Forest enhance traditional bioinformatics</li>
                    <li><strong>Comprehensive Risk Assessment:</strong> Multi-factor scoring for pathogen risk evaluation</li>
                </ul>
            </div>
    """
    
    # Footer
    html_content += f"""
            <div class="footer">
                <p>Generated by Comprehensive Metagenomics Pipeline v2.0</p>
                <p>Report generated on {summary['analysis_date']}</p>
                <p>For questions or support, please refer to the pipeline documentation.</p>
            </div>
        </div>
    </body>
    </html>
    """
    
    # Write HTML file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    logging.info(f"HTML report generated: {output_file}")

def generate_pdf_summary(summary, output_file):
    """Generate executive summary PDF"""
    logging.info("Generating PDF summary")
    
    with PdfPages(output_file) as pdf:
        # Title page
        fig = plt.figure(figsize=(8.5, 11))
        fig.suptitle('Metagenomics Analysis\nExecutive Summary', fontsize=20, fontweight='bold')
        
        # Remove axes
        ax = fig.add_subplot(111)
        ax.axis('off')
        
        # Add summary content
        summary_text = f"""
Analysis Date: {summary['analysis_date']}

DATASET OVERVIEW
• Total Samples: {summary.get('dataset_overview', {}).get('total_samples', 'N/A')}
• Total Features: {summary.get('dataset_overview', {}).get('total_features', 'N/A'):,}
• Data Sparsity: {summary.get('dataset_overview', {}).get('data_sparsity', 0):.1%}

KEY FINDINGS
• Novel Features: {summary.get('key_findings', {}).get('novelty_analysis', {}).get('novel_features_detected', 'N/A'):,}
• High Novelty: {summary.get('key_findings', {}).get('novelty_analysis', {}).get('high_novelty_features', 'N/A'):,}
• Temporal Trends: {summary.get('key_findings', {}).get('temporal_analysis', {}).get('significant_temporal_trends', 'N/A')}

RISK ASSESSMENT
• Overall Risk Level: {summary.get('risk_assessment', {}).get('risk_level', 'LOW')}
• Pathogen Detections: {summary.get('risk_assessment', {}).get('total_pathogen_detections', 'N/A')}
• High-Risk Alerts: {summary.get('risk_assessment', {}).get('high_risk_detections', 'N/A')}

RECOMMENDATIONS
"""
        
        # Add recommendations
        recommendations = summary.get('recommendations', [])
        for i, rec in enumerate(recommendations[:5], 1):  # Top 5 recommendations
            summary_text += f"\n{i}. [{rec['priority']}] {rec['recommendation']}"
        
        ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', fontfamily='monospace')
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    logging.info(f"PDF summary generated: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate comprehensive analysis report')
    parser.add_argument('--abundance-matrix', required=True,
                        help='Abundance matrix file')
    parser.add_argument('--mag-abundance-matrix', required=True,
                        help='MAG abundance matrix file')
    parser.add_argument('--dereplicated-mags', required=True,
                        help='Dereplicated MAGs file')
    parser.add_argument('--vae-embeddings', required=True,
                        help='VAE embeddings file')
    parser.add_argument('--vae-metrics', required=True,
                        help='VAE metrics JSON file')
    parser.add_argument('--isolation-forest-results', required=True,
                        help='Isolation Forest results file')
    parser.add_argument('--pan-novelty-results', required=True,
                        help='Pan-novelty results file')
    parser.add_argument('--novelty-emergence', required=True,
                        help='Novelty emergence patterns file')
    parser.add_argument('--temporal-trends', required=True,
                        help='Temporal trends file')
    parser.add_argument('--pathogen-risk', required=True,
                        help='Pathogen risk assessment file')
    parser.add_argument('--high-risk-alerts', required=True,
                        help='High-risk alerts file')
    parser.add_argument('--metadata', required=True,
                        help='Sample metadata file')
    parser.add_argument('--output-html-report', required=True,
                        help='Output HTML report file')
    parser.add_argument('--output-pdf-summary', required=True,
                        help='Output PDF summary file')
    parser.add_argument('--output-interactive-dashboard', required=True,
                        help='Output interactive dashboard file')
    parser.add_argument('--output-analysis-summary', required=True,
                        help='Output analysis summary JSON file')
    parser.add_argument('--log-file', required=True,
                        help='Log file path')
    
    # Report configuration
    parser.add_argument('--include-executive-summary', action='store_true', default=True,
                        help='Include executive summary')
    parser.add_argument('--include-methodology', action='store_true', default=True,
                        help='Include methodology section')
    parser.add_argument('--include-recommendations', action='store_true', default=True,
                        help='Include recommendations')
    
    args = parser.parse_args()
    
    # Setup logging
    os.makedirs(os.path.dirname(args.log_file), exist_ok=True)
    setup_logging(args.log_file)
    
    try:
        logging.info("Starting comprehensive report generation")
        
        # Collect file paths
        file_paths = {
            'abundance_matrix': args.abundance_matrix,
            'mag_abundance_matrix': args.mag_abundance_matrix,
            'dereplicated_mags': args.dereplicated_mags,
            'vae_embeddings': args.vae_embeddings,
            'vae_metrics': args.vae_metrics,
            'isolation_forest_results': args.isolation_forest_results,
            'pan_novelty_results': args.pan_novelty_results,
            'novelty_emergence': args.novelty_emergence,
            'temporal_trends': args.temporal_trends,
            'pathogen_risk': args.pathogen_risk,
            'high_risk_alerts': args.high_risk_alerts,
            'metadata': args.metadata
        }
        
        # Load analysis results
        results = load_analysis_results(file_paths)
        
        # Generate executive summary
        summary = generate_executive_summary(results)
        
        # Report configuration
        config = {
            'include_executive_summary': args.include_executive_summary,
            'include_methodology': args.include_methodology,
            'include_recommendations': args.include_recommendations
        }
        
        # Create output directories
        for output_file in [args.output_html_report, args.output_pdf_summary,
                           args.output_interactive_dashboard, args.output_analysis_summary]:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Create visualizations
        output_dir = os.path.dirname(args.output_html_report)
        visualization_files = create_summary_visualizations(results, output_dir)
        
        # Generate HTML report
        generate_html_report(results, summary, visualization_files, config, args.output_html_report)
        
        # Generate PDF summary
        generate_pdf_summary(summary, args.output_pdf_summary)
        
        # Save analysis summary JSON
        with open(args.output_analysis_summary, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        
        # Create simple interactive dashboard (placeholder)
        dashboard_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Interactive Dashboard</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; text-align: center; }}
                .redirect {{ background: #f0f8ff; padding: 20px; border-radius: 10px; }}
            </style>
        </head>
        <body>
            <div class="redirect">
                <h1>Interactive Dashboard</h1>
                <p>For full interactive dashboard functionality, please refer to the main HTML report:</p>
                <p><a href="{os.path.basename(args.output_html_report)}">View Comprehensive Report</a></p>
                <p><em>Advanced interactive features will be available in future pipeline versions.</em></p>
            </div>
        </body>
        </html>
        """
        
        with open(args.output_interactive_dashboard, 'w') as f:
            f.write(dashboard_content)
        
        logging.info("Comprehensive report generation completed successfully")
        
        # Print summary
        print("\n" + "="*80)
        print("COMPREHENSIVE REPORT GENERATION SUMMARY")
        print("="*80)
        print(f"Analysis Date: {summary['analysis_date']}")
        print(f"Samples Analyzed: {summary.get('dataset_overview', {}).get('total_samples', 'N/A')}")
        print(f"Features Analyzed: {summary.get('dataset_overview', {}).get('total_features', 'N/A'):,}")
        print(f"Risk Level: {summary.get('risk_assessment', {}).get('risk_level', 'LOW')}")
        print(f"Novel Features: {summary.get('key_findings', {}).get('novelty_analysis', {}).get('high_novelty_features', 'N/A')}")
        print(f"High-Risk Detections: {summary.get('risk_assessment', {}).get('high_risk_detections', 'N/A')}")
        
        print(f"\nOutput Files Generated:")
        print(f"  - HTML Report: {args.output_html_report}")
        print(f"  - PDF Summary: {args.output_pdf_summary}")
        print(f"  - Analysis Summary: {args.output_analysis_summary}")
        print(f"  - Interactive Dashboard: {args.output_interactive_dashboard}")
        
        print("="*80)
        
    except Exception as e:
        logging.error(f"Error in report generation: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()



