#!/usr/bin/env python3
"""
Pathogen Detection and Risk Assessment
Comprehensive pathogen screening with ML-enhanced risk scoring
"""

import sys
import os
import pandas as pd
import numpy as np
import json
import logging
import subprocess
import tempfile
import shutil
from pathlib import Path
from collections import defaultdict
import argparse

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
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

def load_pathogen_databases(database_config):
    """Load and validate pathogen databases"""
    logging.info("Loading pathogen databases")
    
    databases = {}
    
    for db_name, db_path in database_config.items():
        if db_path and os.path.exists(db_path):
            databases[db_name] = db_path
            logging.info(f"Loaded {db_name} database: {db_path}")
        else:
            logging.warning(f"Database not found: {db_name} at {db_path}")
    
    if not databases:
        logging.warning("No pathogen databases available")
    
    return databases

def run_blast_analysis(contigs_files, databases, blast_params, temp_dir):
    """Run BLAST analysis against pathogen databases"""
    logging.info("Running BLAST analysis against pathogen databases")
    
    blast_results = {}
    
    # Combine all contigs into a single file
    combined_contigs = os.path.join(temp_dir, "all_contigs.fasta")
    sample_contig_mapping = {}
    
    with open(combined_contigs, 'w') as outfile:
        for contig_file in contigs_files:
            if os.path.exists(contig_file):
                sample_name = extract_sample_name_from_path(contig_file)
                
                with open(contig_file, 'r') as infile:
                    current_contig = None
                    for line in infile:
                        if line.startswith('>'):
                            # Modify header to include sample name
                            contig_id = line.strip()[1:]
                            new_header = f">{sample_name}_{contig_id}"
                            sample_contig_mapping[f"{sample_name}_{contig_id}"] = {
                                'sample_id': sample_name,
                                'original_contig_id': contig_id
                            }
                            outfile.write(new_header + '\n')
                            current_contig = f"{sample_name}_{contig_id}"
                        else:
                            outfile.write(line)
    
    # Run BLAST against each database
    for db_name, db_path in databases.items():
        logging.info(f"Running BLAST against {db_name}")
        
        blast_output = os.path.join(temp_dir, f"blast_{db_name}.txt")
        
        blast_cmd = [
            "blastn",
            "-query", combined_contigs,
            "-db", db_path,
            "-out", blast_output,
            "-evalue", str(blast_params.get('evalue', 1e-5)),
            "-max_target_seqs", str(blast_params.get('max_targets', 10)),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
            "-num_threads", "8"
        ]
        
        try:
            result = subprocess.run(blast_cmd, capture_output=True, text=True, timeout=3600)
            
            if result.returncode == 0:
                blast_results[db_name] = parse_blast_output(
                    blast_output, 
                    sample_contig_mapping,
                    blast_params
                )
                logging.info(f"BLAST against {db_name} completed successfully")
            else:
                logging.error(f"BLAST against {db_name} failed: {result.stderr}")
                
        except subprocess.TimeoutExpired:
            logging.error(f"BLAST against {db_name} timed out")
        except Exception as e:
            logging.error(f"Error running BLAST against {db_name}: {str(e)}")
    
    return blast_results

def extract_sample_name_from_path(file_path):
    """Extract sample name from contig file path"""
    # Assume structure: .../sample_name/assembly/contigs.fasta
    path_parts = Path(file_path).parts
    for i, part in enumerate(path_parts):
        if part == 'assembly' and i > 0:
            return path_parts[i-1]
    
    # Fallback: use parent directory name
    return Path(file_path).parent.parent.name

def parse_blast_output(blast_file, sample_contig_mapping, blast_params):
    """Parse BLAST output and filter by quality thresholds"""
    logging.info(f"Parsing BLAST output: {blast_file}")
    
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']
    
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        return pd.DataFrame(columns=columns + ['sample_id', 'original_contig_id', 'coverage'])
    
    blast_df = pd.read_csv(blast_file, sep='\t', names=columns)
    
    # Add sample information
    blast_df['sample_id'] = blast_df['qseqid'].map(
        lambda x: sample_contig_mapping.get(x, {}).get('sample_id', 'unknown')
    )
    blast_df['original_contig_id'] = blast_df['qseqid'].map(
        lambda x: sample_contig_mapping.get(x, {}).get('original_contig_id', x)
    )
    
    # Calculate coverage
    blast_df['coverage'] = (blast_df['qend'] - blast_df['qstart'] + 1) / blast_df['length'] * 100
    
    # Apply quality filters
    identity_threshold = blast_params.get('identity_threshold', 70)
    coverage_threshold = blast_params.get('coverage_threshold', 50)
    
    filtered_df = blast_df[
        (blast_df['pident'] >= identity_threshold) &
        (blast_df['coverage'] >= coverage_threshold)
    ].copy()
    
    logging.info(f"BLAST hits before filtering: {len(blast_df)}")
    logging.info(f"BLAST hits after filtering: {len(filtered_df)}")
    
    return filtered_df

def classify_pathogen_hits(blast_results, known_pathogens_file):
    """Classify BLAST hits as potential pathogens"""
    logging.info("Classifying pathogen hits")
    
    # Load known pathogens if available
    known_pathogens = set()
    if known_pathogens_file and os.path.exists(known_pathogens_file):
        with open(known_pathogens_file, 'r') as f:
            known_pathogens = set(line.strip().lower() for line in f if line.strip())
        logging.info(f"Loaded {len(known_pathogens)} known pathogen identifiers")
    
    pathogen_classifications = {}
    
    for db_name, blast_df in blast_results.items():
        if blast_df.empty:
            continue
        
        classifications = []
        
        for _, hit in blast_df.iterrows():
            hit_title = hit['stitle'].lower()
            
            # Pathogen classification logic
            pathogen_score = 0
            pathogen_type = 'unknown'
            confidence = 'low'
            
            # Check against known pathogen list
            if any(pathogen in hit_title for pathogen in known_pathogens):
                pathogen_score += 0.8
                confidence = 'high'
            
            # Database-specific classification
            if db_name == 'virulence_factors':
                pathogen_score += 0.6
                pathogen_type = 'virulent'
            elif db_name == 'antibiotic_resistance':
                pathogen_score += 0.5
                pathogen_type = 'resistant'
            elif db_name == 'pathogen_genomes':
                pathogen_score += 0.7
                pathogen_type = 'pathogenic'
            
            # Quality-based scoring
            identity_score = hit['pident'] / 100
            coverage_score = min(hit['coverage'] / 100, 1.0)
            quality_score = (identity_score + coverage_score) / 2
            
            pathogen_score = pathogen_score * quality_score
            
            # Adjust confidence based on quality
            if quality_score > 0.8:
                if confidence != 'high':
                    confidence = 'medium'
            else:
                confidence = 'low'
            
            # Keyword-based classification
            pathogen_keywords = {
                'virus': ['virus', 'viral', 'phage'],
                'bacteria': ['bacteria', 'bacterial', 'bacteri'],
                'fungus': ['fungus', 'fungi', 'yeast'],
                'parasite': ['parasite', 'parasitic', 'protozoa']
            }
            
            organism_type = 'unknown'
            for org_type, keywords in pathogen_keywords.items():
                if any(keyword in hit_title for keyword in keywords):
                    organism_type = org_type
                    break
            
            classifications.append({
                'sample_id': hit['sample_id'],
                'contig_id': hit['original_contig_id'],
                'database': db_name,
                'hit_id': hit['sseqid'],
                'hit_description': hit['stitle'],
                'identity': hit['pident'],
                'coverage': hit['coverage'],
                'evalue': hit['evalue'],
                'pathogen_score': pathogen_score,
                'pathogen_type': pathogen_type,
                'organism_type': organism_type,
                'confidence': confidence
            })
        
        pathogen_classifications[db_name] = pd.DataFrame(classifications)
    
    return pathogen_classifications

def calculate_risk_scores(pathogen_classifications, abundance_data, novelty_data, risk_weights):
    """Calculate comprehensive risk scores for detected pathogens"""
    logging.info("Calculating pathogen risk scores")
    
    # Combine all pathogen hits
    all_hits = []
    for db_name, df in pathogen_classifications.items():
        if not df.empty:
            all_hits.append(df)
    
    if not all_hits:
        return pd.DataFrame()
    
    combined_hits = pd.concat(all_hits, ignore_index=True)
    
    # Load abundance data
    abundance_df = None
    if abundance_data and os.path.exists(abundance_data):
        abundance_df = pd.read_csv(abundance_data, sep='\t', index_col=0)
    
    # Load novelty data
    novelty_df = None
    if novelty_data and os.path.exists(novelty_data):
        novelty_df = pd.read_csv(novelty_data, sep='\t', index_col=0)
    
    risk_scores = []
    
    # Group by sample and contig
    for (sample_id, contig_id), group in combined_hits.groupby(['sample_id', 'contig_id']):
        
        # Base pathogen score (highest score from multiple hits)
        max_pathogen_score = group['pathogen_score'].max()
        
        # Abundance score
        abundance_score = 0
        if abundance_df is not None:
            if contig_id in abundance_df.index and sample_id in abundance_df.columns:
                abundance = abundance_df.loc[contig_id, sample_id]
                # Normalize abundance score (log scale)
                abundance_score = min(np.log10(abundance + 1) / 6, 1.0)  # Cap at log10(1M)
        
        # Novelty score
        novelty_score = 0
        if novelty_df is not None:
            if contig_id in novelty_df.index:
                novelty_score = novelty_df.loc[contig_id, 'combined_novelty_score'] if 'combined_novelty_score' in novelty_df.columns else 0
        
        # Virulence and resistance scores
        virulence_score = 0
        resistance_score = 0
        
        for _, hit in group.iterrows():
            if hit['database'] == 'virulence_factors':
                virulence_score = max(virulence_score, hit['pathogen_score'])
            elif hit['database'] == 'antibiotic_resistance':
                resistance_score = max(resistance_score, hit['pathogen_score'])
        
        # Calculate weighted risk score
        weighted_risk_score = (
            risk_weights.get('novelty_weight', 0.4) * novelty_score +
            risk_weights.get('abundance_weight', 0.3) * abundance_score +
            risk_weights.get('virulence_weight', 0.2) * virulence_score +
            risk_weights.get('resistance_weight', 0.1) * resistance_score
        )
        
        # Final risk score combines pathogen score and weighted factors
        final_risk_score = 0.6 * max_pathogen_score + 0.4 * weighted_risk_score
        
        # Determine risk level
        risk_level = 'low'
        if final_risk_score >= 0.8:
            risk_level = 'high'
        elif final_risk_score >= 0.6:
            risk_level = 'medium'
        
        # Get best hit information
        best_hit = group.loc[group['pathogen_score'].idxmax()]
        
        risk_scores.append({
            'sample_id': sample_id,
            'contig_id': contig_id,
            'pathogen_score': max_pathogen_score,
            'abundance_score': abundance_score,
            'novelty_score': novelty_score,
            'virulence_score': virulence_score,
            'resistance_score': resistance_score,
            'weighted_risk_score': weighted_risk_score,
            'final_risk_score': final_risk_score,
            'risk_level': risk_level,
            'best_hit_database': best_hit['database'],
            'best_hit_description': best_hit['hit_description'],
            'best_hit_identity': best_hit['identity'],
            'best_hit_coverage': best_hit['coverage'],
            'pathogen_type': best_hit['pathogen_type'],
            'organism_type': best_hit['organism_type'],
            'confidence': best_hit['confidence'],
            'num_database_hits': len(group)
        })
    
    risk_df = pd.DataFrame(risk_scores)
    
    if not risk_df.empty:
        risk_df = risk_df.sort_values('final_risk_score', ascending=False)
        logging.info(f"Calculated risk scores for {len(risk_df)} potential pathogen detections")
    
    return risk_df

def generate_pathogen_alerts(risk_scores, alert_thresholds):
    """Generate pathogen alerts based on risk thresholds"""
    logging.info("Generating pathogen alerts")
    
    alerts = {
        'high_risk': [],
        'medium_risk': [],
        'low_risk': []
    }
    
    if risk_scores.empty:
        return alerts
    
    for _, pathogen in risk_scores.iterrows():
        risk_score = pathogen['final_risk_score']
        
        alert_info = {
            'sample_id': pathogen['sample_id'],
            'contig_id': pathogen['contig_id'],
            'risk_score': risk_score,
            'pathogen_type': pathogen['pathogen_type'],
            'organism_type': pathogen['organism_type'],
            'description': pathogen['best_hit_description'],
            'confidence': pathogen['confidence'],
            'alert_timestamp': pd.Timestamp.now().isoformat()
        }
        
        if risk_score >= alert_thresholds.get('high_risk', 0.8):
            alerts['high_risk'].append(alert_info)
        elif risk_score >= alert_thresholds.get('medium_risk', 0.6):
            alerts['medium_risk'].append(alert_info)
        else:
            alerts['low_risk'].append(alert_info)
    
    logging.info(f"Generated alerts: {len(alerts['high_risk'])} high-risk, "
                f"{len(alerts['medium_risk'])} medium-risk, {len(alerts['low_risk'])} low-risk")
    
    return alerts

def create_pathogen_timeline(risk_scores, metadata):
    """Create temporal timeline of pathogen detections"""
    logging.info("Creating pathogen detection timeline")
    
    if risk_scores.empty:
        return pd.DataFrame()
    
    # Load metadata
    metadata_df = pd.read_csv(metadata, sep='\t')
    
    if 'collection_date' not in metadata_df.columns:
        logging.warning("No collection_date in metadata for timeline analysis")
        return pd.DataFrame()
    
    # Merge with metadata
    timeline_data = risk_scores.merge(
        metadata_df[['sample_id', 'collection_date']], 
        on='sample_id', 
        how='left'
    )
    
    # Convert to datetime
    timeline_data['collection_date'] = pd.to_datetime(timeline_data['collection_date'])
    
    # Sort by date
    timeline_data = timeline_data.sort_values('collection_date')
    
    # Aggregate by date and risk level
    timeline_summary = timeline_data.groupby(['collection_date', 'risk_level']).agg({
        'contig_id': 'count',
        'final_risk_score': 'mean'
    }).reset_index()
    
    timeline_summary.columns = ['collection_date', 'risk_level', 'detection_count', 'mean_risk_score']
    
    return timeline_summary

def create_pathogen_visualizations(risk_scores, pathogen_alerts, timeline_data, output_dir):
    """Create pathogen detection visualizations"""
    logging.info("Creating pathogen detection visualizations")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 1. Risk score distribution
    if not risk_scores.empty:
        axes[0, 0].hist(risk_scores['final_risk_score'], bins=20, alpha=0.7, edgecolor='black')
        axes[0, 0].set_title('Distribution of Pathogen Risk Scores')
        axes[0, 0].set_xlabel('Risk Score')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].axvline(0.6, color='orange', linestyle='--', label='Medium Risk')
        axes[0, 0].axvline(0.8, color='red', linestyle='--', label='High Risk')
        axes[0, 0].legend()
    
    # 2. Risk level distribution
    alert_counts = [len(pathogen_alerts['high_risk']), 
                   len(pathogen_alerts['medium_risk']), 
                   len(pathogen_alerts['low_risk'])]
    risk_labels = ['High Risk', 'Medium Risk', 'Low Risk']
    colors = ['red', 'orange', 'yellow']
    
    axes[0, 1].pie(alert_counts, labels=risk_labels, colors=colors, autopct='%1.1f%%')
    axes[0, 1].set_title('Distribution of Risk Levels')
    
    # 3. Pathogen types
    if not risk_scores.empty:
        pathogen_counts = risk_scores['pathogen_type'].value_counts()
        axes[0, 2].bar(pathogen_counts.index, pathogen_counts.values)
        axes[0, 2].set_title('Pathogen Types Detected')
        axes[0, 2].set_xlabel('Pathogen Type')
        axes[0, 2].set_ylabel('Count')
        axes[0, 2].tick_params(axis='x', rotation=45)
    
    # 4. Organism types
    if not risk_scores.empty:
        organism_counts = risk_scores['organism_type'].value_counts()
        axes[1, 0].bar(organism_counts.index, organism_counts.values)
        axes[1, 0].set_title('Organism Types Detected')
        axes[1, 0].set_xlabel('Organism Type')
        axes[1, 0].set_ylabel('Count')
        axes[1, 0].tick_params(axis='x', rotation=45)
    
    # 5. Risk components correlation
    if not risk_scores.empty and len(risk_scores) > 5:
        risk_components = risk_scores[['pathogen_score', 'abundance_score', 'novelty_score', 
                                     'virulence_score', 'resistance_score']].corr()
        
        sns.heatmap(risk_components, annot=True, cmap='coolwarm', center=0, 
                   square=True, ax=axes[1, 1])
        axes[1, 1].set_title('Risk Component Correlations')
    
    # 6. Timeline if available
    if not timeline_data.empty:
        for risk_level in ['high', 'medium', 'low']:
            level_data = timeline_data[timeline_data['risk_level'] == risk_level]
            if not level_data.empty:
                axes[1, 2].plot(level_data['collection_date'], level_data['detection_count'], 
                               'o-', label=f'{risk_level.title()} Risk', alpha=0.7)
        
        axes[1, 2].set_title('Pathogen Detections Over Time')
        axes[1, 2].set_xlabel('Date')
        axes[1, 2].set_ylabel('Detection Count')
        axes[1, 2].legend()
        axes[1, 2].tick_params(axis='x', rotation=45)
    else:
        axes[1, 2].text(0.5, 0.5, 'No temporal data available', 
                       ha='center', va='center', transform=axes[1, 2].transAxes)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pathogen_detection_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Pathogen detection and risk assessment')
    parser.add_argument('--abundance-matrix', required=True,
                        help='Abundance matrix file')
    parser.add_argument('--novelty-data', required=True,
                        help='Pan-novelty analysis results file')
    parser.add_argument('--contigs', nargs='+', required=True,
                        help='List of contig files from samples')
    parser.add_argument('--metadata', required=True,
                        help='Sample metadata file')
    parser.add_argument('--output-risk-assessment', required=True,
                        help='Output pathogen risk assessment file')
    parser.add_argument('--output-high-risk-alerts', required=True,
                        help='Output high-risk alerts file')
    parser.add_argument('--output-virulence-analysis', required=True,
                        help='Output virulence factor analysis file')
    parser.add_argument('--output-resistance-analysis', required=True,
                        help='Output antibiotic resistance analysis file')
    parser.add_argument('--output-pathogen-timeline', required=True,
                        help='Output pathogen detection timeline file')
    parser.add_argument('--log-file', required=True,
                        help='Log file path')
    
    # Database paths
    parser.add_argument('--virulence-db', 
                        help='Virulence factors database path')
    parser.add_argument('--resistance-db',
                        help='Antibiotic resistance database path')
    parser.add_argument('--pathogen-db',
                        help='Pathogen genomes database path')
    parser.add_argument('--known-pathogens',
                        help='Known pathogens list file')
    
    # Analysis parameters
    parser.add_argument('--evalue', type=float, default=1e-5,
                        help='BLAST e-value threshold')
    parser.add_argument('--identity-threshold', type=float, default=70,
                        help='Minimum identity percentage')
    parser.add_argument('--coverage-threshold', type=float, default=50,
                        help='Minimum coverage percentage')
    parser.add_argument('--novelty-weight', type=float, default=0.4,
                        help='Weight for novelty in risk scoring')
    parser.add_argument('--abundance-weight', type=float, default=0.3,
                        help='Weight for abundance in risk scoring')
    parser.add_argument('--virulence-weight', type=float, default=0.2,
                        help='Weight for virulence in risk scoring')
    parser.add_argument('--resistance-weight', type=float, default=0.1,
                        help='Weight for resistance in risk scoring')
    
    args = parser.parse_args()
    
    # Setup logging
    os.makedirs(os.path.dirname(args.log_file), exist_ok=True)
    setup_logging(args.log_file)
    
    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix="pathogen_detection_")
    
    try:
        logging.info("Starting pathogen detection analysis")
        
        # Configure databases
        database_config = {
            'virulence_factors': args.virulence_db,
            'antibiotic_resistance': args.resistance_db,
            'pathogen_genomes': args.pathogen_db
        }
        
        databases = load_pathogen_databases(database_config)
        
        if not databases:
            logging.warning("No pathogen databases available. Skipping BLAST analysis.")
            # Create empty output files
            for output_file in [args.output_risk_assessment, args.output_high_risk_alerts,
                               args.output_virulence_analysis, args.output_resistance_analysis,
                               args.output_pathogen_timeline]:
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                pd.DataFrame().to_csv(output_file, sep='\t', index=False)
            return
        
        # BLAST parameters
        blast_params = {
            'evalue': args.evalue,
            'identity_threshold': args.identity_threshold,
            'coverage_threshold': args.coverage_threshold,
            'max_targets': 10
        }
        
        # Run BLAST analysis
        blast_results = run_blast_analysis(args.contigs, databases, blast_params, temp_dir)
        
        # Classify pathogen hits
        pathogen_classifications = classify_pathogen_hits(blast_results, args.known_pathogens)
        
        # Risk scoring weights
        risk_weights = {
            'novelty_weight': args.novelty_weight,
            'abundance_weight': args.abundance_weight,
            'virulence_weight': args.virulence_weight,
            'resistance_weight': args.resistance_weight
        }
        
        # Calculate risk scores
        risk_scores = calculate_risk_scores(
            pathogen_classifications,
            args.abundance_matrix,
            args.novelty_data,
            risk_weights
        )
        
        # Generate alerts
        alert_thresholds = {'high_risk': 0.8, 'medium_risk': 0.6, 'low_risk': 0.3}
        pathogen_alerts = generate_pathogen_alerts(risk_scores, alert_thresholds)
        
        # Create pathogen timeline
        timeline_data = create_pathogen_timeline(risk_scores, args.metadata)
        
        # Create output directories
        for output_file in [args.output_risk_assessment, args.output_high_risk_alerts,
                           args.output_virulence_analysis, args.output_resistance_analysis,
                           args.output_pathogen_timeline]:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Save results
        logging.info("Saving pathogen detection results")
        
        # Risk assessment
        if not risk_scores.empty:
            risk_scores.to_csv(args.output_risk_assessment, sep='\t', index=False)
        else:
            pd.DataFrame().to_csv(args.output_risk_assessment, sep='\t', index=False)
        
        # High-risk alerts
        high_risk_df = pd.DataFrame(pathogen_alerts['high_risk'])
        high_risk_df.to_csv(args.output_high_risk_alerts, sep='\t', index=False)
        
        # Virulence analysis
        if 'virulence_factors' in pathogen_classifications:
            pathogen_classifications['virulence_factors'].to_csv(
                args.output_virulence_analysis, sep='\t', index=False
            )
        else:
            pd.DataFrame().to_csv(args.output_virulence_analysis, sep='\t', index=False)
        
        # Resistance analysis
        if 'antibiotic_resistance' in pathogen_classifications:
            pathogen_classifications['antibiotic_resistance'].to_csv(
                args.output_resistance_analysis, sep='\t', index=False
            )
        else:
            pd.DataFrame().to_csv(args.output_resistance_analysis, sep='\t', index=False)
        
        # Timeline
        timeline_data.to_csv(args.output_pathogen_timeline, sep='\t', index=False)
        
        # Create visualizations
        output_dir = os.path.dirname(args.output_risk_assessment)
        create_pathogen_visualizations(risk_scores, pathogen_alerts, timeline_data, output_dir)
        
        logging.info("Pathogen detection analysis completed successfully")
        
        # Print summary
        print("\n" + "="*80)
        print("PATHOGEN DETECTION SUMMARY")
        print("="*80)
        print(f"Databases used: {', '.join(databases.keys())}")
        print(f"Samples analyzed: {len(set(args.contigs))}")
        
        if not risk_scores.empty:
            print(f"Potential pathogen detections: {len(risk_scores)}")
            print(f"High-risk alerts: {len(pathogen_alerts['high_risk'])}")
            print(f"Medium-risk alerts: {len(pathogen_alerts['medium_risk'])}")
            print(f"Low-risk alerts: {len(pathogen_alerts['low_risk'])}")
            
            print(f"\nRisk score statistics:")
            print(f"  Mean: {risk_scores['final_risk_score'].mean():.3f}")
            print(f"  Max: {risk_scores['final_risk_score'].max():.3f}")
            print(f"  Std: {risk_scores['final_risk_score'].std():.3f}")
            
            if len(pathogen_alerts['high_risk']) > 0:
                print(f"\nHigh-risk detections:")
                for alert in pathogen_alerts['high_risk'][:5]:  # Show top 5
                    print(f"  - {alert['sample_id']}: {alert['description'][:50]}... (Risk: {alert['risk_score']:.3f})")
        else:
            print("No potential pathogens detected")
        
        print("="*80)
        
    except Exception as e:
        logging.error(f"Error in pathogen detection: {str(e)}")
        sys.exit(1)
    
    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == "__main__":
    main()

