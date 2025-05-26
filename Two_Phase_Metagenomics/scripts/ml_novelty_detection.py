#!/usr/bin/env python3
"""
🤖 ML-Enhanced Novelty Detection
Part of the Metagenomics Pipeline - Phase 1

Uses DNABert embeddings and Isolation Forest to detect novel genomic patterns
independent of homology-based approaches.
"""

import argparse
import os
import sys
import logging
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset
from transformers import AutoTokenizer, AutoModel
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import json
from tqdm import tqdm
import pickle

# Suppress warnings
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
        description="ML-based novelty detection using DNABert and Isolation Forest"
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
        "--model-name",
        default="zhihan1996/DNA_bert_6",
        help="DNABert model name from HuggingFace"
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=500,
        help="Minimum sequence length for analysis"
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=10000,
        help="Maximum sequence length for analysis"
    )
    parser.add_argument(
        "--max-seq-length",
        type=int,
        default=512,
        help="Maximum sequence length for DNABert"
    )
    parser.add_argument(
        "--overlap-size",
        type=int,
        default=256,
        help="Overlap size for sequence windowing"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Batch size for DNABert inference"
    )
    parser.add_argument(
        "--contamination",
        type=float,
        default=0.1,
        help="Expected contamination rate for Isolation Forest"
    )
    parser.add_argument(
        "--n-estimators",
        type=int,
        default=200,
        help="Number of estimators for Isolation Forest"
    )
    parser.add_argument(
        "--device",
        default="auto",
        help="Device to use: 'cpu', 'cuda', or 'auto'"
    )
    parser.add_argument(
        "--sample-name",
        default="sample",
        help="Sample name for output files"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for CPU operations"
    )
    parser.add_argument(
        "--save-embeddings",
        action="store_true",
        help="Save embeddings for downstream analysis"
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

class DNASequenceDataset(Dataset):
    """Dataset class for DNA sequences"""
    
    def __init__(self, sequences, tokenizer, max_length=512):
        self.sequences = sequences
        self.tokenizer = tokenizer
        self.max_length = max_length
    
    def __len__(self):
        return len(self.sequences)
    
    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        
        # Tokenize sequence
        encoding = self.tokenizer(
            sequence,
            truncation=True,
            padding='max_length',
            max_length=self.max_length,
            return_tensors='pt'
        )
        
        return {
            'input_ids': encoding['input_ids'].flatten(),
            'attention_mask': encoding['attention_mask'].flatten()
        }

class DNABertEmbedder:
    """DNABert-based sequence embedding generator"""
    
    def __init__(self, model_name, device='auto', max_length=512):
        self.model_name = model_name
        self.max_length = max_length
        
        # Set device
        if device == 'auto':
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = torch.device(device)
        
        logging.info(f"Using device: {self.device}")
        
        # Load tokenizer and model
        try:
            logging.info(f"Loading DNABert model: {model_name}")
            self.tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
            self.model = AutoModel.from_pretrained(model_name, trust_remote_code=True)
            self.model.to(self.device)
            self.model.eval()
            
        except Exception as e:
            logging.error(f"Error loading DNABert model: {e}")
            logging.info("Falling back to mock embeddings for testing")
            self.tokenizer = None
            self.model = None
    
    def preprocess_sequence(self, sequence):
        """Preprocess DNA sequence for DNABert"""
        # Convert to uppercase and add spaces between nucleotides
        sequence = sequence.upper()
        # DNABert expects space-separated nucleotides
        spaced_sequence = ' '.join(list(sequence))
        return spaced_sequence
    
    def window_sequence(self, sequence, window_size, overlap_size):
        """Create overlapping windows from long sequences"""
        windows = []
        step = window_size - overlap_size
        
        for i in range(0, len(sequence) - window_size + 1, step):
            window = sequence[i:i + window_size]
            windows.append(window)
        
        # Add final window if sequence doesn't divide evenly
        if len(sequence) % step != 0:
            windows.append(sequence[-window_size:])
        
        return windows
    
    def generate_embeddings(self, sequences, batch_size=32):
        """Generate embeddings for a list of sequences"""
        if self.model is None:
            # Return mock embeddings for testing
            logging.warning("Using mock embeddings - model not loaded")
            return np.random.randn(len(sequences), 768)
        
        # Preprocess sequences
        processed_sequences = [self.preprocess_sequence(seq) for seq in sequences]
        
        # Create dataset and dataloader
        dataset = DNASequenceDataset(processed_sequences, self.tokenizer, self.max_length)
        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
        
        embeddings = []
        
        with torch.no_grad():
            for batch in tqdm(dataloader, desc="Generating embeddings"):
                input_ids = batch['input_ids'].to(self.device)
                attention_mask = batch['attention_mask'].to(self.device)
                
                outputs = self.model(input_ids=input_ids, attention_mask=attention_mask)
                
                # Use pooled output or mean of last hidden states
                if hasattr(outputs, 'pooler_output') and outputs.pooler_output is not None:
                    batch_embeddings = outputs.pooler_output.cpu().numpy()
                else:
                    # Use mean pooling of last hidden states
                    hidden_states = outputs.last_hidden_state
                    # Apply attention mask for proper mean pooling
                    attention_mask_expanded = attention_mask.unsqueeze(-1).expand(hidden_states.size())
                    masked_hidden_states = hidden_states * attention_mask_expanded
                    sum_hidden_states = torch.sum(masked_hidden_states, dim=1)
                    sum_attention_mask = torch.sum(attention_mask, dim=1, keepdim=True)
                    batch_embeddings = (sum_hidden_states / sum_attention_mask).cpu().numpy()
                
                embeddings.append(batch_embeddings)
        
        return np.vstack(embeddings)

class NoveltyDetector:
    """Isolation Forest-based novelty detector"""
    
    def __init__(self, contamination=0.1, n_estimators=200, random_state=42, n_jobs=-1):
        self.contamination = contamination
        self.n_estimators = n_estimators
        self.random_state = random_state
        self.n_jobs = n_jobs
        
        self.isolation_forest = IsolationForest(
            contamination=contamination,
            n_estimators=n_estimators,
            random_state=random_state,
            n_jobs=n_jobs
        )
        
        self.scaler = StandardScaler()
        self.pca = None
        self.is_fitted = False
    
    def fit_predict(self, embeddings, use_pca=True, n_components=50):
        """Fit the model and predict novelty scores"""
        logging.info(f"Training Isolation Forest on {embeddings.shape[0]} samples")
        
        # Standardize embeddings
        embeddings_scaled = self.scaler.fit_transform(embeddings)
        
        # Optional PCA for dimensionality reduction
        if use_pca and embeddings_scaled.shape[1] > n_components:
            self.pca = PCA(n_components=n_components, random_state=self.random_state)
            embeddings_scaled = self.pca.fit_transform(embeddings_scaled)
            logging.info(f"Applied PCA: {embeddings.shape[1]} -> {embeddings_scaled.shape[1]} dimensions")
        
        # Fit and predict
        predictions = self.isolation_forest.fit_predict(embeddings_scaled)
        anomaly_scores = self.isolation_forest.decision_function(embeddings_scaled)
        
        self.is_fitted = True
        
        # Convert to novelty scores (higher = more novel)
        novelty_scores = self.convert_to_novelty_scores(anomaly_scores, predictions)
        
        return novelty_scores, predictions, anomaly_scores
    
    def convert_to_novelty_scores(self, anomaly_scores, predictions):
        """Convert isolation forest scores to novelty scores (0-100)"""
        # Invert scores (more negative = more anomalous = more novel)
        inverted_scores = -anomaly_scores
        
        # Normalize to 0-100 scale
        min_score = np.min(inverted_scores)
        max_score = np.max(inverted_scores)
        
        if max_score > min_score:
            novelty_scores = ((inverted_scores - min_score) / (max_score - min_score)) * 100
        else:
            novelty_scores = np.full_like(inverted_scores, 50.0)
        
        return novelty_scores

class MLNoveltyAnalyzer:
    """Main class for ML-based novelty analysis"""
    
    def __init__(self, config):
        self.config = config
        self.embedder = DNABertEmbedder(
            config['model_name'],
            config['device'],
            config['max_seq_length']
        )
        self.detector = NoveltyDetector(
            contamination=config['contamination'],
            n_estimators=config['n_estimators'],
            n_jobs=config['threads']
        )
        
        self.sequences = {}
        self.embeddings = None
        self.results = {}
    
    def load_sequences(self, assembly_file):
        """Load sequences from assembly file"""
        logging.info(f"Loading sequences from {assembly_file}")
        
        sequences = []
        sequence_ids = []
        
        with open(assembly_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_len = len(record.seq)
                
                if self.config['min_length'] <= seq_len <= self.config['max_length']:
                    sequence = str(record.seq)
                    
                    # Window long sequences
                    if seq_len > self.config['max_seq_length']:
                        windows = self.embedder.window_sequence(
                            sequence,
                            self.config['max_seq_length'],
                            self.config['overlap_size']
                        )
                        
                        for i, window in enumerate(windows):
                            window_id = f"{record.id}_window_{i}"
                            sequences.append(window)
                            sequence_ids.append(window_id)
                            
                            self.sequences[window_id] = {
                                'original_id': record.id,
                                'sequence': window,
                                'length': len(window),
                                'window_index': i,
                                'total_windows': len(windows)
                            }
                    else:
                        sequences.append(sequence)
                        sequence_ids.append(record.id)
                        
                        self.sequences[record.id] = {
                            'original_id': record.id,
                            'sequence': sequence,
                            'length': seq_len,
                            'window_index': 0,
                            'total_windows': 1
                        }
        
        logging.info(f"Loaded {len(sequences)} sequence segments for analysis")
        return sequences, sequence_ids
    
    def analyze_novelty(self, assembly_file):
        """Perform complete novelty analysis"""
        # Load sequences
        sequences, sequence_ids = self.load_sequences(assembly_file)
        
        if not sequences:
            logging.error("No sequences found for analysis")
            return
        
        # Generate embeddings
        logging.info("Generating DNABert embeddings")
        self.embeddings = self.embedder.generate_embeddings(
            sequences, 
            self.config['batch_size']
        )
        
        # Detect novelty
        logging.info("Detecting novel patterns with Isolation Forest")
        novelty_scores, predictions, anomaly_scores = self.detector.fit_predict(self.embeddings)
        
        # Compile results
        self.results = {
            'sequence_ids': sequence_ids,
            'novelty_scores': novelty_scores,
            'predictions': predictions,
            'anomaly_scores': anomaly_scores,
            'sequences': sequences
        }
        
        # Calculate statistics
        self.calculate_statistics()
        
        return self.results
    
    def calculate_statistics(self):
        """Calculate novelty detection statistics"""
        novelty_scores = self.results['novelty_scores']
        predictions = self.results['predictions']
        
        # Calculate score statistics
        stats = {
            'total_sequences': len(novelty_scores),
            'novel_sequences': np.sum(predictions == -1),
            'normal_sequences': np.sum(predictions == 1),
            'novelty_percentage': (np.sum(predictions == -1) / len(predictions)) * 100,
            'mean_novelty_score': np.mean(novelty_scores),
            'median_novelty_score': np.median(novelty_scores),
            'std_novelty_score': np.std(novelty_scores),
            'max_novelty_score': np.max(novelty_scores),
            'min_novelty_score': np.min(novelty_scores),
            'high_novelty_count': np.sum(novelty_scores > 75),
            'medium_novelty_count': np.sum((novelty_scores > 50) & (novelty_scores <= 75)),
            'low_novelty_count': np.sum(novelty_scores <= 50)
        }
        
        self.results['statistics'] = stats
        
        logging.info(f"Novelty analysis complete:")
        logging.info(f"  Total sequences: {stats['total_sequences']}")
        logging.info(f"  Novel sequences: {stats['novel_sequences']} ({stats['novelty_percentage']:.1f}%)")
        logging.info(f"  Mean novelty score: {stats['mean_novelty_score']:.1f}")
    
    def aggregate_window_results(self):
        """Aggregate results for windowed sequences"""
        logging.info("Aggregating results for windowed sequences")
        
        contig_results = {}
        
        for i, seq_id in enumerate(self.results['sequence_ids']):
            seq_info = self.sequences[seq_id]
            original_id = seq_info['original_id']
            
            if original_id not in contig_results:
                contig_results[original_id] = {
                    'window_scores': [],
                    'window_predictions': [],
                    'total_length': 0,
                    'total_windows': seq_info['total_windows']
                }
            
            contig_results[original_id]['window_scores'].append(self.results['novelty_scores'][i])
            contig_results[original_id]['window_predictions'].append(self.results['predictions'][i])
            contig_results[original_id]['total_length'] += seq_info['length']
        
        # Calculate aggregate scores
        aggregated_results = []
        
        for contig_id, data in contig_results.items():
            scores = np.array(data['window_scores'])
            predictions = np.array(data['window_predictions'])
            
            result = {
                'contig_id': contig_id,
                'mean_novelty_score': np.mean(scores),
                'max_novelty_score': np.max(scores),
                'median_novelty_score': np.median(scores),
                'std_novelty_score': np.std(scores),
                'novel_windows': np.sum(predictions == -1),
                'total_windows': len(scores),
                'novelty_percentage': (np.sum(predictions == -1) / len(predictions)) * 100,
                'total_length': data['total_length']
            }
            
            aggregated_results.append(result)
        
        return pd.DataFrame(aggregated_results)

def create_visualization_plots(analyzer, output_dir, sample_name):
    """Create visualization plots"""
    logging.info("Creating visualization plots")
    
    # Create plots directory
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    results = analyzer.results
    
    # 1. Novelty score distribution
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Score histogram
    novelty_scores = results['novelty_scores']
    predictions = results['predictions']
    
    ax1.hist(novelty_scores, bins=50, alpha=0.7, edgecolor='black')
    ax1.axvline(np.mean(novelty_scores), color='red', linestyle='--', 
                label=f'Mean: {np.mean(novelty_scores):.1f}')
    ax1.set_xlabel('ML Novelty Score')
    ax1.set_ylabel('Frequency')
    ax1.set_title(f'{sample_name}: ML Novelty Score Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Score by prediction
    novel_scores = novelty_scores[predictions == -1]
    normal_scores = novelty_scores[predictions == 1]
    
    ax2.hist([normal_scores, novel_scores], bins=30, alpha=0.7, 
             label=['Normal', 'Novel'], color=['blue', 'red'])
    ax2.set_xlabel('ML Novelty Score')
    ax2.set_ylabel('Frequency')
    ax2.set_title(f'{sample_name}: Scores by Prediction')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Score vs sequence length
    seq_lengths = [analyzer.sequences[seq_id]['length'] 
                  for seq_id in results['sequence_ids']]
    
    colors = ['red' if pred == -1 else 'blue' for pred in predictions]
    ax3.scatter(seq_lengths, novelty_scores, c=colors, alpha=0.6, s=20)
    ax3.set_xlabel('Sequence Length (bp)')
    ax3.set_ylabel('ML Novelty Score')
    ax3.set_title(f'{sample_name}: Novelty vs Length')
    ax3.grid(True, alpha=0.3)
    
    # Add legend
    novel_patch = plt.Line2D([0], [0], marker='o', color='w', 
                           markerfacecolor='red', markersize=8, label='Novel')
    normal_patch = plt.Line2D([0], [0], marker='o', color='w', 
                            markerfacecolor='blue', markersize=8, label='Normal')
    ax3.legend(handles=[novel_patch, normal_patch])
    
    # Embedding PCA plot
    if analyzer.embeddings is not None and analyzer.detector.pca is not None:
        # Use existing PCA transformation
        embeddings_pca = analyzer.detector.pca.transform(
            analyzer.detector.scaler.transform(analyzer.embeddings)
        )
        
        scatter = ax4.scatter(embeddings_pca[:, 0], embeddings_pca[:, 1], 
                            c=novelty_scores, cmap='viridis', alpha=0.6, s=20)
        ax4.set_xlabel('PC1')
        ax4.set_ylabel('PC2')
        ax4.set_title(f'{sample_name}: Embedding Space (PCA)')
        plt.colorbar(scatter, ax=ax4, label='Novelty Score')
    else:
        ax4.text(0.5, 0.5, 'PCA not available', ha='center', va='center', 
                transform=ax4.transAxes)
        ax4.set_title('Embedding Visualization (N/A)')
    
    plt.tight_layout()
    plt.savefig(plots_dir / f"{sample_name}_ml_novelty_analysis.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Statistics summary plot
    stats = results['statistics']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Prediction counts
    categories = ['Normal', 'Novel']
    counts = [stats['normal_sequences'], stats['novel_sequences']]
    colors = ['blue', 'red']
    
    wedges, texts, autotexts = ax1.pie(counts, labels=categories, colors=colors, 
                                      autopct='%1.1f%%', startangle=90)
    ax1.set_title(f'{sample_name}: Sequence Classification')
    
    # Novelty score categories
    score_categories = ['Low\n(≤50)', 'Medium\n(50-75)', 'High\n(>75)']
    score_counts = [
        stats['low_novelty_count'],
        stats['medium_novelty_count'],
        stats['high_novelty_count']
    ]
    
    ax2.bar(score_categories, score_counts, color=['green', 'orange', 'red'], 
            alpha=0.7, edgecolor='black')
    ax2.set_ylabel('Number of Sequences')
    ax2.set_title(f'{sample_name}: Novelty Score Categories')
    ax2.grid(True, alpha=0.3)
    
    # Add counts on bars
    for i, count in enumerate(score_counts):
        ax2.text(i, count + max(score_counts) * 0.01, str(count), 
                ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(plots_dir / f"{sample_name}_ml_novelty_summary.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Plots saved to {plots_dir}")

def save_results(analyzer, output_dir, sample_name, save_embeddings=False):
    """Save analysis results"""
    logging.info("Saving results")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    results = analyzer.results
    
    # Save detailed results
    detailed_results = []
    for i, seq_id in enumerate(results['sequence_ids']):
        seq_info = analyzer.sequences[seq_id]
        
        result = {
            'sequence_id': seq_id,
            'original_contig_id': seq_info['original_id'],
            'window_index': seq_info['window_index'],
            'total_windows': seq_info['total_windows'],
            'sequence_length': seq_info['length'],
            'ml_novelty_score': results['novelty_scores'][i],
            'ml_prediction': 'Novel' if results['predictions'][i] == -1 else 'Normal',
            'anomaly_score': results['anomaly_scores'][i]
        }
        detailed_results.append(result)
    
    detailed_df = pd.DataFrame(detailed_results)
    detailed_file = output_path / f"{sample_name}_ml_novelty_detailed.tsv"
    detailed_df.to_csv(detailed_file, sep='\t', index=False)
    
    # Save aggregated results (if sequences were windowed)
    aggregated_df = analyzer.aggregate_window_results()
    aggregated_file = output_path / f"{sample_name}_ml_novelty_contigs.tsv"
    aggregated_df.to_csv(aggregated_file, sep='\t', index=False)
    
    # Save statistics
    stats_file = output_path / f"{sample_name}_ml_novelty_stats.json"
    with open(stats_file, 'w') as f:
        json.dump(results['statistics'], f, indent=2)
    
    # Save embeddings if requested
    if save_embeddings and analyzer.embeddings is not None:
        embeddings_file = output_path / f"{sample_name}_ml_embeddings.npy"
        np.save(embeddings_file, analyzer.embeddings)
        
        # Save sequence mapping
        seq_mapping = {
            'sequence_ids': results['sequence_ids'],
            'original_contigs': [analyzer.sequences[seq_id]['original_id'] 
                               for seq_id in results['sequence_ids']]
        }
        mapping_file = output_path / f"{sample_name}_embedding_mapping.json"
        with open(mapping_file, 'w') as f:
            json.dump(seq_mapping, f, indent=2)
    
    # Save model if trained
    if analyzer.detector.is_fitted:
        model_file = output_path / f"{sample_name}_isolation_forest_model.pkl"
        with open(model_file, 'wb') as f:
            pickle.dump(analyzer.detector, f)
    
    logging.info(f"Results saved to {output_path}")
    
    # Print summary
    stats = results['statistics']
    print(f"\n🤖 ML Novelty Detection Summary for {sample_name}")
    print("=" * 50)
    print(f"Total sequences analyzed: {stats['total_sequences']:,}")
    print(f"Novel sequences detected: {stats['novel_sequences']:,} ({stats['novelty_percentage']:.1f}%)")
    print(f"Mean novelty score: {stats['mean_novelty_score']:.1f}")
    print(f"High novelty sequences (>75): {stats['high_novelty_count']:,}")
    print(f"Medium novelty sequences (50-75): {stats['medium_novelty_count']:,}")
    print(f"Low novelty sequences (≤50): {stats['low_novelty_count']:,}")

def main():
    """Main function"""
    args = parse_arguments()
    setup_logging(args.log_level)
    
    logging.info("Starting ML-based novelty detection analysis")
    
    # Prepare configuration
    config = {
        'model_name': args.model_name,
        'device': args.device,
        'min_length': args.min_length,
        'max_length': args.max_length,
        'max_seq_length': args.max_seq_length,
        'overlap_size': args.overlap_size,
        'batch_size': args.batch_size,
        'contamination': args.contamination,
        'n_estimators': args.n_estimators,
        'threads': args.threads
    }
    
    # Initialize analyzer
    analyzer = MLNoveltyAnalyzer(config)
    
    # Perform analysis
    results = analyzer.analyze_novelty(args.assembly)
    
    if results:
        # Save results
        save_results(analyzer, args.output, args.sample_name, args.save_embeddings)
        
        # Generate plots if requested
        if args.plots:
            create_visualization_plots(analyzer, args.output, args.sample_name)
    
    logging.info("ML novelty detection analysis completed successfully")

if __name__ == "__main__":
    main()


