#!/usr/bin/env python3
"""
VAE Community Embedding for Microbiome Data
Variational Autoencoder for learning low-dimensional representations of microbial communities
"""

import sys
import os
import pandas as pd
import numpy as np
import json
import logging
import pickle
import argparse
from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler
from sklearn.model_selection import train_test_split
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

class MicrobiomeVAE(nn.Module):
    """Variational Autoencoder for microbiome abundance data"""
    
    def __init__(self, input_dim, latent_dim, hidden_dims, dropout_rate=0.2):
        super(MicrobiomeVAE, self).__init__()
        
        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.hidden_dims = hidden_dims
        
        # Encoder
        encoder_layers = []
        prev_dim = input_dim
        
        for hidden_dim in hidden_dims:
            encoder_layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout_rate)
            ])
            prev_dim = hidden_dim
        
        self.encoder = nn.Sequential(*encoder_layers)
        
        # Latent space
        self.fc_mu = nn.Linear(prev_dim, latent_dim)
        self.fc_logvar = nn.Linear(prev_dim, latent_dim)
        
        # Decoder
        decoder_layers = []
        prev_dim = latent_dim
        
        for hidden_dim in reversed(hidden_dims):
            decoder_layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout_rate)
            ])
            prev_dim = hidden_dim
        
        decoder_layers.append(nn.Linear(prev_dim, input_dim))
        self.decoder = nn.Sequential(*decoder_layers)
        
    def encode(self, x):
        """Encode input to latent distribution parameters"""
        h = self.encoder(x)
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        return mu, logvar
    
    def reparameterize(self, mu, logvar):
        """Reparameterization trick"""
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std
    
    def decode(self, z):
        """Decode latent vector to reconstruction"""
        return self.decoder(z)
    
    def forward(self, x):
        """Forward pass through VAE"""
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z)
        return recon, mu, logvar, z

def vae_loss_function(recon_x, x, mu, logvar, beta=1.0):
    """VAE loss function (reconstruction + KL divergence)"""
    # Reconstruction loss (MSE)
    recon_loss = F.mse_loss(recon_x, x, reduction='sum')
    
    # KL divergence loss
    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    
    # Beta-VAE formulation
    total_loss = recon_loss + beta * kl_loss
    
    return total_loss, recon_loss, kl_loss

def preprocess_abundance_data(abundance_matrix, metadata, config):
    """Preprocess abundance data for VAE training"""
    logging.info("Preprocessing abundance data for VAE")
    
    # Log transform if specified
    if config.get('log_transform', True):
        # Add pseudocount to avoid log(0)
        pseudocount = 1e-6
        abundance_processed = np.log(abundance_matrix + pseudocount)
        logging.info("Applied log transformation with pseudocount")
    else:
        abundance_processed = abundance_matrix.copy()
    
    # Filter features by prevalence
    prevalence_threshold = config.get('prevalence_threshold', 0.1)
    prevalence = (abundance_matrix > config.get('min_abundance_threshold', 0.0001)).mean(axis=1)
    high_prevalence_features = prevalence >= prevalence_threshold
    
    abundance_filtered = abundance_processed.loc[high_prevalence_features]
    logging.info(f"Filtered features by prevalence: {len(abundance_filtered)} features retained")
    
    # Transpose for sample-wise analysis (samples as rows, features as columns)
    abundance_transposed = abundance_filtered.T
    
    # Handle missing values
    abundance_transposed = abundance_transposed.fillna(0)
    
    # Scale features
    scaler = StandardScaler()
    abundance_scaled = scaler.fit_transform(abundance_transposed)
    
    logging.info(f"Final data shape for VAE: {abundance_scaled.shape}")
    
    return abundance_scaled, abundance_transposed.index, abundance_filtered.index, scaler

def create_data_loaders(X_train, X_val, batch_size):
    """Create PyTorch data loaders"""
    train_dataset = TensorDataset(torch.FloatTensor(X_train))
    val_dataset = TensorDataset(torch.FloatTensor(X_val))
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    
    return train_loader, val_loader

def train_vae(model, train_loader, val_loader, config, device):
    """Train the VAE model"""
    logging.info("Starting VAE training")
    
    optimizer = optim.Adam(model.parameters(), lr=config['learning_rate'])
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=10, factor=0.8)
    
    num_epochs = config['num_epochs']
    beta = config['beta']
    early_stopping_patience = config['early_stopping_patience']
    
    train_history = {
        'train_loss': [],
        'train_recon_loss': [],
        'train_kl_loss': [],
        'val_loss': [],
        'val_recon_loss': [],
        'val_kl_loss': []
    }
    
    best_val_loss = float('inf')
    patience_counter = 0
    
    for epoch in range(num_epochs):
        # Training phase
        model.train()
        train_loss = 0
        train_recon_loss = 0
        train_kl_loss = 0
        
        for batch_idx, (data,) in enumerate(train_loader):
            data = data.to(device)
            optimizer.zero_grad()
            
            recon_data, mu, logvar, _ = model(data)
            loss, recon_loss, kl_loss = vae_loss_function(recon_data, data, mu, logvar, beta)
            
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
            train_recon_loss += recon_loss.item()
            train_kl_loss += kl_loss.item()
        
        # Validation phase
        model.eval()
        val_loss = 0
        val_recon_loss = 0
        val_kl_loss = 0
        
        with torch.no_grad():
            for data, in val_loader:
                data = data.to(device)
                recon_data, mu, logvar, _ = model(data)
                loss, recon_loss, kl_loss = vae_loss_function(recon_data, data, mu, logvar, beta)
                
                val_loss += loss.item()
                val_recon_loss += recon_loss.item()
                val_kl_loss += kl_loss.item()
        
        # Calculate average losses
        avg_train_loss = train_loss / len(train_loader.dataset)
        avg_train_recon = train_recon_loss / len(train_loader.dataset)
        avg_train_kl = train_kl_loss / len(train_loader.dataset)
        
        avg_val_loss = val_loss / len(val_loader.dataset)
        avg_val_recon = val_recon_loss / len(val_loader.dataset)
        avg_val_kl = val_kl_loss / len(val_loader.dataset)
        
        # Store history
        train_history['train_loss'].append(avg_train_loss)
        train_history['train_recon_loss'].append(avg_train_recon)
        train_history['train_kl_loss'].append(avg_train_kl)
        train_history['val_loss'].append(avg_val_loss)
        train_history['val_recon_loss'].append(avg_val_recon)
        train_history['val_kl_loss'].append(avg_val_kl)
        
        # Learning rate scheduling
        scheduler.step(avg_val_loss)
        
        # Early stopping
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            patience_counter = 0
        else:
            patience_counter += 1
        
        if epoch % 10 == 0:
            logging.info(f"Epoch {epoch}/{num_epochs}: Train Loss={avg_train_loss:.4f}, Val Loss={avg_val_loss:.4f}")
        
        if patience_counter >= early_stopping_patience:
            logging.info(f"Early stopping at epoch {epoch}")
            break
    
    logging.info("VAE training completed")
    return train_history

def generate_embeddings(model, data, device):
    """Generate latent embeddings for all samples"""
    model.eval()
    embeddings = []
    reconstruction_errors = []
    
    with torch.no_grad():
        data_tensor = torch.FloatTensor(data).to(device)
        recon_data, mu, logvar, z = model(data_tensor)
        
        embeddings = z.cpu().numpy()
        
        # Calculate reconstruction error for each sample
        recon_errors = F.mse_loss(recon_data, data_tensor, reduction='none').mean(dim=1)
        reconstruction_errors = recon_errors.cpu().numpy()
    
    return embeddings, reconstruction_errors

def calculate_model_metrics(model, data_loader, device, beta=1.0):
    """Calculate comprehensive model evaluation metrics"""
    model.eval()
    total_loss = 0
    total_recon_loss = 0
    total_kl_loss = 0
    num_samples = 0
    
    all_reconstructions = []
    all_originals = []
    
    with torch.no_grad():
        for data, in data_loader:
            data = data.to(device)
            recon_data, mu, logvar, _ = model(data)
            
            loss, recon_loss, kl_loss = vae_loss_function(recon_data, data, mu, logvar, beta)
            
            total_loss += loss.item()
            total_recon_loss += recon_loss.item()
            total_kl_loss += kl_loss.item()
            num_samples += data.size(0)
            
            all_reconstructions.append(recon_data.cpu().numpy())
            all_originals.append(data.cpu().numpy())
    
    all_reconstructions = np.vstack(all_reconstructions)
    all_originals = np.vstack(all_originals)
    
    # Calculate correlation between original and reconstructed data
    correlations = []
    for i in range(all_originals.shape[0]):
        corr = np.corrcoef(all_originals[i], all_reconstructions[i])[0, 1]
        if not np.isnan(corr):
            correlations.append(corr)
    
    metrics = {
        'average_total_loss': total_loss / num_samples,
        'average_reconstruction_loss': total_recon_loss / num_samples,
        'average_kl_loss': total_kl_loss / num_samples,
        'mean_reconstruction_correlation': np.mean(correlations),
        'std_reconstruction_correlation': np.std(correlations),
        'num_samples_evaluated': num_samples
    }
    
    return metrics

def plot_training_history(train_history, output_dir):
    """Plot training history"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    # Total loss
    axes[0, 0].plot(train_history['train_loss'], label='Train')
    axes[0, 0].plot(train_history['val_loss'], label='Validation')
    axes[0, 0].set_title('Total Loss')
    axes[0, 0].set_xlabel('Epoch')
    axes[0, 0].set_ylabel('Loss')
    axes[0, 0].legend()
    
    # Reconstruction loss
    axes[0, 1].plot(train_history['train_recon_loss'], label='Train')
    axes[0, 1].plot(train_history['val_recon_loss'], label='Validation')
    axes[0, 1].set_title('Reconstruction Loss')
    axes[0, 1].set_xlabel('Epoch')
    axes[0, 1].set_ylabel('Loss')
    axes[0, 1].legend()
    
    # KL loss
    axes[1, 0].plot(train_history['train_kl_loss'], label='Train')
    axes[1, 0].plot(train_history['val_kl_loss'], label='Validation')
    axes[1, 0].set_title('KL Divergence Loss')
    axes[1, 0].set_xlabel('Epoch')
    axes[1, 0].set_ylabel('Loss')
    axes[1, 0].legend()
    
    # Learning curves (smoothed)
    window = 5
    if len(train_history['train_loss']) > window:
        train_smooth = pd.Series(train_history['train_loss']).rolling(window).mean()
        val_smooth = pd.Series(train_history['val_loss']).rolling(window).mean()
        axes[1, 1].plot(train_smooth, label='Train (smoothed)')
        axes[1, 1].plot(val_smooth, label='Validation (smoothed)')
        axes[1, 1].set_title('Smoothed Learning Curves')
        axes[1, 1].set_xlabel('Epoch')
        axes[1, 1].set_ylabel('Loss')
        axes[1, 1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'training_history.png'), dpi=300, bbox_inches='tight')
    plt.close()

def visualize_embeddings(embeddings, sample_ids, metadata, output_dir):
    """Create visualizations of the learned embeddings"""
    try:
        from sklearn.decomposition import PCA
        from sklearn.manifold import TSNE
        import umap
        
        # PCA
        pca = PCA(n_components=2)
        pca_coords = pca.fit_transform(embeddings)
        
        # t-SNE
        tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(embeddings)-1))
        tsne_coords = tsne.fit_transform(embeddings)
        
        # UMAP
        umap_reducer = umap.UMAP(n_components=2, random_state=42)
        umap_coords = umap_reducer.fit_transform(embeddings)
        
        # Create plots
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        # PCA plot
        axes[0].scatter(pca_coords[:, 0], pca_coords[:, 1], alpha=0.7)
        axes[0].set_title(f'PCA (Explained Variance: {pca.explained_variance_ratio_.sum():.3f})')
        axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.3f})')
        axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.3f})')
        
        # t-SNE plot
        axes[1].scatter(tsne_coords[:, 0], tsne_coords[:, 1], alpha=0.7)
        axes[1].set_title('t-SNE')
        axes[1].set_xlabel('t-SNE 1')
        axes[1].set_ylabel('t-SNE 2')
        
        # UMAP plot
        axes[2].scatter(umap_coords[:, 0], umap_coords[:, 1], alpha=0.7)
        axes[2].set_title('UMAP')
        axes[2].set_xlabel('UMAP 1')
        axes[2].set_ylabel('UMAP 2')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'embedding_visualizations.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
    except ImportError as e:
        logging.warning(f"Could not create embedding visualizations: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='VAE community embedding for microbiome data')
    parser.add_argument('--abundance-matrix', required=True,
                        help='Abundance matrix file')
    parser.add_argument('--metadata', required=True,
                        help='Sample metadata file')
    parser.add_argument('--output-embeddings', required=True,
                        help='Output embeddings file (npy format)')
    parser.add_argument('--output-model', required=True,
                        help='Output trained model file')
    parser.add_argument('--output-history', required=True,
                        help='Output training history JSON')
    parser.add_argument('--output-metrics', required=True,
                        help='Output model metrics JSON')
    parser.add_argument('--output-reconstruction-errors', required=True,
                        help='Output reconstruction errors file')
    parser.add_argument('--log-file', required=True,
                        help='Log file path')
    
    # VAE hyperparameters
    parser.add_argument('--latent-dim', type=int, default=32,
                        help='Latent dimension size')
    parser.add_argument('--hidden-dims', nargs='+', type=int, default=[256, 128, 64],
                        help='Hidden layer dimensions')
    parser.add_argument('--batch-size', type=int, default=32,
                        help='Batch size for training')
    parser.add_argument('--learning-rate', type=float, default=0.001,
                        help='Learning rate')
    parser.add_argument('--num-epochs', type=int, default=500,
                        help='Number of training epochs')
    parser.add_argument('--early-stopping-patience', type=int, default=50,
                        help='Early stopping patience')
    parser.add_argument('--beta', type=float, default=1.0,
                        help='Beta parameter for beta-VAE')
    parser.add_argument('--dropout-rate', type=float, default=0.2,
                        help='Dropout rate')
    parser.add_argument('--validation-split', type=float, default=0.2,
                        help='Validation split ratio')
    parser.add_argument('--random-seed', type=int, default=42,
                        help='Random seed')
    
    # Data preprocessing parameters
    parser.add_argument('--log-transform', action='store_true', default=True,
                        help='Apply log transformation')
    parser.add_argument('--min-abundance-threshold', type=float, default=0.0001,
                        help='Minimum abundance threshold')
    parser.add_argument('--prevalence-threshold', type=float, default=0.1,
                        help='Minimum prevalence threshold')
    
    args = parser.parse_args()
    
    # Setup logging
    os.makedirs(os.path.dirname(args.log_file), exist_ok=True)
    setup_logging(args.log_file)
    
    # Set random seeds
    torch.manual_seed(args.random_seed)
    np.random.seed(args.random_seed)
    
    # Check for GPU availability
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f"Using device: {device}")
    
    try:
        logging.info("Starting VAE community embedding analysis")
        
        # Load data
        logging.info(f"Loading abundance matrix from: {args.abundance_matrix}")
        abundance_matrix = pd.read_csv(args.abundance_matrix, sep='\t', index_col=0)
        
        logging.info(f"Loading metadata from: {args.metadata}")
        metadata = pd.read_csv(args.metadata, sep='\t')
        
        # Preprocessing configuration
        preprocessing_config = {
            'log_transform': args.log_transform,
            'min_abundance_threshold': args.min_abundance_threshold,
            'prevalence_threshold': args.prevalence_threshold
        }
        
        # Preprocess data
        processed_data, sample_ids, feature_ids, scaler = preprocess_abundance_data(
            abundance_matrix, metadata, preprocessing_config
        )
        
        # Split data
        X_train, X_val = train_test_split(
            processed_data, 
            test_size=args.validation_split, 
            random_state=args.random_seed
        )
        
        # Create data loaders
        train_loader, val_loader = create_data_loaders(X_train, X_val, args.batch_size)
        
        # Initialize model
        model = MicrobiomeVAE(
            input_dim=processed_data.shape[1],
            latent_dim=args.latent_dim,
            hidden_dims=args.hidden_dims,
            dropout_rate=args.dropout_rate
        ).to(device)
        
        logging.info(f"VAE architecture: {processed_data.shape[1]} -> {args.hidden_dims} -> {args.latent_dim}")
        
        # Training configuration
        training_config = {
            'learning_rate': args.learning_rate,
            'num_epochs': args.num_epochs,
            'early_stopping_patience': args.early_stopping_patience,
            'beta': args.beta
        }
        
        # Train model
        train_history = train_vae(model, train_loader, val_loader, training_config, device)
        
        # Generate embeddings for all samples
        embeddings, reconstruction_errors = generate_embeddings(model, processed_data, device)
        
        # Calculate model metrics
        all_loader = DataLoader(
            TensorDataset(torch.FloatTensor(processed_data)), 
            batch_size=args.batch_size, 
            shuffle=False
        )
        model_metrics = calculate_model_metrics(model, all_loader, device, args.beta)
        
        # Create output directories
        for output_file in [args.output_embeddings, args.output_model, args.output_history, 
                           args.output_metrics, args.output_reconstruction_errors]:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Save embeddings
        np.save(args.output_embeddings, embeddings)
        
        # Save model
        model_state = {
            'model_state_dict': model.state_dict(),
            'model_config': {
                'input_dim': processed_data.shape[1],
                'latent_dim': args.latent_dim,
                'hidden_dims': args.hidden_dims,
                'dropout_rate': args.dropout_rate
            },
            'preprocessing_config': preprocessing_config,
            'scaler': scaler,
            'feature_ids': feature_ids,
            'sample_ids': sample_ids
        }
        with open(args.output_model, 'wb') as f:
            pickle.dump(model_state, f)
        
        # Save training history
        with open(args.output_history, 'w') as f:
            json.dump(train_history, f, indent=2)
        
        # Save model metrics
        with open(args.output_metrics, 'w') as f:
            json.dump(model_metrics, f, indent=2)
        
        # Save reconstruction errors
        reconstruction_df = pd.DataFrame({
            'sample_id': sample_ids,
            'reconstruction_error': reconstruction_errors
        })
        reconstruction_df.to_csv(args.output_reconstruction_errors, sep='\t', index=False)
        
        # Create visualizations
        output_dir = os.path.dirname(args.output_embeddings)
        plot_training_history(train_history, output_dir)
        visualize_embeddings(embeddings, sample_ids, metadata, output_dir)
        
        logging.info("VAE community embedding analysis completed successfully")
        
        # Print summary
        print("\n" + "="*80)
        print("VAE COMMUNITY EMBEDDING SUMMARY")
        print("="*80)
        print(f"Input dimensions: {processed_data.shape}")
        print(f"Latent dimensions: {args.latent_dim}")
        print(f"Training samples: {len(X_train)}")
        print(f"Validation samples: {len(X_val)}")
        print(f"Final training loss: {train_history['train_loss'][-1]:.4f}")
        print(f"Final validation loss: {train_history['val_loss'][-1]:.4f}")
        print(f"Mean reconstruction correlation: {model_metrics['mean_reconstruction_correlation']:.4f}")
        print(f"Embeddings shape: {embeddings.shape}")
        print("="*80)
        
    except Exception as e:
        logging.error(f"Error in VAE community embedding: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()

