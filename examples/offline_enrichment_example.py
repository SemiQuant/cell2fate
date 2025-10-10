#!/usr/bin/env python3
"""
Complete working example of cell2fate pipeline with offline enrichment.

This is a complete, runnable example that follows the Pancreas tutorial structure
but includes offline gene set enrichment analysis.

Based on: https://cell2fate.readthedocs.io/en/latest/notebooks/publication_figures/cell2fate_PancreasWithCC.html#

"""

import cell2fate as c2f
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path


def run_cell2fate_pipeline_with_offline_enrichment(data_path, results_path, data_name='Pancreas_with_cc'):
    """
    Run complete cell2fate pipeline with offline enrichment analysis.
    
    Parameters
    ----------
    data_path : str
        Path to directory containing the data
    results_path : str
        Path to directory for saving results
    data_name : str
        Name of the dataset
    """
    
    print("=== Cell2fate Pipeline with Offline Enrichment ===\n")
    
    # Step 1: Download gene sets for offline use
    print("Step 1: Downloading gene sets for offline use...")
    
    gene_sets_dir = os.path.join(results_path, "gene_sets")
    species = "Mouse"  # Based on the Pancreas tutorial
    
    try:
        downloaded_files = c2f.utils.download_gene_sets(
            output_dir=gene_sets_dir,
            species=species,
            gene_sets=None  # Use default gene sets
        )
        
        print(f"Successfully downloaded {len(downloaded_files)} gene set files:")
        for file_path in downloaded_files:
            print(f"  - {file_path}")
            
    except Exception as e:
        print(f"Error downloading gene sets: {e}")
        print("Continuing with pipeline anyway...\n")
    
    # Step 2: Load and prepare data
    print("\nStep 2: Loading and preparing data...")
    
    # Load the data
    adata = sc.read_h5ad(os.path.join(data_path, f"{data_name}_anndata.h5ad"))
    
    # Remove clusters if needed (empty list means no clusters to remove)
    clusters_to_remove = []
    
    # Get training data
    adata = c2f.utils.get_training_data(
        adata, 
        cells_per_cluster=10**5, 
        cluster_column='clusters',
        remove_clusters=clusters_to_remove,
        min_shared_counts=20, 
        n_var_genes=3000
    )
    
    print(f"Data loaded and processed. Shape: {adata.shape}")
    
    # Step 3: Initialize and train the model
    print("\nStep 3: Initialize and train the model...")
    
    # Setup AnnData
    c2f.Cell2fate_DynamicalModel.setup_anndata(
        adata, 
        spliced_label='spliced', 
        unspliced_label='unspliced'
    )
    
    # Get number of modules
    n_modules = c2f.utils.get_max_modules(adata)
    print(f"Number of modules: {n_modules}")
    
    # Initialize model
    mod = c2f.Cell2fate_DynamicalModel(adata, n_modules=n_modules)
    
    # Train model
    print("Training model...")
    mod.train()
    
    # Step 4: Export posteriors
    print("\nStep 4: Export posteriors...")
    
    adata = mod.export_posterior(adata)
    print("Posteriors exported to AnnData object")
    
    # Step 5: Compute module summary statistics with offline enrichment
    print("\nStep 5: Compute module summary statistics with offline enrichment...")
    
    # Compute module summary statistics
    adata = mod.compute_module_summary_statistics(adata)
    
    # Get top features with offline enrichment
    background = list(adata.var_names)  # Use all genes as background
    
    print("Running offline enrichment analysis...")
    tab, results = mod.get_module_top_features(
        adata=adata,
        background=background,
        species='Mouse',  # or 'Human' depending on your data
        p_adj_cutoff=0.01,
        n_top_genes=100,
        local_gene_sets=gene_sets_dir  # Use offline gene sets
    )
    
    # Display results
    print("\nModule top features with offline enrichment:")
    print(tab)
    
    # Save results
    tab.to_csv(os.path.join(results_path, f"{data_name}_module_top_features.csv"))
    
    # Plot module summary statistics
    mod.plot_module_summary_statistics(
        adata, 
        save=os.path.join(results_path, f"{data_name}_module_summary_stats_plot.pdf")
    )
    
    # Step 6: Additional analyses
    print("\nStep 6: Additional analyses...")
    
    # Compare module activation
    mod.compare_module_activation(
        adata, 
        chosen_modules=[1, 2, 3, 4],
        save=os.path.join(results_path, f"{data_name}_module_activation_comparison.pdf")
    )
    
    # Plot top features
    mod.plot_top_features(
        adata, 
        tab, 
        chosen_modules=[1, 2, 3, 4],
        mode='all genes',
        n_top_features=3,
        save=os.path.join(results_path, f"{data_name}_top_features_plot.pdf")
    )
    
    # Visualize velocities
    mod.plot_velocity_umap_Bergen2020(
        adata,
        use_full_posterior=True,
        save=os.path.join(results_path, f"{data_name}_velocity_umap.pdf")
    )
    
    # Step 7: Save the model and processed data
    print("\nStep 7: Saving model and data...")
    
    # Save model
    mod.save(os.path.join(results_path, f"{data_name}_model"))
    
    # Save processed AnnData
    adata.write(os.path.join(results_path, f"{data_name}_processed.h5ad"))
    
    print("\n=== Pipeline Complete ===")
    print(f"Results saved to: {results_path}")
    print(f"Gene sets available in: {gene_sets_dir}")
    
    return mod, adata, tab, results


def main():
    """Main function with example usage."""
    
    # Example usage - update these paths for your system
    data_path = '/path/to/your/data/'  # Update this path
    results_path = '/path/to/results/'  # Update this path
    data_name = 'Pancreas_with_cc'
    
    # Create results directory if it doesn't exist
    os.makedirs(results_path, exist_ok=True)
    
    print("This is a complete example of the cell2fate pipeline with offline enrichment.")
    print("To run this example:")
    print("1. Update the data_path and results_path variables")
    print("2. Ensure you have the required data file")
    print("3. Uncomment the function call below")
    print("4. Run the script")
    
    # Uncomment the line below when you have actual data:
    # mod, adata, tab, results = run_cell2fate_pipeline_with_offline_enrichment(data_path, results_path, data_name)
    
    print("\nExample code prepared. Update paths and uncomment the function call to run.")


if __name__ == "__main__":
    main()
