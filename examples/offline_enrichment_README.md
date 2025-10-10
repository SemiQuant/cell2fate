# Cell2fate Offline Enrichment Analysis

This directory contains examples for using cell2fate with offline gene set enrichment analysis.

## Overview

The offline enrichment functionality allows you to:
1. Download gene sets from Enrichr for offline use
2. Perform enrichment analysis without internet access
3. Use the same interface as online enrichment

## Quick Start

### 1. Download Gene Sets

Using the command-line interface:

```bash
# Download default gene sets for Human
cell2fate download-genesets --species Human

# Download default gene sets for Mouse
cell2fate download-genesets --species Mouse

# Download to custom directory
cell2fate download-genesets --output-dir /path/to/gene_sets
```

Or programmatically:

```python
from cell2fate.utils import download_gene_sets

# Download gene sets
downloaded_files = download_gene_sets(
    output_dir='gene_sets',
    species='Mouse',
    gene_sets=None  # Use default gene sets
)
```

### 2. Use Offline Enrichment in Full Pipeline

Following the [Pancreas tutorial](https://cell2fate.readthedocs.io/en/latest/notebooks/publication_figures/cell2fate_PancreasWithCC.html#):

```python
import cell2fate as c2f
import scanpy as sc
import os

# Step 1: Download gene sets for offline use
gene_sets_dir = "gene_sets"
species = "Mouse"  # or "Human"

downloaded_files = c2f.utils.download_gene_sets(
    output_dir=gene_sets_dir,
    species=species,
    gene_sets=None
)

# Step 2: Load and prepare data
adata = sc.read_h5ad("your_data.h5ad")
adata = c2f.utils.get_training_data(
    adata, 
    cells_per_cluster=10**5, 
    cluster_column='clusters',
    remove_clusters=[],
    min_shared_counts=20, 
    n_var_genes=3000
)

# Step 3: Initialize and train the model
c2f.Cell2fate_DynamicalModel.setup_anndata(
    adata, 
    spliced_label='spliced', 
    unspliced_label='unspliced'
)

n_modules = c2f.utils.get_max_modules(adata)
mod = c2f.Cell2fate_DynamicalModel(adata, n_modules=n_modules)
mod.train()

# Step 4: Export posteriors
adata = mod.export_posterior(adata)

# Step 5: Compute module summary statistics with offline enrichment
adata = mod.compute_module_summary_statistics(adata)

# Get top features with offline enrichment
background = list(adata.var_names)
tab, results = mod.get_module_top_features(
    adata=adata,
    background=background,
    species='Mouse',
    p_adj_cutoff=0.01,
    n_top_genes=100,
    local_gene_sets=gene_sets_dir  # Use offline gene sets
)

# Display results
print("Module top features with offline enrichment:")
print(tab)

# Plot module summary statistics
mod.plot_module_summary_statistics(adata, save="module_summary_stats_plot.pdf")

# Step 6: Additional analyses
mod.compare_module_activation(
    adata, 
    chosen_modules=[1, 2, 3, 4],
    save="module_activation_comparison.pdf"
)

mod.plot_top_features(
    adata, 
    tab, 
    chosen_modules=[1, 2, 3, 4],
    mode='all genes',
    n_top_features=3,
    save="top_features_plot.pdf"
)

mod.plot_velocity_umap_Bergen2020(
    adata,
    use_full_posterior=True,
    save="velocity_umap.pdf"
)
```

## Files

- `offline_enrichment_example.py`: Basic example script
- `cell2fate_pipeline_with_offline_enrichment.py`: Pipeline integration example
- `complete_pipeline_example.py`: Complete working pipeline example
- `README.md`: This documentation file

## Gene Set Files

The downloaded gene sets are stored as GMT (Gene Matrix Transpose) files:
- `GO_Biological_Process_2021.gmt`: Gene Ontology Biological Process terms
- `GO_Cellular_Component_2021.gmt`: Gene Ontology Cellular Component terms  
- `KEGG_2021_Human.gmt`: KEGG pathways (Human)
- `KEGG_2019_Mouse.gmt`: KEGG pathways (Mouse)

## Notes

- The offline enrichment analysis uses hypergeometric testing with FDR correction
- Results are compatible with the online Enrichr API format
- Gene sets are downloaded from the Maayan Lab Enrichr database
- Requires internet access only for initial download
