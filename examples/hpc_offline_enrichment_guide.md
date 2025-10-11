# Cell2fate Offline Enrichment for HPC Environments

## Quick Start

**Error on compute node?** `ConnectionError: Failed to resolve 'maayanlab.cloud'`

**Fix in 2 steps:**

```bash
# 1. On login node (with internet):
cell2fate download-genesets --species Mouse --output-dir gene_sets
```

```python
# 2. In your compute job script, add this parameter:
tab, results = mod.get_module_top_features(
    adata=adata,
    background=list(adata.var_names),
    species='Mouse',
    local_gene_sets='gene_sets'  # ← Add this
)
```

That's it! Read below for full details.

---

## Problem

When running cell2fate on HPC compute nodes, you may encounter connection errors like:
```
requests.exceptions.ConnectionError: HTTPConnectionPool(host='maayanlab.cloud', port=80): Max retries exceeded
```

This happens because compute nodes typically don't have internet access for security reasons.

## Solution: Two-Step Process

### Step 1: Download Gene Sets (on login node with internet)

Before submitting your job, download the gene sets on a login node that has internet access:

```python
#!/usr/bin/env python3
"""
Run this script on a login node BEFORE submitting your compute job.
This downloads the gene sets that will be used offline on compute nodes.
"""

import cell2fate as c2f
import os

# Define where to save gene sets (accessible from compute nodes)
gene_sets_dir = "/path/to/your/project/gene_sets"  # Update this path
species = "Mouse"  # or "Human"

# Create directory if it doesn't exist
os.makedirs(gene_sets_dir, exist_ok=True)

# Download gene sets
print(f"Downloading {species} gene sets to {gene_sets_dir}...")
downloaded_files = c2f.utils.download_gene_sets(
    output_dir=gene_sets_dir,
    species=species,
    gene_sets=None  # Use defaults, or specify custom list
)

print(f"\nSuccessfully downloaded {len(downloaded_files)} gene set files:")
for file_path in downloaded_files:
    print(f"  - {file_path}")

print("\nGene sets are ready for offline use!")
print(f"Use local_gene_sets='{gene_sets_dir}' in your compute job.")
```

**Using the CLI:**
```bash
# On login node with internet access
cell2fate download-genesets --species Mouse --output-dir /path/to/your/project/gene_sets
```

### Step 2: Run Analysis on Compute Node (offline)

Your compute job script will use the downloaded gene sets:

```python
#!/usr/bin/env python3
"""
This script runs on compute nodes WITHOUT internet access.
It uses the gene sets downloaded in Step 1.
"""

import cell2fate as c2f
import scanpy as sc

# Load your data
adata = sc.read_h5ad("your_data.h5ad")

# Setup and train model (same as before)
c2f.Cell2fate_DynamicalModel.setup_anndata(
    adata, 
    spliced_label='spliced', 
    unspliced_label='unspliced'
)

n_modules = c2f.utils.get_max_modules(adata)
mod = c2f.Cell2fate_DynamicalModel(adata, n_modules=n_modules)
mod.train()

# Export posteriors
adata = mod.export_posterior(adata)
adata = mod.compute_module_summary_statistics(adata)

# Run enrichment with offline gene sets
gene_sets_dir = "/path/to/your/project/gene_sets"  # Same path as Step 1
background = list(adata.var_names)

tab, results = mod.get_module_top_features(
    adata=adata,
    background=background,
    species='Mouse',  # or 'Human'
    p_adj_cutoff=0.01,
    n_top_genes=100,
    local_gene_sets=gene_sets_dir,  # Use offline gene sets
    gene_sets=None  # Use defaults, or specify custom list
)

# Save results
tab.to_csv("module_top_features.csv")
print("Analysis complete!")
```

## Example Job Submission Script

Here's a complete example for Wynton (or similar HPC):

```bash
#!/bin/bash
#$ -cwd
#$ -o job_output.log
#$ -e job_error.log
#$ -l mem_free=150G
#$ -l scratch=150G
#$ -l h_rt=1209600
#$ -pe smp 8

# Your environment setup
source activate cell2fate_env  # or micromamba activate

# Run the analysis (uses offline gene sets)
python my_cell2fate_analysis.py
```

## Customizing Gene Sets

You can specify which gene sets to use with the `gene_sets` parameter:

```python
# Default gene sets (if gene_sets=None):
# - Mouse: ['GO_Biological_Process_2021']
# - Human: ['GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'KEGG_2021_Human']

# Custom gene sets for Mouse
tab, results = mod.get_module_top_features(
    adata=adata,
    background=background,
    species='Mouse',
    local_gene_sets=gene_sets_dir,
    gene_sets=['GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'KEGG_2019_Mouse']
)

# Custom gene sets for Human
tab, results = mod.get_module_top_features(
    adata=adata,
    background=background,
    species='Human',
    local_gene_sets=gene_sets_dir,
    gene_sets=['KEGG_2021_Human']  # Only KEGG pathways
)
```

## Troubleshooting

### Error: No gene set files found

If you see:
```
Error: No gene set files found in /path/to/gene_sets.
Looked for: ['GO_Biological_Process_2021.gmt']
```

**Solution:**
1. Make sure you ran Step 1 to download the gene sets
2. Check that the path is correct and accessible from compute nodes
3. Verify the files exist: `ls /path/to/gene_sets/*.gmt`

### Error: Connection refused / Name resolution error

If you see connection errors on compute nodes:
```
Failed to resolve 'maayanlab.cloud'
```

**Solution:**
This means you're trying to use online enrichment on a node without internet.
Make sure you:
1. Downloaded gene sets beforehand (Step 1)
2. Are passing `local_gene_sets=gene_sets_dir` to `get_module_top_features()`

### File Permissions

Make sure your gene sets directory has proper permissions:
```bash
chmod -R 755 /path/to/gene_sets
```

## Complete Working Example

See `offline_enrichment_example.py` for a complete working example that includes:
- Gene set download
- Model training
- Offline enrichment analysis
- Custom gene set usage

## Available Gene Sets

Common gene sets you can download:

**For Mouse:**
- `GO_Biological_Process_2021`
- `GO_Cellular_Component_2021`
- `GO_Molecular_Function_2021`
- `KEGG_2019_Mouse`

**For Human:**
- `GO_Biological_Process_2021`
- `GO_Cellular_Component_2021`
- `GO_Molecular_Function_2021`
- `KEGG_2021_Human`
- `WikiPathways_2021_Human`
- `Reactome_2022`

## Summary

The key steps are:

1. **On login node (with internet):** Download gene sets using `download_gene_sets()` or the CLI
2. **On compute node (no internet):** Use `local_gene_sets` parameter to use offline enrichment
3. The `gene_sets` parameter controls which gene sets are used, both online and offline

This approach ensures your analysis works on compute nodes without internet access while maintaining full functionality.

