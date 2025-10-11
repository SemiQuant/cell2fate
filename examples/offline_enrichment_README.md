# Cell2fate Offline Enrichment

Use cell2fate enrichment analysis without internet access.

> **🖥️ Running on HPC?** See the comprehensive guide: [hpc_offline_enrichment_guide.md](hpc_offline_enrichment_guide.md)

## Quick Start

**Step 1: Download gene sets** (needs internet)
```bash
python examples/download_gene_sets_for_hpc.py --species Mouse --output-dir gene_sets
```

**Step 2: Use in your analysis** (works offline)
```python
import cell2fate as c2f

# ... your model setup and training ...

# Add local_gene_sets parameter for offline enrichment
tab, results = mod.get_module_top_features(
    adata=adata,
    background=list(adata.var_names),
    species='Mouse',
    local_gene_sets='gene_sets'  # ← Use offline gene sets
)
```

## Files in This Directory

- **`download_gene_sets_for_hpc.py`** - Script to download gene sets for offline use
- **`hpc_offline_enrichment_guide.md`** - Comprehensive guide with troubleshooting
- **`offline_enrichment_example.py`** - Complete working example code

## Custom Gene Sets

Customize which gene sets to use:

```python
# Download specific gene sets
python examples/download_gene_sets_for_hpc.py \
    --species Mouse \
    --gene-sets GO_Biological_Process_2021 KEGG_2019_Mouse

# Use them in analysis
tab, results = mod.get_module_top_features(
    adata=adata,
    background=background,
    species='Mouse',
    local_gene_sets='gene_sets',
    gene_sets=['GO_Biological_Process_2021', 'KEGG_2019_Mouse']
)
```

**Default gene sets:**
- Mouse: `GO_Biological_Process_2021`
- Human: `GO_Biological_Process_2021`, `GO_Cellular_Component_2021`, `KEGG_2021_Human`

**For more details, see:** [hpc_offline_enrichment_guide.md](hpc_offline_enrichment_guide.md)
