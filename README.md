[![Run tests](https://github.com/BayraktarLab/cell2fate/actions/workflows/run_tests.yml/badge.svg)](https://github.com/BayraktarLab/cell2fate/actions/workflows/run_tests.yml)  [![codecov](https://codecov.io/gh/BayraktarLab/cell2fate/graph/badge.svg?token=CCJTK20MA7)](https://codecov.io/gh/BayraktarLab/cell2fate)
[![Documentation Status](https://readthedocs.org/projects/cell2fate/badge/?version=latest)](https://cell2fate.readthedocs.io/en/latest/?badge=latest)

![alt text](https://github.com/BayraktarLab/cell2fate/blob/main/cell2fate_diagram.png?raw=true)

## Usage and Tutorials

Please find our documentation and tutorials [here](https://cell2fate.readthedocs.io/en/latest/).

## Publication figures

Results from all datasets in the [cell2fate preprint](https://www.biorxiv.org/content/10.1101/2023.08.03.551650v1.full.pdf) can be reproduced with [these noteobooks](https://github.com/AlexanderAivazidis/cell2fate_notebooks).

## Installation

We suggest using a separate conda environment for installing cell2fate.

Create a conda environment and install the `cell2fate` package.

```bash
conda create -y -n cell2fate_env python=3.9

conda activate cell2fate_env
pip install git+https://github.com/BayraktarLab/cell2fate
```
If you use a newer GPU with CUDA > 10.2 you will need to install a newer pytorch version into the environment like this:

```bash
conda activate cell2fate_env
pip install torch torchvision --extra-index-url https://download.pytorch.org/whl/cu116
```

To use this environment in a jupyter notebook, add a jupyter kernel for this environment:

```bash
conda activate cell2fate_env
pip install ipykernel
python -m ipykernel install --user --name=cell2fate_env --display-name='Environment (cell2fate_env)'
```

If you do not have conda please install Miniconda first:

```bash
cd /path/to/software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# use prefix /path/to/software/miniconda3
```

Before installing cell2fate and it's dependencies, it could be necessary to make sure that you are creating a fully isolated conda environment by telling python to NOT use user site for installing packages, ideally by adding this line to your `~/.bashrc` file , but this would also work during a terminal session:

```bash
export PYTHONNOUSERSITE="someletters"
```

## HPC Usage and Offline Enrichment

When running cell2fate on HPC systems where compute nodes don't have internet access, you'll need to use offline enrichment analysis. This is a two-step process:

### Step 1: Download Gene Sets (on login node)

Before submitting your compute job, download gene sets on a login node that has internet access:

```bash
# Using the convenient download script
python examples/download_gene_sets_for_hpc.py --species Mouse --output-dir gene_sets

# Or using cell2fate CLI
cell2fate download-genesets --species Mouse --output-dir gene_sets
```

### Step 2: Use Offline Enrichment (on compute node)

In your compute job script, specify the path to downloaded gene sets:

```python
import cell2fate as c2f

# Your analysis code...
mod = c2f.Cell2fate_DynamicalModel(adata, n_modules=n_modules)
mod.train()

# Use offline enrichment
tab, results = mod.get_module_top_features(
    adata=adata,
    background=list(adata.var_names),
    species='Mouse',
    local_gene_sets='gene_sets',  # Path to downloaded gene sets
)
```

For detailed instructions and troubleshooting, see:
- `examples/hpc_offline_enrichment_guide.md` - Comprehensive HPC guide
- `examples/offline_enrichment_README.md` - General offline enrichment documentation
- `examples/offline_enrichment_example.py` - Working code example
