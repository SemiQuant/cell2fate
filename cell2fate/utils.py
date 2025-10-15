from typing import Optional
import torch
import numpy as np
import pyro
import pandas as pd
import scanpy as sc
import random
import scvelo as scv
from numpy.linalg import norm
import scipy
from scipy.sparse import csr_matrix
from numpy import inner
import matplotlib.pyplot as plt
from contextlib import contextmanager
import seaborn as sns
import os,sys
import cell2fate as c2f
import click
import requests
from pathlib import Path

import torch

def robust_optimization(mod, save_dir, max_epochs = [200, 400], lr = [0.01, 0.01], use_gpu = True):
    """
    Perform robust optimization of the model.

    Parameters
    ----------
    mod
        Model to optimize.
    save_dir
        The directory to save the optimized model.
    max_epochs
        List of maximum epochs for the first and second optimization runs. Default is [200, 400].
    lr
        List of learning rates for the first and second optimization runs. Default is [0.01, 0.01].

    Returns
    -------
    PyroBaseModuleClass
        The optimized model.
    """
    
    n_modules = mod.module.model.n_modules
    adata = mod.adata
    print('First optimization run.')
    mod.train(use_gpu=use_gpu, max_epochs = max_epochs[0], lr = lr[0])
    sample_kwarg = {"num_samples": 1, "batch_size" : 1000,
                     "use_gpu" : use_gpu, 'return_samples': False}
    mod.adata = mod.export_posterior(mod.adata, sample_kwargs=sample_kwarg)
    t_c = np.argsort(np.array(mod.samples['post_sample_means']['t_c']).flatten())/len(np.array(mod.samples['post_sample_means']['t_c']))
    t_c_reversed = -1*(t_c - np.max(t_c))
    print('Second optimization run.')
    del mod
    mod1 = c2f.Cell2fate_DynamicalModel(adata, n_modules = n_modules, init_vals = {'t_c': torch.tensor(t_c).reshape([len(t_c), 1, 1])})
    mod1.train(use_gpu=use_gpu, max_epochs = max_epochs[1], lr = lr[1])
    history1 = mod1.history
    mod1.save(save_dir+'c2f_model', overwrite=True)
    mod1.adata.write(save_dir+"c2f_model_anndata.h5ad")
    del mod1
    print('Third optimization run.')
    mod2 = c2f.Cell2fate_DynamicalModel(adata, n_modules = n_modules, init_vals = {'t_c': torch.tensor(t_c_reversed).reshape([len(t_c_reversed), 1, 1])})
    del adata
    mod2.train(use_gpu=use_gpu, max_epochs = max_epochs[1], lr = lr[1])
    history2 = mod2.history

    iter_start=0
    iter_end=-1

    fig, ax = plt.subplots(1,2, figsize = (15,5))

    iter_end = len(history1["elbo_train"])

    ax[0].plot(
        history1["elbo_train"].index[iter_start:iter_end],
        np.array(history1["elbo_train"].values.flatten())[iter_start:iter_end],
        label="Original Direction",
    )
    ax[0].plot(
        history2["elbo_train"].index[iter_start:iter_end],
        np.array(history2["elbo_train"].values.flatten())[iter_start:iter_end],
        label="Reversed Direction",
    )
    ax[0].legend()
    ax[0].set_xlim(0, len(history1["elbo_train"]))
    ax[0].set_xlabel("Training epochs")
    ax[0].set_ylabel("-ELBO loss")

    ax[1].plot(
        history1["elbo_train"].index[iter_start:iter_end],
        np.array(history1["elbo_train"].values.flatten())[iter_start:iter_end],
        label="Original Direction",
    )
    ax[1].plot(
        history2["elbo_train"].index[iter_start:iter_end],
        np.array(history2["elbo_train"].values.flatten())[iter_start:iter_end],
        label="Reversed Direction",
    )
    ax[1].legend()
    ax[1].set_xlim(np.round(0.8*len(history1["elbo_train"])), len(history1["elbo_train"]))
    ax[1].set_xlabel("Training epochs")
    ax[1].set_ylabel("-ELBO loss")
    plt.tight_layout()
    plt.show()
    
    print('-ELBO second run:', np.mean(np.array(history1['elbo_train'][-40:])))
    print('-ELBO third run:', np.mean(np.array(history2['elbo_train'][-40:])))
    
    if np.mean(np.array(history1['elbo_train'][-40:])) > np.mean(np.array(history2['elbo_train'][-40:])):
        return mod2
    else:
        del mod2
        adata = sc.read_h5ad(save_dir+"c2f_model_anndata.h5ad")
        mod1 = c2f.Cell2fate_DynamicalModel.load(save_dir+'c2f_model', adata)
        return mod1 

def get_max_modules(adata):
    """
    Lightning module task to train Pyro scvi-tools modules.

    Parameters
    ----------
    adata
        AnnData object containing single-cell RNA sequencing data.

    Returns
    -------
    Int
        The maximal number of modules determined through Leiden clustering.
    """
    print('Leiden clustering ...')
    adata_copy = adata.copy()
    adata_copy.X = adata_copy.layers['unspliced'] + adata_copy.layers['spliced']
    sc.pp.normalize_total(adata_copy, target_sum=1e4)
    sc.pp.log1p(adata_copy)
    sc.pp.highly_variable_genes(adata_copy, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_copy = adata_copy[:, adata_copy.var.highly_variable]
    sc.pp.scale(adata_copy, max_value=10)
    sc.pp.neighbors(adata_copy)
    adata_copy = sc.tl.leiden(adata_copy, resolution = 0.75, copy = True)
    print('Number of Leiden Clusters: ' + str(len(np.unique(adata_copy.obs['leiden']))))
    print('Maximal Number of Modules: ' + str(int(1.15*np.round(len(np.unique(adata_copy.obs['leiden']))))))
    n_modules = int(np.round(1.15*len(np.unique(adata_copy.obs['leiden']))))
    del adata_copy
    return n_modules

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

def multiplot_from_generator(g, num_columns, figsize_for_one_row=None, row_size = 15):
    # Copied from: https://aaronmcdaid.com/blog.posts/multiplot_from_generator/
    # on 28/04/2022
    # Plots multiple plots side by side in a jupyter notebook
        next(g)
        # default to 15-inch rows, with square subplots
        if figsize_for_one_row is None:
            figsize_for_one_row = (row_size, row_size/num_columns)
        try:
            while True:
                # call plt.figure once per row
                plt.figure(figsize=figsize_for_one_row)
                for col in range(num_columns):
                    ax = plt.subplot(1, num_columns, col+1)
                    next(g)
        except StopIteration:
            pass

def compute_velocity_graph_Bergen2020(adata, n_neighbours = None, full_posterior = True, spliced_key = 'Ms'):
    """
    Computes a "velocity graph" similar to the method in:
    "Bergen et al. (2020), Generalizing RNA velocity to transient cell states through dynamical modeling"
    
    Parameters
    ----------
    adata
        anndata object with velocity information in ``adata.layers['velocity']`` (expectation value) or 
        ``adata.uns['velocity_posterior']`` (full posterior). Also normalized spliced counts in adata.layers['spliced_norm'].
    n_neighbours
        how many nearest neighbours to consider (all non nearest neighbours have edge weights set to 0)
        if not specified, 10% of the total number of cells is used.
    full_posterior
        whether to use full posterior to compute velocity graph (otherwise expectation value is used)  
    Returns
    -------
    Array
        Velocity graph
    """
    M = len(adata.obs_names)
    if not n_neighbours:
        n_neighbours = int(np.round(M*0.05, 0))
    scv.pp.neighbors(adata, n_neighbors = n_neighbours)
    adata.obsp['binary'] = adata.obsp['connectivities'] != 0
    distances = []
    velocities = []
    cosines = []
    transition_probabilities = []
    matrices = []
    if full_posterior:
        for i in range(M):
            distances += [adata.layers[spliced_key][adata.obsp['binary'].toarray()[i,:],:] - adata.layers[spliced_key][i,:].flatten()]
            velocities += [adata.uns['velocity_posterior'][:,i,:]]
            cosines += [inner(distances[i], velocities[i])/(norm(distances[i])*norm(velocities[i]))]
            transition_probabilities += [np.exp(2*cosines[i])]
            transition_probabilities[i] = transition_probabilities[i]/np.sum(transition_probabilities[i], axis = 0)
            matrices += [csr_matrix((np.mean(np.array(transition_probabilities[i]), axis = 1),
                    (np.repeat(i, len(transition_probabilities[i])), np.where(adata.obsp['binary'][i,:].toarray())[1])),
                                   shape=(M, M))]
    else:
        for i in range(M):
            distances += [adata.layers[spliced_key][adata.obsp['binary'].toarray()[i,:],:] - adata.layers[spliced_key][i,:].flatten()]
            velocities += [adata.layers['velocity'][i,:].reshape(1,len(adata.var_names))]
            cosines += [inner(distances[i], velocities[i])/(norm(distances[i])*norm(velocities[i]))]
            transition_probabilities += [np.exp(2*cosines[i])]
            transition_probabilities[i] = transition_probabilities[i]/np.sum(transition_probabilities[i], axis = 0)
            matrices += [csr_matrix((np.mean(np.array(transition_probabilities[i]), axis = 1),
                    (np.repeat(i, len(transition_probabilities[i])), np.where(adata.obsp['binary'][i,:].toarray())[1])),
                                   shape=(M, M))]
    return sum(matrices)

def plot_velocity_umap_Bergen2020(adata, use_full_posterior = True, n_neighbours = None,
                                  plotting_kwargs = None, save = False, spliced_key = 'Ms'):
    """
    Visualizes RNAvelocity with arrows on a UMAP, using the method introduced in 
    "Bergen et al. (2020), Generalizing RNA velocity to transient cell states through dynamical modeling"
    The method computes a "velocity graph" before plotting (see the referenced paper for details),
    unless such a graph is already available. The graph is based on the full velocity posterior distribution if available,
    otherwise it is based on the velocity expectation values.
    Velocity is expected in ``adata.layers['velocity']`` or ``adata.uns['velocity_posterior']`` (for the full posterior)
    and the graph is saved/expected in ``adata.layers['velocity_graph']``.
    
    Parameters
    ----------
    adata
        AnnData object with velocity information
    use_full_posterior
        Use full posterior to compute velocity graph (if available)
    plotting_kwargs        
    """
    if 'velocity_graph' in adata.uns.keys():
        print('Using existing velocity graph')
        scv.pl.velocity_embedding_stream(adata, basis='umap', **plotting_kwargs)
    else:
        if use_full_posterior and 'velocity_posterior' in adata.uns.keys():
            print('Using full velocity posterior to calculate velocity graph')
            adata.uns['velocity_graph'] = compute_velocity_graph_Bergen2020(adata,
                                                                            full_posterior = True,
                                                                            n_neighbours = n_neighbours,
                                                                            spliced_key = spliced_key)          
        elif use_full_posterior and 'velocity_posterior' not in adata.uns.keys():
            print('Full velocity posterior not found, using expectation value to calculate velocity graph')
            adata.uns['velocity_graph'] = compute_velocity_graph_Bergen2020(adata,
                                                                            full_posterior = False,
                                                                            n_neighbours = n_neighbours,
                                                                            spliced_key = spliced_key)
        elif not use_full_posterior:
            print('Using velocity expectation value to calculate velocity graph.')
            adata.uns['velocity_graph'] = compute_velocity_graph_Bergen2020(adata,
                                                                            full_posterior = False,
                                                                            n_neighbours = n_neighbours,
                                                                            spliced_key = spliced_key)
        
        scv.pl.velocity_embedding_stream(adata, basis='umap', save = save, **plotting_kwargs)

def get_training_data(adata, remove_clusters = None, cells_per_cluster = 100,
                         cluster_column = 'clusters', min_shared_counts = 10, n_var_genes = 2000):
    """
    Reduces and anndata object to the most relevant cells and genes for understanding the differentiation trajectories
    in the data.
    
    Parameters
    ----------
    adata
        AnnData Object
    remove_clusters
        Names of clusters to be removed
    cells_per_cluster
        How many cells to keep per cluster. For Louvain clustering with resolution = 1, keeping more than 300 cells
        per cluster does not provide much extra information.
    cluster_column
        Name of the column in ``adata.obs`` that contains cluster names
    min_shared_counts
        Minimum number of spliced+unspliced counts across all cells for a gene to be retained
    n_var_genes
        Number of top variable genes to retain
        
    Returns
    -------
    AnnData
        AnnData object reduced to the most informative cells and genes
    """
    random.seed(a=1)
    adata = adata[[c not in remove_clusters for c in adata.obs[cluster_column]], :]
    # Restrict samples per cell type:
    N = cells_per_cluster
    unique_celltypes = np.unique(adata.obs[cluster_column])
    index = []
    for i in range(len(unique_celltypes)):
        if adata.obs[cluster_column].value_counts()[unique_celltypes[i]] > N:
            subset = np.where(adata.obs[cluster_column] == unique_celltypes[i])[0]
            subset = random.sample(list(subset), N)
        else:
            subset = np.where(adata.obs[cluster_column] == unique_celltypes[i])[0]
        index += list(subset)
    adata = adata[index,:]
    print('Keeping at most ' + str(N) + ' cells per cluster')
    scv.pp.filter_genes(adata, min_shared_counts=min_shared_counts)
    sc.pp.normalize_total(adata, target_sum=1e4)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=n_var_genes)
    if scipy.sparse.issparse(adata.layers['spliced']):
        adata.layers['spliced'] = np.array(adata.layers['spliced'].toarray(), dtype=np.float32)
    if scipy.sparse.issparse(adata.layers['unspliced']):
        adata.layers['unspliced'] = np.array(adata.layers['unspliced'].toarray(), dtype=np.float32)
    return adata

def G_a(mu, sd):
    """
    Converts mean and standard deviation for a Gamma distribution into the shape parameter.

    Parameters
    ----------
    mu
        The mean of the Gamma distribution.
    sd
        The standard deviation of the Gamma distribution.

    Returns
    -------
    Float
        The shape parameter of the Gamma distribution.
    """
    return mu**2/sd**2

def G_b(mu, sd):
    """
    Converts mean and standard deviation for a Gamma distribution into the scale parameter.

    Parameters
    ----------
    mu
        The mean of the Gamma distribution.
    sd
        The standard deviation of the Gamma distribution.

    Returns
    -------
    Float
        The scale parameter of the Gamma distribution.
    """
    
    return mu/sd**2

def mu_alpha(alpha_new, alpha_old, tau, lam):
    '''
    Calculates transcription rate as a function of new target transcription rate,
    old transcription rate at changepoint, time since changepoint, and rate of exponential change process.

    Parameters
    ----------
    alpha_new
        The new target transcription rate.
    alpha_old
        The old transcription rate at the changepoint.
    tau
        Time since the changepoint.
    lam
        Rate of the exponential change process.

    Returns
    -------
    Float
        The calculated transcription rate.
    '''

    return (alpha_new - alpha_old) * (1 - torch.exp(-lam*tau)) + alpha_old

def mu_mRNA_continuousAlpha(alpha, beta, gamma, tau, u0, s0, delta_alpha, lam):
    '''
    Calculates the expected value of spliced and unspliced counts as a function of rates, latent time, initial states,
    difference to transcription rate in the previous state, and rate of exponential change process between states.

    Parameters
    ----------
    alpha
        Transcription rate.
    beta
        Splicing rate.
    gamma
        Degradation rate.
    tau
        Latent time.
    u0
        Initial unspliced count.
    s0
        Initial spliced count.
    delta_alpha
        Difference to transcription rate in the previous state.
    lam
        Rate of exponential change process between states.

    Returns
    -------
    Torch.Tensor
        A tensor containing the expected value of unspliced and spliced counts.
    '''

    mu_u = u0*torch.exp(-beta*tau) + (alpha/beta)* (1 - torch.exp(-beta*tau)) + delta_alpha/(beta-lam+10**(-5))*(torch.exp(-beta*tau) - torch.exp(-lam*tau))
    mu_s = (s0*torch.exp(-gamma*tau) + 
    alpha/gamma * (1 - torch.exp(-gamma*tau)) +
    (alpha - beta * u0)/(gamma - beta+10**(-5)) * (torch.exp(-gamma*tau) - torch.exp(-beta*tau)) +
    (delta_alpha*beta)/((beta - lam+10**(-5))*(gamma - beta+10**(-5))) * (torch.exp(-beta*tau) - torch.exp(-gamma*tau))-
    (delta_alpha*beta)/((beta - lam+10**(-5))*(gamma - lam+10**(-5))) * (torch.exp(-lam*tau) - torch.exp(-gamma*tau)))

    return torch.stack([mu_u, mu_s], axis = -1)

def mu_mRNA_continousAlpha_globalTime_twoStates(alpha_ON, alpha_OFF, beta, gamma, lam_gi, T_c, T_gON, T_gOFF, Zeros):
    '''
    Calculates the expected value of spliced and unspliced counts as a function of rates,
    global latent time, initial states, and global switch times between two states.

    Parameters
    ----------
    alpha_ON
        Transcription rate when the gene is turned ON.
    alpha_OFF
        Transcription rate when the gene is turned OFF.
    beta
        Splicing rate.
    gamma
        Degradation rate.
    lam_gi
        Parameters for the rate of exponential change process between states.
    T_c
        Global latent time.
    T_gON
        Global switch time when the gene is turned ON.
    T_gOFF
        Global switch time when the gene is turned OFF.
    Zeros
        Tensor with zeros used for initialization.

    Returns
    -------
    Torch.Tensor
        The expected value of unspliced and spliced counts.
    '''
    
    n_cells = T_c.shape[-2]
    n_genes = alpha_ON.shape[-1]
    tau = torch.clip(T_c - T_gON, min = 10**(-5))
    t0 = T_gOFF - T_gON
    # Transcription rate in each cell for each gene:
    boolean = (tau < t0).reshape(n_cells, 1)
    alpha_cg = alpha_ON*boolean + alpha_OFF*~boolean
    # Time since changepoint for each cell and gene:
    tau_cg = tau*boolean + (tau - t0)*~boolean
    # Initial condition for each cell and gene:
    lam_g = ~boolean*lam_gi[:,1] + boolean*lam_gi[:,0]
    initial_state = mu_mRNA_continuousAlpha(alpha_ON, beta, gamma, t0,
                                                       Zeros, Zeros, alpha_ON-alpha_OFF, lam_gi[:,0])
    initial_alpha = mu_alpha(alpha_ON, alpha_OFF, t0, lam_gi[:,0])
    u0_g = 10**(-5) + ~boolean*initial_state[:,:,0]
    s0_g = 10**(-5) + ~boolean*initial_state[:,:,1]
    delta_alpha = ~boolean*initial_alpha*(-1) + boolean*alpha_ON*(1)
    alpha_0 = alpha_OFF + ~boolean*initial_alpha
    # Unspliced and spliced count variance for each gene in each cell:
    # TODO: For PyTorch 2.0+ compatibility, replace torch.clip with torch.clamp
    # torch.clip was deprecated and removed in PyTorch 2.0
    # CHANGE: mu_RNAvelocity = torch.clip(mu_mRNA_continuousAlpha(alpha_cg, beta, gamma, tau_cg,
    #                                                              u0_g, s0_g, delta_alpha, lam_g), min = 10**(-5))
    # TO:     mu_RNAvelocity = torch.clamp(mu_mRNA_continuousAlpha(alpha_cg, beta, gamma, tau_cg,
    #                                                              u0_g, s0_g, delta_alpha, lam_g), min = 10**(-5))
    mu_RNAvelocity = torch.clip(mu_mRNA_continuousAlpha(alpha_cg, beta, gamma, tau_cg,
                                                         u0_g, s0_g, delta_alpha, lam_g), min = 10**(-5))
    return mu_RNAvelocity

def download_gene_sets(output_dir='gene_sets', species='Human', gene_sets=None):
    """
    Download gene sets from Enrichr for offline use.
    
    Parameters
    ----------
    output_dir
        Directory to save gene set files. Defaults to 'gene_sets'.
    species
        Species for gene sets ('Human' or 'Mouse'). Defaults to 'Human'.
    gene_sets
        List of gene set names to download. If None, downloads default sets.
        
    Returns
    -------
    List
        List of downloaded file paths.
    """
    if gene_sets is None:
        if species == 'Human':
            gene_sets = ['GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'KEGG_2021_Human']
        elif species == 'Mouse':
            # Use HDSigDB_Mouse_2021 as default for Mouse gene sets
            gene_sets = ['HDSigDB_Mouse_2021']
        else:
            raise ValueError(f"Unsupported species: {species}. Use 'Human' or 'Mouse'.")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    downloaded_files = []
    
    with click.progressbar(gene_sets, label='Downloading gene sets') as bar:
        for gene_set in bar:
            url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={gene_set}&species={species.lower()}"
            output_file = output_path / f"{gene_set}.gmt"
            
            try:
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                
                with open(output_file, 'w') as f:
                    f.write(response.text)
                
                downloaded_files.append(str(output_file))
                click.echo(f"Downloaded: {gene_set}")
                
            except requests.RequestException as e:
                click.echo(f"Error downloading {gene_set}: {e}", err=True)
                continue
    
    click.echo(f"Downloaded {len(downloaded_files)} gene sets to {output_dir}")
    return downloaded_files

def convert_human_to_mouse_gene_symbols(gene_symbols):
    """
    Convert Human gene symbols to Mouse gene symbols.
    
    Human gene symbols are typically in ALL CAPS (e.g., 'NQO1', 'CEBPB')
    Mouse gene symbols are typically in Title Case (e.g., 'Nqo1', 'Cebpb')
    
    Parameters
    ----------
    gene_symbols : list
        List of Human gene symbols
        
    Returns
    -------
    list
        List of Mouse gene symbols
    """
    mouse_symbols = []
    for symbol in gene_symbols:
        if symbol.isupper() and len(symbol) > 1:
            # Convert ALL CAPS to Title Case
            mouse_symbol = symbol.capitalize()
            mouse_symbols.append(mouse_symbol)
        else:
            # Keep as is if not ALL CAPS
            mouse_symbols.append(symbol)
    return mouse_symbols

def parse_gmt_file(gmt_file, convert_to_mouse=False):
    """
    Parse a GMT (Gene Matrix Transpose) file and return gene sets.
    
    Parameters
    ----------
    gmt_file
        Path to GMT file.
    convert_to_mouse : bool
        If True, convert Human gene symbols to Mouse gene symbols.
        
    Returns
    -------
    Dict
        Dictionary mapping gene set names to lists of genes.
    """
    gene_sets = {}
    
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                gene_set_name = parts[0]
                description = parts[1]
                genes = parts[2:]
                
                if convert_to_mouse:
                    genes = convert_human_to_mouse_gene_symbols(genes)
                
                gene_sets[gene_set_name] = genes
    
    return gene_sets

def offline_enrichment(gene_list, background, gene_sets_dict, p_adj_cutoff=0.01):
    """
    Perform offline enrichment analysis using local gene sets.
    
    Parameters
    ----------
    gene_list
        List of genes to test for enrichment.
    background
        List of background genes.
    gene_sets_dict
        Dictionary mapping gene set names to lists of genes.
    p_adj_cutoff
        Adjusted p-value cutoff for significance.
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with enrichment results.
    """
    import pandas as pd
    from scipy.stats import hypergeom
    from statsmodels.stats.multitest import multipletests
    
    results = []
    
    for gene_set_name, gene_set_genes in gene_sets_dict.items():
        # Convert to sets for efficient operations
        gene_list_set = set(gene_list)
        background_set = set(background)
        gene_set_set = set(gene_set_genes)
        
        # Find overlaps
        genes_in_set = gene_list_set.intersection(gene_set_set)
        genes_in_background = background_set.intersection(gene_set_set)
        
        # Skip if no genes in background
        if len(genes_in_background) == 0:
            continue
            
        # Calculate hypergeometric test
        N = len(background_set)  # Total genes in background
        K = len(genes_in_background)  # Genes in gene set that are in background
        n = len(gene_list)  # Genes in query list
        x = len(genes_in_set)  # Genes in query list that are in gene set
        
        if K == 0 or n == 0:
            continue
            
        # Hypergeometric test
        p_value = hypergeom.sf(x-1, N, K, n)
        
        # Calculate odds ratio
        if x == 0 or (K-x) == 0:
            odds_ratio = 0
        else:
            odds_ratio = (x * (N - K - n + x)) / ((n - x) * (K - x))
        
        results.append({
            'Term': gene_set_name,
            'Overlap': f"{x}/{len(gene_list)}",
            'P-value': p_value,
            'Adjusted P-value': p_value,  # Will be adjusted later
            'Odds Ratio': odds_ratio,
            'Combined Score': x * np.log10(p_value) if p_value > 0 else 0,
            'Genes': ', '.join(sorted(genes_in_set))
        })
    
    if not results:
        return pd.DataFrame()
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Adjust p-values
    if len(df) > 0:
        _, p_adjusted, _, _ = multipletests(df['P-value'], method='fdr_bh')
        df['Adjusted P-value'] = p_adjusted
    
    # Sort by adjusted p-value
    df = df.sort_values('Adjusted P-value')
    
    return df
