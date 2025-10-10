#!/usr/bin/env python3
"""
Command-line interface for cell2fate utilities.
"""

import click
from cell2fate.utils import download_gene_sets


@click.group()
def cli():
    """Cell2fate command-line utilities."""
    pass


@cli.command()
@click.option('--output-dir', '-o', default='gene_sets', 
              help='Directory to save gene set files (default: gene_sets)')
@click.option('--species', '-s', type=click.Choice(['Human', 'Mouse']), default='Human',
              help='Species for gene sets (default: Human)')
@click.option('--gene-sets', '-g', multiple=True,
              help='Specific gene sets to download (can be used multiple times)')
def download_genesets(output_dir, species, gene_sets):
    """Download gene sets from Enrichr for offline use.
    
    This command downloads gene set files (.gmt) that can be used for offline
    enrichment analysis when internet access is not available.
    
    Examples:
    
        # Download default gene sets for Human
        cell2fate download-genesets --species Human
        
        # Download default gene sets for Mouse
        cell2fate download-genesets --species Mouse
        
        # Download specific gene sets
        cell2fate download-genesets --gene-sets GO_Biological_Process_2021 --gene-sets KEGG_2021_Human
        
        # Download to custom directory
        cell2fate download-genesets --output-dir /path/to/gene_sets
    """
    gene_sets_list = list(gene_sets) if gene_sets else None
    
    try:
        downloaded_files = download_gene_sets(
            output_dir=output_dir,
            species=species,
            gene_sets=gene_sets_list
        )
        
        if downloaded_files:
            click.echo(f"\nSuccessfully downloaded {len(downloaded_files)} gene set files:")
            for file_path in downloaded_files:
                click.echo(f"  - {file_path}")
            click.echo(f"\nTo use these gene sets offline, pass the directory path to get_module_top_features:")
            click.echo(f"  model.get_module_top_features(adata, background, local_gene_sets='{output_dir}')")
        else:
            click.echo("No gene sets were downloaded.", err=True)
            
    except Exception as e:
        click.echo(f"Error downloading gene sets: {e}", err=True)
        raise click.Abort()


if __name__ == '__main__':
    cli()
