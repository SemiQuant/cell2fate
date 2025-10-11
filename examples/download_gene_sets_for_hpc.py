#!/usr/bin/env python3
"""
Download Gene Sets for Offline Enrichment Analysis

Run this script on a login node with internet access BEFORE submitting
your compute job. This downloads gene sets from Enrichr that will be
used offline on compute nodes without internet access.

Usage:
    python download_gene_sets_for_hpc.py --species Mouse --output-dir gene_sets
    
Or using cell2fate CLI:
    cell2fate download-genesets --species Mouse --output-dir gene_sets

Author: SemiQuant (JasonLimberis@ucsf.edu)
"""

import argparse
import os
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description='Download gene sets for offline enrichment analysis in HPC environments',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download Mouse gene sets
  python download_gene_sets_for_hpc.py --species Mouse
  
  # Download Human gene sets to custom directory
  python download_gene_sets_for_hpc.py --species Human --output-dir /wynton/home/mylab/user/gene_sets
  
  # Download specific gene sets
  python download_gene_sets_for_hpc.py --species Mouse --gene-sets GO_Biological_Process_2021 KEGG_2019_Mouse

After downloading, use in your analysis:
  tab, results = mod.get_module_top_features(
      adata=adata,
      background=background,
      species='Mouse',
      local_gene_sets='gene_sets',  # Path to downloaded gene sets
      gene_sets=None  # Or specify custom list
  )
        """
    )
    
    parser.add_argument(
        '--species',
        type=str,
        required=True,
        choices=['Mouse', 'Human'],
        help='Species for gene sets (Mouse or Human)'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default='gene_sets',
        help='Directory to save gene sets (default: gene_sets)'
    )
    
    parser.add_argument(
        '--gene-sets',
        nargs='*',
        default=None,
        help='Specific gene sets to download (default: use species defaults)'
    )
    
    parser.add_argument(
        '--force',
        action='store_true',
        help='Re-download even if files already exist'
    )
    
    args = parser.parse_args()
    
    # Import cell2fate
    try:
        import cell2fate as c2f
    except ImportError:
        print("Error: cell2fate package not found.")
        print("Please install it first: pip install cell2fate")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("Cell2fate Gene Set Downloader for HPC")
    print("=" * 70)
    print(f"\nSpecies: {args.species}")
    print(f"Output directory: {output_dir.absolute()}")
    
    # Show default gene sets
    if args.gene_sets is None:
        if args.species == 'Mouse':
            default_sets = ['GO_Biological_Process_2021']
        else:  # Human
            default_sets = ['GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'KEGG_2021_Human']
        print(f"Gene sets (default): {default_sets}")
    else:
        print(f"Gene sets (custom): {args.gene_sets}")
    
    # Check if files already exist
    if not args.force:
        gene_sets_to_check = args.gene_sets if args.gene_sets else default_sets
        existing_files = []
        for gs in gene_sets_to_check:
            gmt_file = output_dir / f"{gs}.gmt"
            if gmt_file.exists():
                existing_files.append(gmt_file.name)
        
        if existing_files:
            print(f"\nWarning: {len(existing_files)} gene set file(s) already exist:")
            for f in existing_files:
                print(f"  - {f}")
            response = input("Re-download? (y/N): ")
            if response.lower() != 'y':
                print("Skipping download. Use --force to re-download automatically.")
                sys.exit(0)
    
    # Download gene sets
    print("\nDownloading gene sets...")
    print("This may take a few minutes depending on your internet connection.")
    
    try:
        downloaded_files = c2f.utils.download_gene_sets(
            output_dir=str(output_dir),
            species=args.species,
            gene_sets=args.gene_sets
        )
        
        print("\n" + "=" * 70)
        print("SUCCESS!")
        print("=" * 70)
        print(f"\nDownloaded {len(downloaded_files)} gene set file(s):")
        for file_path in downloaded_files:
            file_size = Path(file_path).stat().st_size / 1024  # KB
            print(f"  ✓ {Path(file_path).name} ({file_size:.1f} KB)")
        
        print(f"\nGene sets saved to: {output_dir.absolute()}")
        print("\nNext steps:")
        print("1. Submit your compute job")
        print("2. In your Python script, use:")
        print(f"   local_gene_sets='{output_dir.absolute()}'")
        print("\nExample code:")
        print("  tab, results = mod.get_module_top_features(")
        print("      adata=adata,")
        print("      background=background,")
        print(f"      species='{args.species}',")
        print(f"      local_gene_sets='{output_dir.absolute()}',")
        print("      gene_sets=None  # Or specify custom list")
        print("  )")
        print("\nFor more information, see examples/hpc_offline_enrichment_guide.md")
        
    except Exception as e:
        print("\n" + "=" * 70)
        print("ERROR!")
        print("=" * 70)
        print(f"\nFailed to download gene sets: {e}")
        print("\nPossible causes:")
        print("  - No internet connection")
        print("  - Enrichr service is down")
        print("  - Network firewall blocking access")
        print("\nPlease check your internet connection and try again.")
        sys.exit(1)

if __name__ == "__main__":
    main()

