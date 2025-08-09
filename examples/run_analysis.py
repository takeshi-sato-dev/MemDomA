#!/usr/bin/env python
"""
Run MemDomA analysis on your data.
This is a template - modify the parameters for your system.
"""

import os
import sys
import argparse
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from memdoma_package import MemDomA


def main():
    """Main function for running membrane domain analysis."""
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='MemDomA - Analyze membrane domains from MD simulation'
    )
    parser.add_argument('psf', help='Path to PSF topology file')
    parser.add_argument('xtc', help='Path to XTC trajectory file')
    parser.add_argument(
        '--start',
        type=int,
        default=0,
        help='Starting frame (default: 0)'
    )
    parser.add_argument(
        '--stop',
        type=int,
        default=1000,
        help='Stopping frame (default: 1000)'
    )
    parser.add_argument(
        '--step',
        type=int,
        default=10,
        help='Frame step size (default: 10)'
    )
    parser.add_argument(
        '--output',
        default='results',
        help='Output directory (default: results)'
    )
    parser.add_argument(
        '--no-animation',
        action='store_true',
        help='Skip animation generation'
    )
    parser.add_argument(
        '--no-plots',
        action='store_true',
        help='Skip plot generation'
    )
    parser.add_argument(
        '--lipids',
        nargs='+',
        default=["DPSM", "DIPC", "DPG3"],
        help='Lipid types to analyze (default: DPSM DIPC DPG3)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 60)
    print("MemDomA - Membrane Domain Analysis")
    print("=" * 60)
    print(f"PSF file: {args.psf}")
    print(f"XTC file: {args.xtc}")
    print(f"Frames: {args.start} to {args.stop} (step {args.step})")
    print(f"Lipids: {', '.join(args.lipids)}")
    print(f"Output directory: {args.output}")
    print("=" * 60)
    
    # Initialize analysis
    analysis = MemDomA(
        psf_file=args.psf,
        xtc_file=args.xtc,
        start=args.start,
        stop=args.stop,
        step=args.step,
        use_multiprocessing=False  # Safer for MDAnalysis
    )
    
    # Set up analysis
    # Default protein definitions - MODIFY FOR YOUR SYSTEM
    protein_definitions = {
        "Protein 1": "segid PROA and resid 65:76",
        "Protein 2": "segid PROB and resid 65:76",
        "Protein 3": "segid PROC and resid 65:76",
        "Protein 4": "segid PROD and resid 65:76"
    }
    
    print("\nSetting up analysis...")
    analysis.setup(
        lipid_types=args.lipids,
        protein_definitions=protein_definitions
    )
    
    # Run analysis
    print("\nRunning analysis...")
    analysis.run_analysis()
    
    # Generate visualizations
    print("\nGenerating outputs...")
    analysis.generate_visualizations(
        output_dir=args.output,
        animation=not args.no_animation,
        plots=not args.no_plots
    )
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print("=" * 60)
    
    # List generated files
    print("\nGenerated files:")
    output_files = [
        'domain_statistics.csv',
        'domain_area_time_series.png',
        'domain_area_time_series.pdf',
        'protein_domain_distribution.png',
        'protein_domain_distribution.pdf'
    ]
    
    if not args.no_animation:
        output_files.append('membrane_domain_analysis.mp4')
    
    for filename in output_files:
        filepath = os.path.join(args.output, filename)
        if os.path.exists(filepath):
            size = os.path.getsize(filepath) / 1024
            print(f"  ✓ {filename} ({size:.1f} KB)")
    
    print(f"\nAnalysis results saved to: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()