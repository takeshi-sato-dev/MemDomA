#!/usr/bin/env python
"""
Test MemDomA with included sample data.
This script validates all features of the package.
"""

import os
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from memdoma_package import MemDomA


def test_memdoma():
    """Complete test of MemDomA functionality."""
    
    print("=" * 60)
    print("MemDomA Test Suite")
    print("=" * 60)
    
    # Output directory
    output_dir = 'test_output'
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n1. Initializing MemDomA with test data...")
    analysis = MemDomA(
        psf_file='test_data/test_system.psf',
        xtc_file='test_data/test_trajectory.xtc',
        start=0,
        stop=8,  # Test data has 8 frames
        step=1,
        use_multiprocessing=False  # Single process for testing
    )
    print("   ✓ Initialized")
    
    print("\n2. Setting up analysis...")
    # Configuration for test data
    lipid_types = ["DPSM", "DOPC", "DOPS", "DPG3"]
    protein_definitions = {
        "Protein A": "segid PROA and backbone",
        "Protein B": "segid PROB and backbone"
    }
    
    analysis.setup(
        lipid_types=lipid_types,
        protein_definitions=protein_definitions
    )
    print("   ✓ Setup complete")
    
    print("\n3. Running analysis...")
    analysis.run_analysis()
    print("   ✓ Analysis complete")
    
    print("\n4. Generating visualizations...")
    print("   - Creating statistical plots...")
    print("   - Creating animation...")
    analysis.generate_visualizations(
        output_dir=output_dir,
        animation=True,
        plots=True
    )
    print("   ✓ Visualizations complete")
    
    print("\n" + "=" * 60)
    print("Test Results")
    print("=" * 60)
    
    # Check generated files
    expected_files = [
        'domain_statistics.csv',
        'domain_area_time_series.png',
        'domain_area_time_series.pdf',
        'domain_area_separate_plots.png',
        'domain_area_separate_plots.pdf',
        'protein_domain_distribution.png',
        'protein_domain_distribution.pdf',
        'membrane_domain_analysis.mp4'
    ]
    
    print("\nGenerated files:")
    missing_files = []
    for filename in expected_files:
        filepath = os.path.join(output_dir, filename)
        if os.path.exists(filepath):
            size = os.path.getsize(filepath) / 1024
            print(f"  ✓ {filename:<40} ({size:>8.1f} KB)")
        else:
            print(f"  ✗ {filename:<40} (NOT FOUND)")
            missing_files.append(filename)
    
    # Summary
    print("\n" + "=" * 60)
    if not missing_files:
        print("✓ All tests passed successfully!")
    else:
        print("⚠ Some files were not generated:")
        for f in missing_files:
            print(f"  - {f}")
    
    print(f"\nResults saved in: {os.path.abspath(output_dir)}/")
    print("=" * 60)
    
    return len(missing_files) == 0


if __name__ == "__main__":
    success = test_memdoma()
    sys.exit(0 if success else 1)