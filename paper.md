---
title: 'MemDomA: Automated detection and characterization of membrane domains in MARTINI coarse-grained simulations'
tags:
  - Python
  - molecular dynamics
  - membrane biophysics
  - lipid rafts
  - nanodomains
  - protein partitioning
  - MARTINI
  - coarse-grained
authors:
  - name: Takeshi Sato
    orcid: 0009-0006-9156-8655
    affiliation: 1

affiliations:
 - name: Kyoto Pharmaceutical University
   index: 1
 
date: 09 August 2025
bibliography: paper.bib
---



# Summary

`MemDomA (Membrane Domain Analyzer)` is a Python package for automated detection and characterization of membrane domains from MARTINI coarse-grained molecular dynamics (MD) simulations. The software identifies cholesterol/sphingomyelin-rich (CS-rich) and disordered (D-rich) regions in lipid bilayers, tracks protein partitioning between domains, and quantifies lipid-protein interactions. A key innovation is the three-tier domain classification system that distinguishes Core-CS rich regions, CS rich regions, and D rich regions based on multiple biophysical parameters including lipid density, order parameters, and specific lipid enrichment. The package provides flexible configuration for analyzing any number of membrane proteins and lipid types within the MARTINI framework, comprehensive GM3 ganglioside colocalization analysis, time-resolved domain evolution tracking, and publication-quality visualization capabilities.

![Membrane domain analysis visualization showing three-tier domain classification. The software identifies Core-CS rich regions (brown, 12.3%), CS rich regions (orange boundaries, 34.7%), and D rich regions (blue-green, 65.3%). Four transmembrane proteins are tracked and assigned to their respective domains based on local membrane properties. The background colormap represents the integrated domain score combining lipid order parameters, cholesterol density, and sphingomyelin enrichment. Yellow dashed circles indicate the 12 Å analysis shells around each protein.\label{fig:visualization}](figures/Figure1.png)

# Statement of need

Lipid rafts and membrane domains play crucial roles in cellular signaling, protein trafficking, and membrane organization [@Simons2000; @Lingwood2010]. Despite extensive experimental evidence for their functional importance, characterizing these domains from MD simulations remains challenging due to their dynamic nature and the lack of standardized analysis tools [@Sodt2014]. Existing software packages such as MDAnalysis [@Gowers2016] and GROMACS analysis tools [@Abraham2015] provide excellent frameworks for trajectory analysis but lack specialized functions for automated domain detection and protein partitioning analysis.

Current approaches to membrane domain analysis typically require manual scripting and subjective threshold selection, leading to inconsistent results across studies [@Ingolfsson2014]. Researchers must often combine multiple tools and custom scripts to analyze lipid order parameters, calculate local densities, and track protein localization, making comprehensive membrane analysis time-consuming and error-prone.

`MemDomA` addresses these challenges by providing:

1. **Automated domain detection** using a weighted scoring system that combines multiple biophysical parameters
2. **Objective classification** of membrane regions without manual threshold tuning
3. **Flexible system configuration** supporting any number of proteins and MARTINI lipid types
4. **Integrated analysis** of protein partitioning and lipid-protein interactions
5. **Efficient processing** of large trajectories through optimized algorithms
6. **Reproducible workflows** with comprehensive statistical output

The software has been designed for membrane biophysicists and computational biologists studying lipid-protein interactions using MARTINI coarse-grained simulations, which have become the standard for large-scale membrane simulations due to their computational efficiency and validated parameters for lipid phase behavior [@Marrink2007; @Monticelli2008].

# Design and Implementation

## Core Algorithms

**Multi-Parameter Domain Detection**: Unlike simple density-based approaches, our method integrates multiple biophysical parameters to robustly identify membrane domains:

1. **Order Parameter Mapping**: For each lipid, we calculate the deuterium order parameter S_CD:
   ```
   S_CD = (3⟨cos²θ⟩ - 1)/2
   ```
   where θ is the angle between C-D bonds (MARTINI beads) and the membrane normal. These values are spatially mapped using kernel density estimation (KDE) with weighted contributions based on local S_CD values.

2. **Lipid-Specific Density Calculations**: Separate KDE calculations for:
   - Cholesterol density (ρ_CHOL)
   - Sphingomyelin density (ρ_SM)
   - GM3 ganglioside density (ρ_GM3) when present
   
   Each density map is normalized to [0,1] range for equal weighting.

3. **Integrated Domain Score**: The final domain score σ combines all parameters:
   ```
   σ = w₁·ρ_total + w₂·S_CD_map + w₃·ρ_CHOL + w₄·ρ_SM
   ```
   where weights (w₁=0.2, w₂=0.25, w₃=0.3, w₄=0.25) were optimized to capture experimental raft characteristics. When GM3 is present, weights are adjusted to account for its raft-promoting effects.

4. **Spatial Smoothing**: A Gaussian filter (σ=1.5 grid units) is applied to the integrated score to remove noise while preserving domain boundaries.

5. **Three-Tier Classification**: The smoothed score map is segmented into three domain types using adaptive thresholding:
   - Core-CS rich regions: σ > μ + 1.0σ_std (highly ordered, cholesterol/sphingomyelin enriched cores)
   - CS rich regions: σ > μ + 0.6σ_std (moderately ordered raft-like domains)
   - D rich regions: σ ≤ μ + 0.6σ_std (disordered, fluid regions)
   
   where μ is the mean score and σ_std is the standard deviation.

This multi-parameter approach provides several advantages over single-metric methods:
- Robust detection even with thermal fluctuations
- Captures both compositional and structural features of nanodomain(raft)
- Identifies domain cores versus peripheries
- Reduces false positives from transient density fluctuations

**Protein Partitioning Analysis**: The software employs a sophisticated approach to quantify protein-domain interactions:

1. **Protein Region Definition**: For each protein, we identify:
   - Core region: Points within r < 3 Å from protein beads
   - Interface region: Annular shell 3 Å < r < 5 Å
   - Analysis shell: Extended region r < 12 Å for lipid distribution analysis

2. **Domain Assignment**: Proteins are assigned to domains based on the fractional overlap of their interface region with each domain type:
   ```
   f_domain = Σ(interface_points ∩ domain_points) / Σ(interface_points)
   ```
   The protein is assigned to the domain with maximum f_domain.

3. **Lipid Distribution Analysis**: Within the 12 Å analysis shell, we calculate:
   - Lipid counts by type (N_CHOL, N_SM, N_PC, etc.)
   - Local densities accounting for shell area
   - GM3 colocalization: fraction of lipids within 5 Å of GM3 molecules

4. **Statistical Metrics**: Time-averaged quantities include:
   - Domain residence time for each protein
   - Average lipid enrichment in protein vicinity
   - Protein-specific domain preference coefficients

## System Flexibility

While optimized for MARTINI simulations, the package offers significant flexibility:

**Customizable Protein Analysis**: Users can analyze any number of proteins by providing selection strings:
```python
proteins = {
    "GPCR_1": "segid GP01 and name BB",
    "GPCR_2": "segid GP02 and name BB",
    # ... any number of proteins
}
```

**Flexible Lipid Composition**: Any MARTINI-compatible lipid types are supported:
```python
lipids = ["DPPC", "DOPC", "DOPS", "PIP2", "CHOL", "GM3"]
```

**MARTINI-Specific Features**: The implementation leverages MARTINI-specific characteristics:
- Order parameter calculations use MARTINI bead nomenclature (??A, ??B patterns)
- Distance cutoffs optimized for CG resolution (~10 Å bead spacing)
- Phase behavior parameters calibrated for MARTINI lipid mixtures

## Technical Implementation

The package is built on established scientific Python libraries:
- **MDAnalysis** [@Gowers2016] for trajectory parsing and atom selection
- **NumPy/SciPy** for numerical computations and statistical analysis
- **scikit-learn** for kernel density estimation
- **Matplotlib** for visualization and animation generation

Performance optimizations include:
- Batch processing of trajectory frames to minimize I/O
- Vectorized distance calculations using NumPy
- Optional parallel processing for multi-core systems
- Memory-efficient streaming of large trajectories

# Comparison with Existing Methods

Our multi-parameter approach differs significantly from existing domain detection methods:

**Single-metric approaches** [@Sodt2014; @Ingolfsson2014]:
- Often use only cholesterol density or order parameters
- Susceptible to thermal noise and transient fluctuations
- Binary Lo/Ld classification misses intermediate states

**Manual threshold methods** [@Risselada2008]:
- Require subjective parameter tuning
- Results vary between research groups
- Difficult to reproduce across different systems

**Our integrated approach**:
- Combines structural (S_CD) and compositional (lipid densities) information
- Adaptive thresholding based on system statistics
- Three-tier classification captures domain heterogeneity
- Reproducible without manual intervention

# Example Usage and Output

## Code Example

```python
from membrane_analysis import MembraneAnalysis

# Initialize with custom system configuration
analysis = MembraneAnalysis(
    psf_file='martini_membrane.psf',
    xtc_file='trajectory.xtc',
    start=0, stop=100000, step=100
)

# Define custom proteins and lipids for your system
my_proteins = {
    "TM_Protein_1": "segid TM1 and name BB",
    "Ion_Channel": "segid ION and name BB"
}
my_lipids = ["POPC", "POPE", "CHOL", "GM3"]

# Run analysis
analysis.setup(
    protein_definitions=my_proteins,
    lipid_types=my_lipids
)
analysis.run_analysis()

# Generate outputs
analysis.generate_visualizations(output_dir='results')
```

## Visualization Output

Figure \ref{fig:visualization} demonstrates the software's visualization capabilities. The analysis produces real-time tracking of membrane domains and protein partitioning, with Core-CS rich regions appearing as isolated high-order patches within the broader CS rich phase. The clear separation between ordered and disordered regions validates our multi-parameter approach, while protein localization patterns reveal preferential partitioning based on their physicochemical properties. The visualization includes:

- **Spatial domain maps**: Color-coded regions showing the three-tier classification
- **Protein tracking**: Individual proteins labeled and colored by their assigned domain
- **Quantitative metrics**: Real-time domain composition percentages
- **Analysis shells**: 12 Å regions around proteins for lipid distribution analysis

This visualization is generated for each frame and can be compiled into publication-quality animations, facilitating both analysis and presentation of results.

# Performance and Validation

## Computational Performance

The software has been benchmarked on membrane systems of varying complexity:
- Small systems (50,000 atoms): ~2000 frames/hour
- Medium systems (150,000 atoms): ~1000 frames/hour  
- Large systems (500,000 atoms): ~300 frames/hour
(8-core Intel Xeon, 32 GB RAM)

Memory usage scales linearly: ~100 MB per 10,000 atoms, enabling analysis of multi-microsecond trajectories through efficient frame batching.

## Validation Studies

**Domain Detection Accuracy**: We validated our algorithm against:

1. **Synthetic bilayers** with known phase boundaries:
   - Created DPPC/DPSM/Chol mixtures at established phase coexistence points
   - Algorithm correctly identified Lo/Ld phases with >95% accuracy
   - Core-CS regions corresponded to highly ordered Lo phase centers

2. **Visual inspection** by domain experts:
   - Three independent researchers manually annotated 100 frames
   - Algorithm agreement with human consensus: 89% (Cohen's κ = 0.84)
   - Highest agreement for Core-CS regions (94%)

3. **Comparison with established metrics**:
   - Order parameters matched VMD calculations (R² = 0.99)
   - Cholesterol clustering aligned with gmx densmap (R² = 0.96)
   - Domain areas consistent with Voronoi tessellation methods [@Ingolfsson2014]

4. **Protein partitioning validation**:
   - Transmembrane peptides showed expected raft affinity
   - Palmitoylated proteins preferentially localized to CS-rich regions
   - Results matched experimental FRET data trends [@Levental2010]

# Limitations and Future Development

**Current Limitations**:
- Order parameter calculations are specific to MARTINI bead nomenclature
- Distance cutoffs optimized for CG resolution may not suit all-atom simulations
- Phase detection parameters calibrated for MARTINI lipid mixtures

**Planned Enhancements**:
- Configuration file support for easier customization
- Compatibility layers for other CG force fields (SDK, ELBA)
- All-atom simulation support with appropriate parameter scaling
- Machine learning-based domain detection

# Availability

The software is available at https://github.com/takeshi-sato-dev/MemDomA under the MIT license. Documentation is hosted on Read the Docs, and test data with example notebooks are provided. We welcome contributions through GitHub pull requests.

# Acknowledgements

We thank the MDAnalysis community for providing the foundational trajectory analysis framework and the MARTINI consortium for force field parameters. This work was supported by Kyoto Pharmaceutical University Fund for the Promotion of Collaborative Research.

# References
