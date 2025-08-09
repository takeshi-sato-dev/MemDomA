"""
Protein region analysis functions.
"""

import numpy as np
from scipy.spatial import ConvexHull, distance
from typing import Tuple, Dict
import MDAnalysis as mda

from ..utils.constants import HELIX_RADIUS, INTERFACE_WIDTH, GRID_SPACING


class ProteinAnalyzer:
    """Handles protein-specific analysis operations."""
    
    @staticmethod
    def calculate_helix_regions(protein: mda.AtomGroup) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate core and interface regions around protein helices.
        
        Parameters
        ----------
        protein : mda.AtomGroup
            Protein atoms to analyze
            
        Returns
        -------
        tuple
            Core points and interface points as numpy arrays
        """
        positions = protein.positions[:, :2]
        hull = ConvexHull(positions)
        boundary_points = positions[hull.vertices]
        
        # Calculate grid bounds
        x_min, y_min = np.min(boundary_points, axis=0) - (HELIX_RADIUS + INTERFACE_WIDTH)
        x_max, y_max = np.max(boundary_points, axis=0) + (HELIX_RADIUS + INTERFACE_WIDTH)
        
        # Create grid
        x_grid = np.arange(x_min, x_max, GRID_SPACING)
        y_grid = np.arange(y_min, y_max, GRID_SPACING)
        XX, YY = np.meshgrid(x_grid, y_grid)
        
        grid_points = np.column_stack((XX.ravel(), YY.ravel()))
        
        # Calculate distances to protein
        distances_to_protein = distance.cdist(grid_points, positions[:, :2])
        min_distances = np.min(distances_to_protein, axis=1)
        
        # Identify regions
        core_points = grid_points[min_distances < HELIX_RADIUS]
        interface_points = grid_points[
            (min_distances >= HELIX_RADIUS) & 
            (min_distances < HELIX_RADIUS + INTERFACE_WIDTH)
        ]
        
        return core_points, interface_points
    
    @staticmethod
    def update_protein_statistics(protein_stats: Dict,
                                 coordinates_df: 'pd.DataFrame',
                                 protein_regions: Dict,
                                 domain_info: Dict,
                                 lipid_stats: Dict,
                                 proteins: Dict[str, mda.AtomGroup]) -> None:
        """
        Update protein statistics based on domain location and lipid distribution.
        
        Parameters
        ----------
        protein_stats : dict
            Dictionary to update with statistics
        coordinates_df : pd.DataFrame
            DataFrame with coordinate information
        protein_regions : dict
            Dictionary of protein regions
        domain_info : dict
            Domain information including masks
        lipid_stats : dict
            Lipid distribution statistics
        proteins : dict
            Dictionary of protein AtomGroups
        """
        for protein_name in proteins.keys():
            protein_stats[protein_name]['total'] += 1
            
            if domain_info is not None and protein_name in protein_regions:
                # Create domain mask
                x_grid, y_grid = np.mgrid[0:100:1, 0:100:1]
                surrounding_mask = np.zeros_like(domain_info['cs_rich'], dtype=bool)
                interface_points = protein_regions[protein_name]['interface']
                
                # Calculate indices for interface points
                idx_x = (interface_points[:, 0] * 100 / coordinates_df['x'].max()).astype(int)
                idx_y = (interface_points[:, 1] * 100 / coordinates_df['y'].max()).astype(int)
                valid_indices = (
                    (idx_x >= 0) & (idx_x < 100) & 
                    (idx_y >= 0) & (idx_y < 100)
                )
                
                # Apply valid indices
                valid_x = idx_x[valid_indices]
                valid_y = idx_y[valid_indices]
                if len(valid_x) > 0:
                    surrounding_mask[valid_x, valid_y] = True
                    
                    # Calculate domain occupancy
                    if 'core_cs_rich' in domain_info:
                        core_cs_fraction = np.sum(domain_info['core_cs_rich'] & surrounding_mask) / np.sum(surrounding_mask)
                        cs_fraction = np.sum(domain_info['cs_rich'] & surrounding_mask) / np.sum(surrounding_mask)
                        d_fraction = np.sum(domain_info['d_rich'] & surrounding_mask) / np.sum(surrounding_mask)
                        
                        # Determine dominant domain
                        fractions = [core_cs_fraction, cs_fraction, d_fraction]
                        domains = ['core_cs', 'cs', 'd']
                        domain = domains[np.argmax(fractions)]
                    else:
                        # Fallback to CS/D classification
                        cs_fraction = (np.sum(domain_info['cs_rich'] & surrounding_mask) / 
                                     np.sum(surrounding_mask))
                        domain = 'cs' if cs_fraction > 0.5 else 'd'
                    
                    # Update statistics
                    if domain in protein_stats[protein_name]:
                        stats = protein_stats[protein_name][domain]
                        stats['frames'] += 1
                        
                        # Update lipid statistics
                        if protein_name in lipid_stats:
                            curr_lipid_stats = lipid_stats[protein_name]
                            for lipid_type in ['CHOL', 'DIPC', 'DPSM', 'DPG3']:
                                if lipid_type in curr_lipid_stats:
                                    values = curr_lipid_stats[lipid_type]
                                    stats[f'{lipid_type}_sum'] += values['count']
                                    stats[f'{lipid_type}_density_sum'] += values['density']
                                    
                                    # Update GM3 colocalization if present
                                    if 'gm3_colocalization' in values:
                                        if f'{lipid_type}_gm3_coloc_sum' not in stats:
                                            stats[f'{lipid_type}_gm3_coloc_sum'] = 0.0
                                        stats[f'{lipid_type}_gm3_coloc_sum'] += values['gm3_colocalization']
    
    @staticmethod
    def initialize_protein_statistics(protein_names: list) -> Dict:
        """
        Initialize protein statistics dictionary.
        
        Parameters
        ----------
        protein_names : list
            List of protein names
            
        Returns
        -------
        dict
            Initialized statistics dictionary
        """
        protein_stats = {}
        
        for name in protein_names:
            protein_stats[name] = {
                'total': 0,
                'core_cs': {
                    'frames': 0,
                    'CHOL_sum': 0, 'CHOL_density_sum': 0,
                    'DIPC_sum': 0, 'DIPC_density_sum': 0,
                    'DPSM_sum': 0, 'DPSM_density_sum': 0,
                    'DPG3_sum': 0, 'DPG3_density_sum': 0
                },
                'cs': {
                    'frames': 0,
                    'CHOL_sum': 0, 'CHOL_density_sum': 0,
                    'DIPC_sum': 0, 'DIPC_density_sum': 0,
                    'DPSM_sum': 0, 'DPSM_density_sum': 0,
                    'DPG3_sum': 0, 'DPG3_density_sum': 0
                },
                'd': {
                    'frames': 0,
                    'CHOL_sum': 0, 'CHOL_density_sum': 0,
                    'DIPC_sum': 0, 'DIPC_density_sum': 0,
                    'DPSM_sum': 0, 'DPSM_density_sum': 0,
                    'DPG3_sum': 0, 'DPG3_density_sum': 0
                }
            }
        
        return protein_stats