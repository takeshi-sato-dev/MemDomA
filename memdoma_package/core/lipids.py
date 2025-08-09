"""
Lipid analysis functions including order parameter calculations.
"""

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.spatial import distance
from typing import Dict, Optional, Tuple, Any
import MDAnalysis as mda

from ..utils.constants import (
    CHOL_SHELL_RADIUS, CHAIN_A_PATTERN, CHAIN_B_PATTERN,
    GRID_SPACING
)


class LipidAnalyzer:
    """Handles lipid-specific analysis operations."""
    
    @staticmethod
    def calculate_order_parameter(residue,  # mda.Residueを削除
                                 chainA: str = CHAIN_A_PATTERN,
                                 chainB: str = CHAIN_B_PATTERN) -> float:
        """
        Calculate order parameter for a lipid residue.
        
        Parameters
        ----------
        residue : MDAnalysis.Residue
            Lipid residue to analyze
        chainA : str
            Selection pattern for chain A atoms
        chainB : str
            Selection pattern for chain B atoms
            
        Returns
        -------
        float
            Order parameter S_CD, or NaN if calculation fails
        """
        chain_atomsA = residue.atoms.select_atoms(chainA)
        chain_atomsB = residue.atoms.select_atoms(chainB)
        
        if len(chain_atomsA) < 2 or len(chain_atomsB) < 2:
            return np.nan

        def get_chain_vector(chain_atoms):
            vectors = np.diff(chain_atoms.positions, axis=0)
            norms = np.linalg.norm(vectors, axis=1)
            norms[norms == 0] = np.finfo(float).eps
            unit_vectors = vectors / norms[:, np.newaxis]
            return unit_vectors[:, 2]
        
        cos_thetasA = get_chain_vector(chain_atomsA)
        cos_thetasB = get_chain_vector(chain_atomsB)

        S_CD_A = (3 * np.nanmean(cos_thetasA ** 2) - 1) / 2
        S_CD_B = (3 * np.nanmean(cos_thetasB ** 2) - 1) / 2

        S_CD = (S_CD_A + S_CD_B) / 2
        return S_CD if np.isfinite(S_CD) else np.nan
    
    @staticmethod
    def analyze_lipid_distribution(protein_com: np.ndarray,
                                  selections: Dict[str, Any],  # mda.AtomGroupも同様
                                  box_dimensions: np.ndarray) -> Dict[str, Dict[str, float]]:
        """
        Analyze lipid distribution around a protein.
        
        Parameters
        ----------
        protein_com : np.ndarray
            Center of mass of the protein (x, y)
        selections : dict
            Dictionary of lipid selections
        box_dimensions : np.ndarray
            Box dimensions for PBC
            
        Returns
        -------
        dict
            Lipid statistics including counts, densities, and colocalization
        """
        shell_area = np.pi * (CHOL_SHELL_RADIUS**2) / 100
        
        # Initialize statistics for all lipid types
        lipid_stats = {
            lipid_type: {
                'count': 0,
                'density': 0.0,
                'gm3_colocalization': 0.0
            }
            for lipid_type in ['CHOL', 'DIPC', 'DPSM', 'DPG3']
        }
        
        # Calculate GM3 positions first
        gm3_positions = None
        if 'DPG3' in selections:
            gm3_positions = np.array([
                res.atoms.center_of_mass()[:2]
                for res in selections['DPG3'].residues
            ])
        
        # Analyze each lipid type
        for lipid_type, selection in selections.items():
            if lipid_type in lipid_stats:
                lipid_positions = np.array([
                    res.atoms.center_of_mass()[:2]
                    for res in selection.residues
                ])
                
                if len(lipid_positions) > 0:
                    # Apply periodic boundary conditions
                    dx = lipid_positions[:, 0] - protein_com[0]
                    dy = lipid_positions[:, 1] - protein_com[1]
                    
                    dx -= box_dimensions[0] * np.round(dx / box_dimensions[0])
                    dy -= box_dimensions[1] * np.round(dy / box_dimensions[1])
                    
                    distances = np.sqrt(dx**2 + dy**2)
                    in_shell = distances <= CHOL_SHELL_RADIUS
                    lipid_count = np.sum(in_shell)
                    
                    stats = lipid_stats[lipid_type]
                    stats['count'] = lipid_count
                    stats['density'] = lipid_count / shell_area if shell_area > 0 else 0.0
                    
                    # GM3 colocalization analysis
                    if gm3_positions is not None and len(gm3_positions) > 0 and lipid_type != 'DPG3':
                        lipid_pos_in_shell = lipid_positions[in_shell]
                        if len(lipid_pos_in_shell) > 0:
                            gm3_distances = distance.cdist(lipid_pos_in_shell, gm3_positions)
                            colocalization = np.mean(np.any(gm3_distances <= 5.0, axis=1))
                            stats['gm3_colocalization'] = colocalization
        
        return lipid_stats
    
    @staticmethod
    def calculate_op_kde(lipid_data: pd.DataFrame,
                        exclude_points: Optional[np.ndarray] = None,
                        bw_method: float = 0.15) -> Optional[gaussian_kde]:
        """
        Calculate KDE for order parameters.
        
        Parameters
        ----------
        lipid_data : pd.DataFrame
            DataFrame with x, y, and S_CD columns
        exclude_points : np.ndarray, optional
            Points to exclude from KDE
        bw_method : float
            Bandwidth method for KDE
            
        Returns
        -------
        gaussian_kde or None
            KDE object or None if insufficient data
        """
        positions = lipid_data[['x', 'y']].values
        op_values = lipid_data['S_CD'].values
        
        # Filter valid data
        valid_mask = np.isfinite(op_values)
        positions = positions[valid_mask]
        op_values = op_values[valid_mask]
        
        if exclude_points is not None:
            mask = np.ones(len(positions), dtype=bool)
            distances = distance.cdist(positions, exclude_points)
            mask &= np.all(distances > GRID_SPACING, axis=1)
            positions = positions[mask]
            op_values = op_values[mask]
        
        if len(positions) < 2:
            return None
        
        # Check for finite values
        if not np.all(np.isfinite(positions)):
            finite_mask = np.all(np.isfinite(positions), axis=1)
            positions = positions[finite_mask]
            op_values = op_values[finite_mask]
            
            if len(positions) < 2:
                return None
        
        try:
            weights = op_values - np.min(op_values) + 1e-10
            if not np.all(np.isfinite(weights)):
                return None
                
            return gaussian_kde(positions.T, weights=weights, bw_method=bw_method)
        except Exception as e:
            print(f"Error in calculate_op_kde: {str(e)}")
            return None
    
    @staticmethod
    def extract_lipid_positions(selections: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """
        Extract 2D positions for all lipid types.
        
        Parameters
        ----------
        selections : dict
            Dictionary of lipid selections
            
        Returns
        -------
        dict
            Dictionary mapping lipid types to position arrays
        """
        selection_positions = {}
        
        for lipid_type, selection in selections.items():
            positions = []
            for res in selection.residues:
                com = res.atoms.center_of_mass()
                if not isinstance(com, np.ndarray):
                    com = np.array(com)
                if len(com.shape) == 1:
                    com = com.reshape(1, -1)
                positions.append(com[:, :2])
            
            if positions:
                selection_positions[lipid_type] = np.vstack(positions)
            else:
                selection_positions[lipid_type] = np.array([]).reshape(0, 2)
        
        return selection_positions
