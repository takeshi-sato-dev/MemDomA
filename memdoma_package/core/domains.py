"""
Domain identification and analysis functions.
"""

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter
from typing import Dict, Optional, Tuple

from ..utils.constants import (
    SMOOTHING_SIGMA, CORE_THRESHOLD_STD, CS_THRESHOLD_STD
)


class DomainAnalyzer:
    """Handles membrane domain identification and analysis."""
    
    @staticmethod
    def calculate_domain_info(density: np.ndarray,
                            lipid_data: pd.DataFrame,
                            lipid_positions: Dict[str, np.ndarray],
                            dimensions: np.ndarray,
                            x_grid: np.ndarray,
                            y_grid: np.ndarray) -> Dict:
        """Calculate membrane domain characteristics."""
        # Prepare grid points
        positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
        
        # Calculate cholesterol density
        if 'CHOL' in lipid_positions and len(lipid_positions['CHOL']) > 0:
            kde = gaussian_kde(lipid_positions['CHOL'].T)
            chol_density = kde(positions).reshape(x_grid.shape)
        else:
            chol_density = np.zeros_like(density)
        
        # Calculate sphingomyelin density
        if 'DPSM' in lipid_positions and len(lipid_positions['DPSM']) > 0:
            kde = gaussian_kde(lipid_positions['DPSM'].T)
            sm_density = kde(positions).reshape(x_grid.shape)
        else:
            sm_density = np.zeros_like(density)
        
        # Calculate GM3 (DPG3) density
        has_gm3 = False
        if 'DPG3' in lipid_positions and len(lipid_positions['DPG3']) > 0:
            kde = gaussian_kde(lipid_positions['DPG3'].T)
            gm3_density = kde(positions).reshape(x_grid.shape)
            has_gm3 = True
        else:
            gm3_density = np.zeros_like(density)
        
        # Grid order parameters
        order_params = DomainAnalyzer._grid_order_parameters(
            lipid_data, dimensions, density.shape
        )
        
        # Normalize maps
        density_norm = DomainAnalyzer._normalize_array(density)
        order_norm = DomainAnalyzer._normalize_array(order_params)
        chol_norm = DomainAnalyzer._normalize_array(chol_density)
        sm_norm = DomainAnalyzer._normalize_array(sm_density)
        
        # Set weights based on GM3 presence
        if has_gm3:
            weights = {'density': 0.2, 'order': 0.25, 'chol': 0.3, 'sm': 0.25}
        else:
            weights = {'density': 0.2, 'order': 0.25, 'chol': 0.3, 'sm': 0.25}
        
        # Calculate integrated score
        cs_rich_score = (
            weights['density'] * density_norm + 
            weights['order'] * order_norm + 
            weights['chol'] * chol_norm + 
            weights['sm'] * sm_norm
        )
        
        # Apply spatial smoothing
        cs_rich_score_smooth = gaussian_filter(cs_rich_score, sigma=SMOOTHING_SIGMA)
        
        # Determine thresholds
        mean_score = np.mean(cs_rich_score_smooth)
        std_score = np.std(cs_rich_score_smooth)
        core_threshold = mean_score + CORE_THRESHOLD_STD * std_score
        cs_threshold = mean_score + CS_THRESHOLD_STD * std_score
        
        # Identify domains
        core_cs_rich = cs_rich_score_smooth > core_threshold
        cs_rich = cs_rich_score_smooth > cs_threshold
        d_rich = ~cs_rich
        
        # Calculate domain statistics
        domain_stats = DomainAnalyzer._calculate_domain_statistics(
            core_cs_rich, cs_rich, d_rich,
            order_params, chol_density, sm_density, gm3_density, has_gm3
        )
        
        return {
            'core_cs_rich': core_cs_rich,
            'cs_rich': cs_rich,
            'd_rich': d_rich,
            'domain_stats': domain_stats,
            'parameters': {
                'weights': weights,
                'smoothing_sigma': SMOOTHING_SIGMA,
                'core_threshold': core_threshold,
                'cs_threshold': cs_threshold
            },
            'raw_data': {
                'density_norm': density_norm,
                'order_norm': order_norm,
                'chol_norm': chol_norm,
                'sm_norm': sm_norm,
                'cs_score': cs_rich_score_smooth
            }
        }
    
    @staticmethod
    def _normalize_array(arr: np.ndarray) -> np.ndarray:
        """Normalize array to [0, 1] range."""
        arr_min = np.min(arr)
        arr_max = np.max(arr)
        if arr_max - arr_min > 0:
            return (arr - arr_min) / (arr_max - arr_min)
        return np.zeros_like(arr)
    
    @staticmethod
    def _grid_order_parameters(lipid_data: pd.DataFrame,
                              dimensions: np.ndarray,
                              shape: Tuple[int, int]) -> np.ndarray:
        """Grid order parameters from lipid data."""
        order_params = np.zeros(shape)
        
        for idx, row in lipid_data.iterrows():
            if not np.isnan(row['S_CD']):
                x_idx = int((row['x'] % dimensions[0]) * shape[0] / dimensions[0])
                y_idx = int((row['y'] % dimensions[1]) * shape[1] / dimensions[1])
                if 0 <= x_idx < shape[0] and 0 <= y_idx < shape[1]:
                    order_params[x_idx, y_idx] = row['S_CD']
        
        return order_params
    
    @staticmethod
    def _calculate_domain_statistics(core_cs_rich: np.ndarray,
                                    cs_rich: np.ndarray,
                                    d_rich: np.ndarray,
                                    order_params: np.ndarray,
                                    chol_density: np.ndarray,
                                    sm_density: np.ndarray,
                                    gm3_density: np.ndarray,
                                    has_gm3: bool) -> Dict:
        """Calculate statistics for each domain."""
        domain_stats = {
            'area_fraction_core_cs': np.sum(core_cs_rich) / core_cs_rich.size,
            'area_fraction_cs': np.sum(cs_rich) / cs_rich.size,
            'area_fraction_d': np.sum(d_rich) / d_rich.size,
            
            'mean_order_core_cs': np.mean(order_params[core_cs_rich]) if np.any(core_cs_rich) else 0,
            'mean_order_cs': np.mean(order_params[cs_rich]) if np.any(cs_rich) else 0,
            'mean_order_d': np.mean(order_params[d_rich]) if np.any(d_rich) else 0,
            
            'mean_chol_core_cs': np.mean(chol_density[core_cs_rich]) if np.any(core_cs_rich) else 0,
            'mean_chol_cs': np.mean(chol_density[cs_rich]) if np.any(cs_rich) else 0,
            'mean_chol_d': np.mean(chol_density[d_rich]) if np.any(d_rich) else 0,
            
            'mean_sm_core_cs': np.mean(sm_density[core_cs_rich]) if np.any(core_cs_rich) else 0,
            'mean_sm_cs': np.mean(sm_density[cs_rich]) if np.any(cs_rich) else 0,
            'mean_sm_d': np.mean(sm_density[d_rich]) if np.any(d_rich) else 0,
        }
        
        # Add GM3 statistics if present
        if has_gm3:
            gm3_stats = {
                'mean_gm3_core_cs': np.mean(gm3_density[core_cs_rich]) if np.any(core_cs_rich) else 0,
                'mean_gm3_cs': np.mean(gm3_density[cs_rich]) if np.any(cs_rich) else 0,
                'mean_gm3_d': np.mean(gm3_density[d_rich]) if np.any(d_rich) else 0,
            }
            domain_stats.update(gm3_stats)
        
        return domain_stats
