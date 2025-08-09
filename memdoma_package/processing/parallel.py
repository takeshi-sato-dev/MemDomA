"""
Parallel processing utilities for frame batch processing.
"""

import numpy as np
import pandas as pd
import multiprocessing as mp
from typing import List, Dict, Tuple, Any
from functools import partial
from tqdm import tqdm
import psutil
from scipy.ndimage import gaussian_filter

from ..core.lipids import LipidAnalyzer
from ..core.proteins import ProteinAnalyzer
from ..core.domains import DomainAnalyzer
from ..utils.constants import DEFAULT_CPU_FRACTION


class ParallelProcessor:
    """Handles parallel processing of trajectory frames."""
    
    def __init__(self, cpu_fraction: float = DEFAULT_CPU_FRACTION):
        """Initialize parallel processor."""
        self.cpu_fraction = cpu_fraction
        self.num_cores = self._calculate_cores()
        
    def _calculate_cores(self) -> int:
        """Calculate number of CPU cores to use."""
        total_cores = mp.cpu_count()
        used_cores = max(1, int(total_cores * self.cpu_fraction))
        print(f"Available CPU cores: {total_cores}")
        print(f"Using {used_cores} CPU cores for processing")
        return used_cores
    
    @staticmethod
    def process_frame_batch(batch_frames: List[int],
                          u: Any,
                          selections: Dict,
                          proteins: Dict) -> List[Tuple]:
        """Process a batch of frames."""
        results = []
        lipid_analyzer = LipidAnalyzer()
        protein_analyzer = ProteinAnalyzer()
        
        for frame_number in batch_frames:
            u.trajectory[frame_number]
            
            # Extract lipid positions
            selection_positions = lipid_analyzer.extract_lipid_positions(selections)
            
            # Initialize coordinate arrays
            coordinates = {
                'type': [], 'resname': [], 'resid': [],
                'x': [], 'y': [], 'S_CD': []
            }
            
            # Process lipids
            for lipid_type, positions in selection_positions.items():
                if len(positions) > 0:
                    n_residues = len(positions)
                    coordinates['type'].extend(['Lipid'] * n_residues)
                    coordinates['resname'].extend([lipid_type] * n_residues)
                    coordinates['resid'].extend(range(n_residues))
                    coordinates['x'].extend(positions[:, 0].flatten())
                    coordinates['y'].extend(positions[:, 1].flatten())
                    
                    if lipid_type != 'CHOL':
                        scd_values = np.array([
                            lipid_analyzer.calculate_order_parameter(res)
                            for res in selections[lipid_type].residues
                        ])
                    else:
                        scd_values = np.full(n_residues, np.nan)
                    coordinates['S_CD'].extend(scd_values.flatten())
            
            # Add proteins
            for protein_name, protein in proteins.items():
                if len(protein) > 0:
                    com = protein.center_of_mass()[:2]
                    coordinates['type'].append(protein_name)
                    coordinates['resname'].append('PROTEIN')
                    coordinates['resid'].append(protein.residues.resids[0])
                    coordinates['x'].append(com[0])
                    coordinates['y'].append(com[1])
                    coordinates['S_CD'].append(np.nan)
            
            df = pd.DataFrame(coordinates)
            
            # Calculate domains
            density_map = None
            domain_info = None
            
            lipid_data = df[df['type'] == 'Lipid']
            if not lipid_data.empty:
                kde = lipid_analyzer.calculate_op_kde(lipid_data)
                if kde is not None:
                    x_grid, y_grid = np.mgrid[0:u.dimensions[0]:100j, 0:u.dimensions[1]:100j]
                    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
                    
                    density_map = gaussian_filter(
                        kde(positions).reshape(x_grid.shape),
                        sigma=2.0
                    )
                    
                    domain_info = DomainAnalyzer.calculate_domain_info(
                        density_map, lipid_data, selection_positions,
                        u.dimensions, x_grid, y_grid
                    )
            
            # Process protein regions and stats
            protein_regions = {}
            lipid_stats = {}
            
            for protein_name, protein in proteins.items():
                if len(protein) > 0:
                    core_points, interface_points = protein_analyzer.calculate_helix_regions(protein)
                    protein_regions[protein_name] = {
                        'core': core_points,
                        'interface': interface_points
                    }
                    
                    protein_com = protein.center_of_mass()[:2]
                    lipid_stats[protein_name] = lipid_analyzer.analyze_lipid_distribution(
                        protein_com, selections, u.dimensions
                    )
            
            results.append((df, protein_regions, lipid_stats, density_map, domain_info))
        
        return results
    
    def process_trajectory(self,
                         frame_batches: List[List[int]],
                         u: Any,
                         selections: Dict,
                         proteins: Dict,
                         protein_stats: Dict) -> List[Tuple]:
        """Process all frame batches in parallel."""
        all_frames_data = []
        protein_analyzer = ProteinAnalyzer()
        
        with mp.Pool(self.num_cores) as pool:
            process_batch_partial = partial(
                self.process_frame_batch,
                u=u,
                selections=selections,
                proteins=proteins
            )
            
            num_batches = len(frame_batches)
            with tqdm(total=num_batches,
                     desc=f"Processing batches using {self.num_cores} cores",
                     unit="batch") as pbar:
                
                for batch_results in pool.imap(process_batch_partial, frame_batches):
                    for frame_data in batch_results:
                        coordinates_df, protein_regions, lipid_stats, density_map, domain_info = frame_data
                        
                        # Update statistics
                        protein_analyzer.update_protein_statistics(
                            protein_stats, coordinates_df, protein_regions,
                            domain_info, lipid_stats, proteins
                        )
                        
                        all_frames_data.append(frame_data)
                    
                    pbar.update(1)
                    process = psutil.Process()
                    pbar.set_postfix({'memory': f'{process.memory_info().rss / 1024 / 1024:.1f}MB'})
        
        return all_frames_data
