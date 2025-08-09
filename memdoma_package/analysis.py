"""
Main analysis class for membrane domain analysis - Fixed version.
"""

import os
import numpy as np
from typing import Dict, List, Optional, Tuple
from tqdm import tqdm

from .core.trajectory import TrajectoryHandler
from .core.lipids import LipidAnalyzer
from .core.proteins import ProteinAnalyzer
from .core.domains import DomainAnalyzer
from .processing.parallel import ParallelProcessor
from .visualization.plots import PlotGenerator
from .visualization.animation import AnimationGenerator
from .utils.io import save_statistics
from .utils.constants import (
    DEFAULT_START, DEFAULT_STOP, DEFAULT_STEP,
    BATCH_SIZE, STANDARD_LIPIDS
)


class MemDomA:
    """Main class for membrane domain analysis."""
    
    def __init__(self,
                 psf_file: str,
                 xtc_file: str,
                 start: int = DEFAULT_START,
                 stop: int = DEFAULT_STOP,
                 step: int = DEFAULT_STEP,
                 use_multiprocessing: bool = False):  # 追加: デフォルトでFalse
        """Initialize membrane analysis."""
        self.psf_file = psf_file
        self.xtc_file = xtc_file
        self.start = start
        self.stop = stop
        self.step = step
        self.use_multiprocessing = use_multiprocessing  # 追加
        
        # Initialize handlers
        self.trajectory_handler = TrajectoryHandler(psf_file, xtc_file)
        self.lipid_analyzer = LipidAnalyzer()
        self.protein_analyzer = ProteinAnalyzer()
        self.domain_analyzer = DomainAnalyzer()
        self.parallel_processor = ParallelProcessor()
        self.plot_generator = PlotGenerator()
        self.animation_generator = AnimationGenerator()
        
        # Data storage
        self.universe = None
        self.leaflet = None
        self.selections = None
        self.proteins = None
        self.protein_stats = None
        self.all_frames_data = None
        
    def setup(self,
              lipid_types: Optional[List[str]] = None,
              protein_definitions: Optional[Dict[str, str]] = None):
        """Set up the analysis by loading trajectory and making selections."""
        # Load trajectory
        self.universe = self.trajectory_handler.load_universe()
        
        # Identify leaflets
        self.leaflet = self.trajectory_handler.identify_lipid_leaflets()
        
        # Set default lipid types if not provided
        if lipid_types is None:
            lipid_types = STANDARD_LIPIDS
        
        # Select lipids
        self.selections = self.trajectory_handler.select_lipids_and_cholesterol(
            lipid_types, self.leaflet
        )
        
        # Set default protein definitions if not provided
        if protein_definitions is None:
            protein_definitions = {
                "Protein 1": "segid PROA and resid 65:76",
                "Protein 2": "segid PROB and resid 65:76",
                "Protein 3": "segid PROC and resid 65:76",
                "Protein 4": "segid PROD and resid 65:76"
            }
        
        # Select proteins
        self.proteins = self.trajectory_handler.select_proteins(protein_definitions)
        
        # Initialize protein statistics
        self.protein_stats = self.protein_analyzer.initialize_protein_statistics(
            list(self.proteins.keys())
        )
        
        print("Analysis setup complete.")
    
    def run_analysis(self):
        """Run the main analysis pipeline."""
        if self.universe is None:
            raise RuntimeError("Analysis not set up. Call setup() first.")
        
        # Create frame batches
        num_frames = int((self.stop - self.start) / self.step)
        frame_batches = self._create_frame_batches(num_frames)
        
        print(f"Total frames to process: {num_frames}")
        print(f"Processing in {len(frame_batches)} batches")
        
        # シングルプロセスモードで処理
        print("Using single process mode (safer for MDAnalysis)...")
        self.all_frames_data = self._process_trajectory_single(
            frame_batches,
            self.universe,
            self.selections,
            self.proteins,
            self.protein_stats
        )
        
        print("Analysis complete.")
    
    def _process_trajectory_single(self, frame_batches, u, selections, proteins, protein_stats):
        """Process trajectory in single process mode."""
        all_frames_data = []
        
        num_batches = len(frame_batches)
        with tqdm(total=num_batches,
                 desc="Processing batches (single process)",
                 unit="batch") as pbar:
            
            for batch in frame_batches:
                batch_results = self.parallel_processor.process_frame_batch(
                    batch, u, selections, proteins
                )
                
                for frame_data in batch_results:
                    coordinates_df, protein_regions, lipid_stats, density_map, domain_info = frame_data
                    
                    # Update statistics
                    self.protein_analyzer.update_protein_statistics(
                        protein_stats, coordinates_df, protein_regions,
                        domain_info, lipid_stats, proteins
                    )
                    
                    all_frames_data.append(frame_data)
                
                pbar.update(1)
        
        return all_frames_data
    
    def generate_visualizations(self,
                              output_dir: str = ".",
                              animation: bool = True,
                              plots: bool = True):
        """Generate all visualizations."""
        if self.all_frames_data is None:
            raise RuntimeError("Analysis not run. Call run_analysis() first.")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Save statistics
        stats_file = os.path.join(output_dir, 'domain_statistics.csv')
        save_statistics(self.protein_stats, stats_file)
        
        if plots:
            print("Generating plots...")
            
            # Domain area time series
            plot_file = os.path.join(output_dir, 'domain_area_time_series')
            self.plot_generator.plot_domain_area_time_series(
                self.all_frames_data,
                self.start,
                self.step,
                self.universe,
                plot_file
            )
            
            # Protein domain distribution
            plot_file = os.path.join(output_dir, 'protein_domain_distribution')
            self.plot_generator.plot_protein_domain_distribution(
                self.protein_stats,
                self.start,
                self.step,
                plot_file
            )
        
        if animation:
            print("Generating animation...")
            animation_file = os.path.join(output_dir, 'membrane_domain_analysis.mp4')
            self.animation_generator.create_animation(
                self.all_frames_data,
                self.universe,
                self.start,
                self.step,
                animation_file
            )
        
        print(f"Visualizations saved to {output_dir}")
    
    def _create_frame_batches(self, num_frames: int) -> List[List[int]]:
        """Create batches of frame numbers for parallel processing."""
        frame_batches = []
        for i in range(0, num_frames, BATCH_SIZE):
            batch = [
                self.start + frame * self.step
                for frame in range(i, min(i + BATCH_SIZE, num_frames))
            ]
            frame_batches.append(batch)
        return frame_batches
    
    def get_results(self) -> Dict:
        """Get analysis results."""
        return {
            'protein_stats': self.protein_stats,
            'frames_data': self.all_frames_data,
            'parameters': {
                'start': self.start,
                'stop': self.stop,
                'step': self.step
            }
        }
