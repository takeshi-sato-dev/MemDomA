"""
Static plot generation for membrane analysis.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import List, Dict, Optional, Tuple
import MDAnalysis as mda

from ..utils.constants import FIGURE_SIZE


class PlotGenerator:
    """Generates static plots for membrane analysis results."""
    
    def plot_domain_area_time_series(self,
                                    all_frames_data: List[Tuple],
                                    start_frame: int,
                                    step: int,
                                    universe: mda.Universe,
                                    output_prefix: str):
        """
        Plot domain area fraction time series.
        
        Parameters
        ----------
        all_frames_data : list
            List of frame data tuples
        start_frame : int
            Starting frame number
        step : int
            Frame step size
        universe : mda.Universe
            MDAnalysis universe for time information
        output_prefix : str
            Prefix for output files
        """
        # Collect time and area fraction data
        times = []
        core_cs_area_fractions = []
        cs_area_fractions = []
        d_area_fractions = []
        
        for i, (_, _, _, _, domain_info) in enumerate(all_frames_data):
            if domain_info is not None and 'domain_stats' in domain_info:
                stats = domain_info['domain_stats']
                frame_number = start_frame + i * step
                
                # Get simulation time
                try:
                    universe.trajectory[frame_number]
                    time_ns = universe.trajectory.time / 1000  # ps to ns
                except Exception:
                    time_ns = frame_number * 0.02  # Fallback
                
                times.append(time_ns)
                
                # Collect area fractions
                if 'area_fraction_core_cs' in stats:
                    core_cs_area_fractions.append(stats['area_fraction_core_cs'] * 100)
                else:
                    core_cs_area_fractions.append(0)
                
                cs_area_fractions.append(stats['area_fraction_cs'] * 100)
                d_area_fractions.append(stats['area_fraction_d'] * 100)
        
        # Create main plot
        plt.figure(figsize=(12, 6))
        
        plt.plot(times, core_cs_area_fractions, 'r-', linewidth=2, label='Core-CS rich')
        plt.plot(times, cs_area_fractions, 'orange', linewidth=2, label='CS rich Region')
        plt.plot(times, d_area_fractions, 'c-', linewidth=2, label='D rich Region')
        
        # Add moving average
        window_size = min(25, len(times) // 10) if len(times) > 50 else 5
        if window_size > 1:
            valid_times = times[window_size-1:]
            plt.plot(valid_times, self._moving_average(core_cs_area_fractions, window_size),
                    'r--', linewidth=1.5, alpha=0.7, label='Core-CS (Moving Avg)')
            plt.plot(valid_times, self._moving_average(cs_area_fractions, window_size),
                    'orange', linestyle='--', linewidth=1.5, alpha=0.7, label='CS rich (Moving Avg)')
            plt.plot(valid_times, self._moving_average(d_area_fractions, window_size),
                    'c--', linewidth=1.5, alpha=0.7, label='D rich (Moving Avg)')
        
        plt.xlabel('Simulation Time (ns)', fontsize=14)
        plt.ylabel('Domain Area Fraction (%)', fontsize=14)
        plt.title('Membrane Domain Area Fraction vs. Time', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.legend(loc='center right')
        plt.ylim(0, 100)
        plt.tight_layout()
        
        # Save plots
        plt.savefig(f'{output_prefix}.png', dpi=300)
        plt.savefig(f'{output_prefix}.pdf')
        plt.close()
        
        # Create separate plots
        self._plot_separate_domains(times, core_cs_area_fractions, 
                                   cs_area_fractions, d_area_fractions,
                                   window_size, output_prefix)
        
        # Save data to CSV
        self._save_time_series_data(times, core_cs_area_fractions,
                                   cs_area_fractions, d_area_fractions,
                                   output_prefix)
        
        print(f"Domain area time series plots saved as '{output_prefix}.png/.pdf'")
    
    def plot_protein_domain_distribution(self,
                                        protein_stats: Dict,
                                        start_frame: int,
                                        step: int,
                                        output_prefix: str):
        """
        Plot protein domain distribution bar chart.
        
        Parameters
        ----------
        protein_stats : dict
            Protein statistics dictionary
        start_frame : int
            Starting frame number
        step : int
            Frame step size
        output_prefix : str
            Prefix for output files
        """
        plt.figure(figsize=(12, 8))
        
        # Extract data for each protein
        proteins = list(protein_stats.keys())
        x = np.arange(len(proteins))
        width = 0.25
        
        # Calculate percentages
        core_cs_percentages = []
        cs_percentages = []
        d_percentages = []
        
        for protein_name, domains in protein_stats.items():
            total = domains['total']
            
            core_cs_frames = domains['core_cs']['frames'] if 'core_cs' in domains else 0
            cs_frames = domains['cs']['frames'] if 'cs' in domains else 0
            d_frames = domains['d']['frames'] if 'd' in domains else 0
            
            core_cs_pct = (core_cs_frames / total * 100) if total > 0 else 0
            cs_pct = (cs_frames / total * 100) if total > 0 else 0
            d_pct = (d_frames / total * 100) if total > 0 else 0
            
            core_cs_percentages.append(core_cs_pct)
            cs_percentages.append(cs_pct)
            d_percentages.append(d_pct)
        
        # Create bar chart
        plt.bar(x - width, core_cs_percentages, width, label='Core-CS rich', color='red', alpha=0.7)
        plt.bar(x, cs_percentages, width, label='CS rich Region', color='orange', alpha=0.7)
        plt.bar(x + width, d_percentages, width, label='D rich Region', color='cyan', alpha=0.7)
        
        plt.xlabel('Protein', fontsize=14)
        plt.ylabel('Time in Domain (%)', fontsize=14)
        plt.title('Protein Domain Distribution', fontsize=16)
        plt.xticks(x, proteins, fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend()
        plt.grid(axis='y', alpha=0.3)
        
        # Add value labels
        for i, v in enumerate(core_cs_percentages):
            plt.text(i - width, v + 1, f"{v:.1f}%", ha='center', fontsize=10)
        for i, v in enumerate(cs_percentages):
            plt.text(i, v + 1, f"{v:.1f}%", ha='center', fontsize=10)
        for i, v in enumerate(d_percentages):
            plt.text(i + width, v + 1, f"{v:.1f}%", ha='center', fontsize=10)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}.png', dpi=300)
        plt.savefig(f'{output_prefix}.pdf')
        plt.close()
        
        print(f"Protein domain distribution plot saved as '{output_prefix}.png/.pdf'")
    
    def _moving_average(self, data: List[float], window_size: int) -> np.ndarray:
        """Calculate moving average."""
        return np.convolve(data, np.ones(window_size)/window_size, mode='valid')
    
    def _plot_separate_domains(self,
                              times: List[float],
                              core_cs: List[float],
                              cs: List[float],
                              d: List[float],
                              window_size: int,
                              output_prefix: str):
        """Create separate plots for each domain."""
        plt.figure(figsize=(12, 10))
        
        valid_times = times[window_size-1:] if window_size > 1 else times
        
        # Core-CS rich
        plt.subplot(3, 1, 1)
        plt.plot(times, core_cs, 'r-', linewidth=2)
        if window_size > 1:
            plt.plot(valid_times, self._moving_average(core_cs, window_size),
                    'k--', linewidth=1.5)
        plt.title('Core-CS rich Area Fraction vs. Time', fontsize=14)
        plt.ylabel('Area Fraction (%)', fontsize=12)
        plt.grid(True, alpha=0.3)
        
        # CS rich
        plt.subplot(3, 1, 2)
        plt.plot(times, cs, color='orange', linewidth=2)
        if window_size > 1:
            plt.plot(valid_times, self._moving_average(cs, window_size),
                    'k--', linewidth=1.5)
        plt.title('CS rich Region Area Fraction vs. Time', fontsize=14)
        plt.ylabel('Area Fraction (%)', fontsize=12)
        plt.grid(True, alpha=0.3)
        
        # D rich
        plt.subplot(3, 1, 3)
        plt.plot(times, d, 'c-', linewidth=2)
        if window_size > 1:
            plt.plot(valid_times, self._moving_average(d, window_size),
                    'k--', linewidth=1.5)
        plt.title('D rich Region Area Fraction vs. Time', fontsize=14)
        plt.xlabel('Simulation Time (ns)', fontsize=12)
        plt.ylabel('Area Fraction (%)', fontsize=12)
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix.replace("time_series", "separate_plots")}.png', dpi=300)
        plt.savefig(f'{output_prefix.replace("time_series", "separate_plots")}.pdf')
        plt.close()
    
    def _save_time_series_data(self,
                              times: List[float],
                              core_cs: List[float],
                              cs: List[float],
                              d: List[float],
                              output_prefix: str):
        """Save time series data to CSV."""
        df = pd.DataFrame({
            'Time (ns)': times,
            'Core-CS rich (%)': core_cs,
            'CS rich Region (%)': cs,
            'D rich Region (%)': d
        })
        csv_file = f'{output_prefix}.csv'
        df.to_csv(csv_file, index=False)
        print(f"Time series data saved as '{csv_file}'")