"""
Animation generation for membrane analysis.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.ndimage import gaussian_filter
from typing import List, Dict, Tuple
import MDAnalysis as mda

from ..utils.constants import (
    FIGURE_SIZE, DEFAULT_FPS, DEFAULT_DPI, DEFAULT_BITRATE,
    CHOL_SHELL_RADIUS
)


class AnimationGenerator:
    """Generates animations for membrane analysis results."""
    
    def create_animation(self,
                        all_frames_data: List[Tuple],
                        universe: mda.Universe,
                        start_frame: int,
                        step: int,
                        output_file: str):
        """
        Create animation of membrane domain analysis.
        
        Parameters
        ----------
        all_frames_data : list
            List of frame data tuples
        universe : mda.Universe
            MDAnalysis universe
        start_frame : int
            Starting frame number
        step : int
            Frame step size
        output_file : str
            Path to output animation file
        """
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        fig.patch.set_facecolor('black')
        
        num_frames = len(all_frames_data)
        print(f"Creating animation with {num_frames} frames...")
        
        def update_plot(frame):
            """Update plot for each frame."""
            coordinates_df, protein_regions, lipid_stats, density_map, domain_info = all_frames_data[frame]
            ax.clear()
            self._plot_membrane_phase_map(
                ax=ax,
                protein_coords=coordinates_df,
                box_size=universe.dimensions[:2],
                frame_number=start_frame + frame * step,
                density=density_map,
                domain_info=domain_info,
                protein_regions=protein_regions,
                lipid_stats=lipid_stats
            )
        
        # Create animation
        anim = FuncAnimation(
            fig, update_plot, frames=num_frames,
            repeat=False, interval=25, cache_frame_data=False
        )
        
        # Save animation
        print(f"Saving animation to {output_file}...")
        anim.save(
            output_file,
            writer='ffmpeg',
            fps=DEFAULT_FPS,
            dpi=DEFAULT_DPI,
            bitrate=DEFAULT_BITRATE
        )
        
        plt.close(fig)
        print(f"Animation saved successfully as '{output_file}'")
    
    def _plot_membrane_phase_map(self,
                                ax,
                                protein_coords,
                                box_size,
                                frame_number,
                                density,
                                domain_info,
                                protein_regions,
                                lipid_stats):
        """
        Plot a single frame of the membrane phase map.
        
        Parameters
        ----------
        ax : matplotlib axis
            Axis to plot on
        protein_coords : pd.DataFrame
            Protein coordinates
        box_size : tuple
            Box dimensions
        frame_number : int
            Current frame number
        density : np.ndarray
            Density map
        domain_info : dict
            Domain information
        protein_regions : dict
            Protein region data
        lipid_stats : dict
            Lipid statistics
        """
        ax.set_facecolor('black')
        x_min, x_max = 0, box_size[0]
        y_min, y_max = 0, box_size[1]
        x_grid, y_grid = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]
        
        # Plot density map
        if density is not None:
            ax.imshow(density.T,
                     extent=[x_min, x_max, y_min, y_max],
                     cmap='viridis',
                     origin='lower',
                     alpha=0.8)
        
        # Plot domain boundaries
        if domain_info is not None:
            self._plot_domain_boundaries(ax, x_grid, y_grid, domain_info)
        
        # Plot DPG3 lipids
        dpg3_data = protein_coords[protein_coords['resname'] == 'DPG3']
        if not dpg3_data.empty:
            ax.scatter(dpg3_data['x'], dpg3_data['y'],
                      color='white', s=10, alpha=0.5)
        
        # Plot proteins
        for protein_name, regions in protein_regions.items():
            self._plot_protein(ax, protein_name, regions, protein_coords,
                             domain_info, box_size, lipid_stats)
        
        # Set plot properties
        ax.set(title='Membrane Domain Analysis',
               xlabel='X (Å)',
               ylabel='Y (Å)',
               xlim=(x_min, x_max),
               ylim=(y_min, y_max))
        
        ax.title.set_color('white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.tick_params(colors='white')
        
        # Add legend
        self._add_legend(ax, domain_info)
        
        # Add frame information
        self._add_frame_info(ax, frame_number, domain_info)
    
    def _plot_domain_boundaries(self, ax, x_grid, y_grid, domain_info):
        """Plot domain boundaries and fills."""
        if 'core_cs_rich' in domain_info:
            # Core-CS rich regions
            smoothed_core = gaussian_filter(domain_info['core_cs_rich'].astype(float), sigma=1.0)
            ax.contourf(x_grid, y_grid, smoothed_core > 0.5,
                       levels=[0.5, 1.5],
                       colors=['red'],
                       alpha=0.3)
            
            # CS rich regions
            smoothed_cs = gaussian_filter(domain_info['cs_rich'].astype(float), sigma=1.0)
            ax.contourf(x_grid, y_grid, smoothed_cs > 0.5,
                       levels=[0.5, 1.5],
                       colors=['orange'],
                       alpha=0.15)
            
            # Boundaries
            ax.contour(x_grid, y_grid, smoothed_core,
                      levels=[0.5],
                      colors='yellow',
                      linewidths=1.5,
                      linestyles='dashed',
                      alpha=0.8)
            
            ax.contour(x_grid, y_grid, smoothed_cs,
                      levels=[0.5],
                      colors='orange',
                      linewidths=1.5,
                      alpha=0.7)
        else:
            # Fallback to simple CS/D domains
            smoothed_cs_domain = gaussian_filter(
                domain_info['cs_rich'].astype(float), sigma=1.0
            )
            
            ax.contourf(x_grid, y_grid, smoothed_cs_domain > 0.5,
                       levels=[0.5, 1.5],
                       colors=['orange'],
                       alpha=0.2)
            
            ax.contour(x_grid, y_grid, smoothed_cs_domain,
                      levels=[0.5],
                      colors='orange',
                      linewidths=2,
                      alpha=0.8)
    
    def _plot_protein(self, ax, protein_name, regions, protein_coords,
                     domain_info, box_size, lipid_stats):
        """Plot individual protein with its regions."""
        # Plot core and interface
        ax.scatter(regions['core'][:, 0], regions['core'][:, 1],
                  color='white', s=10, alpha=0.3)
        ax.scatter(regions['interface'][:, 0], regions['interface'][:, 1],
                  color='cyan', s=10, alpha=0.3)
        
        protein_data = protein_coords[protein_coords['type'] == protein_name]
        
        if not protein_data.empty and domain_info is not None:
            # Determine protein domain
            protein_color = self._determine_protein_color(
                regions['interface'], domain_info, box_size
            )
            
            # Draw analysis shell
            circle = patches.Circle(
                (protein_data['x'].iloc[0], protein_data['y'].iloc[0]),
                CHOL_SHELL_RADIUS,
                fill=False,
                color='yellow',
                linestyle='--',
                alpha=0.5
            )
            ax.add_patch(circle)
            
            # Draw protein
            ax.scatter(protein_data['x'], protein_data['y'],
                      color=protein_color,
                      s=200,
                      edgecolors='black',
                      linewidths=2)
            
            # Add label
            core_center = np.mean(regions['core'], axis=0)
            ax.text(core_center[0], core_center[1] + 10,
                   protein_name,
                   color='white',
                   fontsize=12,
                   ha='center',
                   va='bottom',
                   fontweight='bold')
    
    def _determine_protein_color(self, interface_points, domain_info, box_size):
        """Determine protein color based on domain location."""
        surrounding_mask = np.zeros_like(domain_info['cs_rich'], dtype=bool)
        
        idx_x = (interface_points[:, 0] * 100 / box_size[0]).astype(int)
        idx_y = (interface_points[:, 1] * 100 / box_size[1]).astype(int)
        valid_indices = (
            (idx_x >= 0) & (idx_x < 100) &
            (idx_y >= 0) & (idx_y < 100)
        )
        
        if np.any(valid_indices):
            surrounding_mask[idx_x[valid_indices], idx_y[valid_indices]] = True
            
            if 'core_cs_rich' in domain_info:
                core_cs_fraction = np.sum(domain_info['core_cs_rich'] & surrounding_mask) / np.sum(surrounding_mask)
                cs_fraction = np.sum(domain_info['cs_rich'] & surrounding_mask) / np.sum(surrounding_mask)
                d_fraction = np.sum(domain_info['d_rich'] & surrounding_mask) / np.sum(surrounding_mask)
                
                fractions = [core_cs_fraction, cs_fraction, d_fraction]
                domains = ['core_cs', 'cs', 'd']
                domain = domains[np.argmax(fractions)]
                
                if domain == 'core_cs':
                    return 'red'
                elif domain == 'cs':
                    return 'orange'
                else:
                    return 'cyan'
            else:
                cs_fraction = np.sum(domain_info['cs_rich'] & surrounding_mask) / np.sum(surrounding_mask)
                return 'orange' if cs_fraction > 0.5 else 'cyan'
        
        return 'gray'
    
    def _add_legend(self, ax, domain_info):
        """Add legend to the plot."""
        legend_elements = [
            plt.Line2D([], [], color='yellow', linewidth=1.5, linestyle='dashed',
                      alpha=0.8, label='Core-CS rich Region'),
        ]
        
        if domain_info and 'core_cs_rich' in domain_info:
            legend_elements.append(
                plt.Line2D([], [], color='orange', linewidth=1.5,
                          alpha=0.7, label='CS rich Region')
            )
        
        legend_elements.extend([
            plt.scatter([], [], color='white', alpha=0.3, s=10, label='Protein Core'),
            plt.scatter([], [], color='cyan', alpha=0.3, s=10, label='Protein Interface'),
            plt.scatter([], [], color='red', s=200, edgecolors='black',
                       linewidths=2, label='Protein in Core-CS'),
            plt.scatter([], [], color='orange', s=200, edgecolors='black',
                       linewidths=2, label='Protein in CS'),
            plt.scatter([], [], color='cyan', s=200, edgecolors='black',
                       linewidths=2, label='Protein in D'),
            plt.Line2D([], [], color='yellow', linestyle='--',
                      alpha=0.5, label=f'Analysis Shell ({CHOL_SHELL_RADIUS}Å)')
        ])
        
        ax.legend(handles=legend_elements, loc='upper right',
                 fontsize=8, bbox_to_anchor=(1.15, 1))
    
    def _add_frame_info(self, ax, frame_number, domain_info):
        """Add frame information text."""
        frame_info = [f'Frame: {frame_number}']
        
        if domain_info and 'domain_stats' in domain_info:
            stats = domain_info['domain_stats']
            
            if 'area_fraction_core_cs' in stats:
                frame_info.append(f"Core-CS rich: {stats['area_fraction_core_cs']*100:.1f}%")
            
            frame_info.extend([
                f"CS rich region: {stats['area_fraction_cs']*100:.1f}%",
                f"D rich region: {stats['area_fraction_d']*100:.1f}%"
            ])
        
        for i, text in enumerate(frame_info):
            if text:
                ax.text(0.02, 0.98 - i*0.04, text,
                       transform=ax.transAxes,
                       color='white',
                       fontsize=10,
                       verticalalignment='top')