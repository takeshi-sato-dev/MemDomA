"""
Trajectory handling and MDAnalysis utilities.
"""

import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy as np
from typing import Dict, Tuple, Optional


class TrajectoryHandler:
    """Handles trajectory loading and basic operations."""
    
    def __init__(self, psf_file: str, xtc_file: str):
        """
        Initialize trajectory handler.
        
        Parameters
        ----------
        psf_file : str
            Path to PSF topology file
        xtc_file : str
            Path to XTC trajectory file
        """
        self.psf_file = psf_file
        self.xtc_file = xtc_file
        self.universe = None
        self.leaflet0 = None
        
    def load_universe(self) -> mda.Universe:
        """Load trajectory into MDAnalysis Universe."""
        print("Loading trajectory...")
        self.universe = mda.Universe(self.psf_file, self.xtc_file)
        print("Trajectory loaded successfully.")
        return self.universe
    
    def identify_lipid_leaflets(self) -> mda.AtomGroup:
        """
        Identify lipid leaflets using LeafletFinder.
        
        Returns
        -------
        mda.AtomGroup
            Atom group containing the upper leaflet
        """
        if self.universe is None:
            raise RuntimeError("Universe not loaded. Call load_universe() first.")
            
        print("Identifying lipid leaflets...")
        L = LeafletFinder(self.universe, "name GL1 GL2 AM1 AM2 ROH GM1 GM2")
        cutoff = L.update(10)
        self.leaflet0 = L.groups(0)
        print("Lipid leaflets identified.")
        return self.leaflet0
    
    def select_lipids_and_cholesterol(self, 
                                     lipid_types: list, 
                                     leaflet: Optional[mda.AtomGroup] = None) -> Dict[str, mda.AtomGroup]:
        """
        Select lipids including DPG3 and cholesterol.
        
        Parameters
        ----------
        lipid_types : list
            List of lipid residue names to select
        leaflet : mda.AtomGroup, optional
            Leaflet to select from. If None, uses self.leaflet0
            
        Returns
        -------
        dict
            Dictionary mapping lipid types to AtomGroups
        """
        if leaflet is None:
            leaflet = self.leaflet0
            
        if leaflet is None:
            raise RuntimeError("Leaflet not identified. Call identify_lipid_leaflets() first.")
            
        selections = {}
        
        # Select regular lipids
        for resname in lipid_types:
            selection = leaflet.select_atoms(f"resname {resname}")
            selections[resname] = selection
            
        # Separate cholesterol selection
        selections['CHOL'] = leaflet.select_atoms("resname CHOL")
        
        return selections
    
    def select_proteins(self, protein_definitions: Dict[str, str]) -> Dict[str, mda.AtomGroup]:
        """
        Select protein regions based on selection strings.
        
        Parameters
        ----------
        protein_definitions : dict
            Dictionary mapping protein names to selection strings
            
        Returns
        -------
        dict
            Dictionary mapping protein names to AtomGroups
        """
        if self.universe is None:
            raise RuntimeError("Universe not loaded. Call load_universe() first.")
            
        proteins = {}
        for name, selection_string in protein_definitions.items():
            proteins[name] = self.universe.select_atoms(selection_string)
            
        return proteins
    
    def get_frame_data(self, frame_number: int) -> Tuple[np.ndarray, float]:
        """
        Get frame data including box dimensions and time.
        
        Parameters
        ----------
        frame_number : int
            Frame number to retrieve
            
        Returns
        -------
        tuple
            Box dimensions and time for the frame
        """
        if self.universe is None:
            raise RuntimeError("Universe not loaded. Call load_universe() first.")
            
        self.universe.trajectory[frame_number]
        box_dimensions = self.universe.dimensions[:2]
        time = self.universe.trajectory.time
        
        return box_dimensions, time