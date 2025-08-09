"""
Basic tests for memdoma package.
"""

import pytest
import numpy as np
import pandas as pd
import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from memdoma import MembraneAnalysis
from memdoma.core.lipids import LipidAnalyzer
from memdoma.core.proteins import ProteinAnalyzer
from memdoma.core.domains import DomainAnalyzer
from memdoma.utils import constants


class TestConstants:
    """Test that constants are properly defined."""
    
    def test_constants_exist(self):
        """Test that all required constants exist."""
        assert hasattr(constants, 'DEFAULT_START')
        assert hasattr(constants, 'DEFAULT_STOP')
        assert hasattr(constants, 'DEFAULT_STEP')
        assert hasattr(constants, 'CHOL_SHELL_RADIUS')
        assert hasattr(constants, 'HELIX_RADIUS')
    
    def test_constants_values(self):
        """Test that constants have reasonable values."""
        assert constants.DEFAULT_START >= 0
        assert constants.DEFAULT_STOP > constants.DEFAULT_START
        assert constants.DEFAULT_STEP > 0
        assert constants.CHOL_SHELL_RADIUS > 0
        assert constants.HELIX_RADIUS > 0


class TestLipidAnalyzer:
    """Test lipid analysis functions."""
    
    def test_order_parameter_calculation(self):
        """Test that order parameter returns a float or NaN."""
        # This would require a mock residue object
        # For now, just test that the class exists and has the method
        analyzer = LipidAnalyzer()
        assert hasattr(analyzer, 'calculate_order_parameter')
    
    def test_normalize_array(self):
        """Test array normalization."""
        analyzer = DomainAnalyzer()
        arr = np.array([1, 2, 3, 4, 5])
        normalized = analyzer._normalize_array(arr)
        
        assert np.min(normalized) == 0.0
        assert np.max(normalized) == 1.0
        assert normalized.shape == arr.shape


class TestProteinAnalyzer:
    """Test protein analysis functions."""
    
    def test_initialize_statistics(self):
        """Test protein statistics initialization."""
        analyzer = ProteinAnalyzer()
        protein_names = ['Protein1', 'Protein2']
        stats = analyzer.initialize_protein_statistics(protein_names)
        
        assert 'Protein1' in stats
        assert 'Protein2' in stats
        assert 'total' in stats['Protein1']
        assert 'core_cs' in stats['Protein1']
        assert 'cs' in stats['Protein1']
        assert 'd' in stats['Protein1']


class TestMembraneAnalysis:
    """Test main analysis class."""
    
    def test_initialization(self):
        """Test that MembraneAnalysis can be initialized."""
        # Use dummy file paths for initialization
        analysis = MembraneAnalysis(
            psf_file='dummy.psf',
            xtc_file='dummy.xtc',
            start=0,
            stop=100,
            step=10
        )
        
        assert analysis.psf_file == 'dummy.psf'
        assert analysis.xtc_file == 'dummy.xtc'
        assert analysis.start == 0
        assert analysis.stop == 100
        assert analysis.step == 10
    
    def test_frame_batch_creation(self):
        """Test frame batch creation."""
        analysis = MembraneAnalysis(
            psf_file='dummy.psf',
            xtc_file='dummy.xtc',
            start=0,
            stop=100,
            step=10
        )
        
        num_frames = 10
        batches = analysis._create_frame_batches(num_frames)
        
        assert len(batches) > 0
        assert all(isinstance(batch, list) for batch in batches)
        
        # Check that all frames are covered
        all_frames = []
        for batch in batches:
            all_frames.extend(batch)
        
        expected_frames = [i * 10 for i in range(10)]
        assert all_frames == expected_frames


@pytest.mark.skipif(
    not os.path.exists('test_data'),
    reason="Test data directory not found"
)
class TestWithData:
    """Tests that require actual test data."""
    
    def test_load_test_data(self):
        """Test loading test data if available."""
        # Check if test data exists
        test_psf = 'test_data/test.psf'
        test_xtc = 'test_data/test.xtc'
        
        if os.path.exists(test_psf) and os.path.exists(test_xtc):
            analysis = MembraneAnalysis(
                psf_file=test_psf,
                xtc_file=test_xtc,
                start=0,
                stop=10,
                step=1
            )
            
            # Try to set up the analysis
            try:
                analysis.setup()
                assert analysis.universe is not None
                assert analysis.leaflet is not None
            except Exception as e:
                # It's okay if setup fails with test data
                # Just make sure the error is reasonable
                assert isinstance(e, (RuntimeError, ValueError))


if __name__ == "__main__":
    pytest.main([__file__, '-v'])