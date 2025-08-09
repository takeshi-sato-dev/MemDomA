"""Core analysis modules."""
from .trajectory import TrajectoryHandler
from .lipids import LipidAnalyzer
from .proteins import ProteinAnalyzer
from .domains import DomainAnalyzer

__all__ = ['TrajectoryHandler', 'LipidAnalyzer', 'ProteinAnalyzer', 'DomainAnalyzer']
