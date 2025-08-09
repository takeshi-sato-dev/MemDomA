"""
MemDomA - Membrane Domain Analysis Package

A comprehensive toolkit for analyzing membrane domains in MARTINI MD simulations.
"""

__version__ = "0.1.0"

from .analysis import MemDomA
from .core.trajectory import TrajectoryHandler
from .core.lipids import LipidAnalyzer
from .core.proteins import ProteinAnalyzer
from .core.domains import DomainAnalyzer

__all__ = [
    'MemDomA',
    'TrajectoryHandler',
    'LipidAnalyzer',
    'ProteinAnalyzer',
    'DomainAnalyzer'
]
