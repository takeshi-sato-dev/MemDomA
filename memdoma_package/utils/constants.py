"""
Constants and configuration for membrane analysis.
"""

# Analysis parameters
DEFAULT_START = 20000
DEFAULT_STOP = 80000
DEFAULT_STEP = 100

# Spatial parameters
RADIUS = 10
MAX_GRID_SIZE = 100
HELIX_RADIUS = 3.0
INTERFACE_WIDTH = 2.0
GRID_SPACING = 1.0
CHOL_SHELL_RADIUS = 12.0

# Processing parameters
BATCH_SIZE = 50
DEFAULT_CPU_FRACTION = 0.75

# Domain analysis parameters
SMOOTHING_SIGMA = 1.5
CORE_THRESHOLD_STD = 1.0
CS_THRESHOLD_STD = 0.6

# Visualization parameters
DEFAULT_FPS = 5
DEFAULT_DPI = 200
DEFAULT_BITRATE = 1800
FIGURE_SIZE = (12, 10)

# Lipid types
STANDARD_LIPIDS = ["DPSM", "DIPC", "DPG3"]
ANALYSIS_LIPIDS = ["CHOL", "DIPC", "DPSM", "DPG3"]

# Chain selection patterns
CHAIN_A_PATTERN = 'name ??A'
CHAIN_B_PATTERN = 'name ??B'