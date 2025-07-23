"""
Spatial Needleman-Wunsch: A Deterministic Dynamic Programming Framework for 3D Molecular Docking

This package provides deterministic, interpretable molecular docking using principles
from sequence alignment extended to three-dimensional molecular space.

Key Components:
    - SpatialDocking: Main docking algorithm class
    - scoring: Chemical compatibility matrices and scoring functions
    - visualization: 3D plotting and result analysis
    - boltzmann: Ensemble docking with thermodynamic weighting
    - pareto: Multi-objective optimization
    - conformers: Flexible ligand sampling
    - adaptive: Machine learning-based scoring optimization

Example:
    >>> from spatial_docking import SpatialDocking
    >>> from spatial_docking.scoring import default_compatibility_matrix
    >>> 
    >>> docker = SpatialDocking(grid_spacing=0.5)
    >>> score, pose = docker.align(cavity, molecule, 
    ...                           compatibility_matrix=default_compatibility_matrix())
    >>> docker.visualize(cavity, molecule, pose)
"""

__version__ = "1.0.0"
__author__ = "Joshua D. Curry"
__email__ = "jcurry3428@smail.pcd.edu"
__license__ = "MIT"

# Core imports
from .core import SpatialDocking, Voxel
from .alignment import spatial_alignment, calculate_placement_score
from .scoring import default_compatibility_matrix
from .visualization import visualize_docking, plot_score_heatmap
from .utils import load_molecule, load_cavity, voxelize_atoms

# Advanced features
from .boltzmann import boltzmann_ensemble_docking
from .pareto import pareto_optimal_docking
from .conformers import generate_torsion_states, ConformationalDockingState
from .adaptive import AdaptiveScoringFunction

__all__ = [
    # Core classes and functions
    'SpatialDocking',
    'Voxel',
    'spatial_alignment',
    'calculate_placement_score',
    
    # Scoring and compatibility
    'default_compatibility_matrix',

    
    # Visualization
    'visualize_docking',
    'plot_score_heatmap',
    
    # Utilities
    'load_molecule',
    'load_cavity',
    'voxelize_atoms',
    
    # Advanced features
    'boltzmann_ensemble_docking',
    'pareto_optimal_docking',
    'generate_torsion_states',
    'ConformationalDockingState',
    'AdaptiveScoringFunction',
]

# Package metadata
__package_name__ = "spatial-needleman-wunsch"
__description__ = "A deterministic dynamic programming framework for 3D molecular docking"
__url__ = "https://github.com/JDCurry/spatial-needleman-wunsch"