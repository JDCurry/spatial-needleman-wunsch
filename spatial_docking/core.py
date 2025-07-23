"""
Core classes and functions for Spatial Needleman-Wunsch docking.

This module provides the main SpatialDocking class that orchestrates
the entire docking workflow, from voxelization to alignment to visualization.
"""

import numpy as np
from typing import Dict, Tuple, List, Optional, Any
from dataclasses import dataclass
import logging

from .alignment import spatial_alignment, calculate_placement_score
from .scoring import default_compatibility_matrix
from .visualization import visualize_docking
from .utils import voxelize_atoms

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class Voxel:
    """
    Represents a single voxel in 3D space with chemical properties.
    
    Attributes:
        position: (x, y, z) coordinates in Angstroms
        property_type: Chemical property ('hydrophobic', 'polar', 'charged_pos', etc.)
        accessibility: Solvent accessibility score [0.0, 1.0]
        flexibility: Conformational flexibility [0.0, 1.0]
        nearby_residue: Optional nearby protein residue identifier
    """
    position: Tuple[float, float, float]
    property_type: str
    accessibility: float = 1.0
    flexibility: float = 0.5
    nearby_residue: Optional[str] = None
    
    def __post_init__(self):
        """Validate voxel parameters."""
        if not 0.0 <= self.accessibility <= 1.0:
            raise ValueError(f"Accessibility must be between 0.0 and 1.0, got {self.accessibility}")
        if not 0.0 <= self.flexibility <= 1.0:
            raise ValueError(f"Flexibility must be between 0.0 and 1.0, got {self.flexibility}")


@dataclass
class DockingResult:
    """
    Container for docking results.
    
    Attributes:
        score: Optimal docking score
        translation: Best (x, y, z) translation vector
        rotation: Optional rotation quaternion (for future 6D extension)
        cavity_grid: Original cavity voxel grid
        molecule_grid: Original molecule voxel grid
        compatibility_matrix: Scoring matrix used
        metadata: Additional information about the docking run
    """
    score: float
    translation: Tuple[float, float, float]
    rotation: Optional[Tuple[float, float, float, float]] = None
    cavity_grid: Optional[Dict] = None
    molecule_grid: Optional[Dict] = None
    compatibility_matrix: Optional[Dict] = None
    metadata: Optional[Dict] = None


class SpatialDocking:
    """
    Main class for spatial Needleman-Wunsch molecular docking.
    
    This class orchestrates the entire docking workflow:
    1. Voxelization of cavity and molecule
    2. Spatial alignment using dynamic programming
    3. Result analysis and visualization
    
    Args:
        grid_spacing: Voxel size in Angstroms (default: 0.5)
        max_translation: Maximum translation range to search
        gap_penalty: Penalty for unfilled cavity space
        verbose: Enable detailed logging
    """
    
    def __init__(self, 
                 grid_spacing: float = 0.5,
                 max_translation: int = 10,
                 gap_penalty: float = -1.0,
                 verbose: bool = False):
        
        self.grid_spacing = grid_spacing
        self.max_translation = max_translation
        self.gap_penalty = gap_penalty
        self.verbose = verbose
        
        if verbose:
            logger.setLevel(logging.DEBUG)
        
        logger.info(f"Initialized SpatialDocking with grid_spacing={grid_spacing}Å")
    
    def align(self, 
              cavity_grid: Dict[Tuple[float, float, float], Voxel],
              molecule_grid: Dict[Tuple[float, float, float], Voxel],
              compatibility_matrix: Optional[Dict] = None) -> DockingResult:
        """
        Perform spatial alignment of molecule within cavity.
        
        Args:
            cavity_grid: Voxelized protein cavity
            molecule_grid: Voxelized ligand molecule
            compatibility_matrix: Chemical compatibility scoring matrix
            
        Returns:
            DockingResult containing optimal score and pose
        """
        
        if compatibility_matrix is None:
            compatibility_matrix = default_compatibility_matrix()
            logger.info("Using default compatibility matrix")
        
        logger.info(f"Aligning molecule ({len(molecule_grid)} voxels) "
                   f"in cavity ({len(cavity_grid)} voxels)")
        
        # Perform spatial alignment
        score, translation = spatial_alignment(
            cavity_grid=cavity_grid,
            molecule_grid=molecule_grid,
            compatibility_matrix=compatibility_matrix,
            max_translation=self.max_translation
        )
        
        logger.info(f"Optimal docking score: {score:.3f}")
        logger.info(f"Best translation: {translation}")
        
        # Create result object
        result = DockingResult(
            score=score,
            translation=translation,
            cavity_grid=cavity_grid,
            molecule_grid=molecule_grid,
            compatibility_matrix=compatibility_matrix,
            metadata={
                'grid_spacing': self.grid_spacing,
                'max_translation': self.max_translation,
                'gap_penalty': self.gap_penalty,
                'cavity_voxels': len(cavity_grid),
                'molecule_voxels': len(molecule_grid)
            }
        )
        
        return result
    
    def align_multiple(self,
                      cavity_grid: Dict[Tuple[float, float, float], Voxel],
                      molecules: List[Dict[Tuple[float, float, float], Voxel]],
                      compatibility_matrix: Optional[Dict] = None) -> List[DockingResult]:
        """
        Align multiple molecules against the same cavity.
        
        Args:
            cavity_grid: Voxelized protein cavity
            molecules: List of voxelized ligand molecules
            compatibility_matrix: Chemical compatibility scoring matrix
            
        Returns:
            List of DockingResult objects, one per molecule
        """
        
        results = []
        for i, molecule_grid in enumerate(molecules):
            logger.info(f"Docking molecule {i+1}/{len(molecules)}")
            result = self.align(cavity_grid, molecule_grid, compatibility_matrix)
            results.append(result)
        
        # Sort by score (best first)
        results.sort(key=lambda x: x.score, reverse=True)
        
        logger.info(f"Completed docking of {len(molecules)} molecules")
        logger.info(f"Best score: {results[0].score:.3f}")
        
        return results
    
    def visualize(self, 
                  result: DockingResult,
                  save_path: Optional[str] = None,
                  show_plot: bool = True) -> None:
        """
        Visualize docking result in 3D.
        
        Args:
            result: DockingResult from align() method
            save_path: Optional path to save figure
            show_plot: Whether to display the plot
        """
        
        if result.cavity_grid is None or result.molecule_grid is None:
            raise ValueError("Result must contain cavity_grid and molecule_grid for visualization")
        
        visualize_docking(
            cavity_grid=result.cavity_grid,
            molecule_grid=result.molecule_grid,
            offset=result.translation,
            save_path=save_path,
            show_plot=show_plot
        )
    
    def score_pose(self,
                   cavity_grid: Dict[Tuple[float, float, float], Voxel],
                   molecule_grid: Dict[Tuple[float, float, float], Voxel],
                   translation: Tuple[float, float, float],
                   compatibility_matrix: Optional[Dict] = None) -> float:
        """
        Score a specific molecular pose without optimization.
        
        Args:
            cavity_grid: Voxelized protein cavity
            molecule_grid: Voxelized ligand molecule  
            translation: (x, y, z) translation to score
            compatibility_matrix: Chemical compatibility scoring matrix
            
        Returns:
            Compatibility score for the given pose
        """
        
        if compatibility_matrix is None:
            compatibility_matrix = default_compatibility_matrix()
        
        return calculate_placement_score(
            cavity_grid=cavity_grid,
            molecule_grid=molecule_grid,
            offset=translation,
            compatibility_matrix=compatibility_matrix
        )
    
    def get_binding_site_analysis(self, result: DockingResult) -> Dict[str, Any]:
        """
        Analyze binding site interactions for a docking result.
        
        Args:
            result: DockingResult from align() method
            
        Returns:
            Dictionary with detailed binding site analysis
        """
        
        if result.cavity_grid is None or result.molecule_grid is None:
            raise ValueError("Result must contain cavity_grid and molecule_grid for analysis")
        
        analysis = {
            'total_score': result.score,
            'translation': result.translation,
            'interactions': [],
            'interaction_counts': {},
            'residue_contributions': {}
        }
        
        # Analyze individual voxel interactions
        for mol_pos, mol_voxel in result.molecule_grid.items():
            translated_pos = tuple(mol_pos[i] + result.translation[i] for i in range(3))
            
            if translated_pos in result.cavity_grid:
                cav_voxel = result.cavity_grid[translated_pos]
                
                # Record interaction
                interaction = {
                    'molecule_voxel': mol_pos,
                    'cavity_voxel': translated_pos,
                    'molecule_type': mol_voxel.property_type,
                    'cavity_type': cav_voxel.property_type,
                    'residue': cav_voxel.nearby_residue,
                    'score': result.compatibility_matrix.get(
                        (mol_voxel.property_type, cav_voxel.property_type), 0.0
                    )
                }
                analysis['interactions'].append(interaction)
                
                # Count interaction types
                pair = (mol_voxel.property_type, cav_voxel.property_type)
                analysis['interaction_counts'][pair] = analysis['interaction_counts'].get(pair, 0) + 1
                
                # Sum residue contributions
                if cav_voxel.nearby_residue:
                    res = cav_voxel.nearby_residue
                    analysis['residue_contributions'][res] = analysis['residue_contributions'].get(res, 0) + interaction['score']
        
        return analysis


def create_synthetic_cavity(size: Tuple[int, int, int] = (10, 10, 8),
                           pattern: str = 'checkerboard') -> Dict[Tuple[float, float, float], Voxel]:
    """
    Create a synthetic protein cavity for testing.
    
    Args:
        size: (x, y, z) dimensions in voxels
        pattern: Pattern type ('checkerboard', 'random', 'layered')
        
    Returns:
        Dictionary mapping positions to Voxel objects
    """
    
    cavity = {}
    
    for x in range(size[0]):
        for y in range(size[1]):
            for z in range(size[2]):
                pos = (float(x), float(y), float(z))
                
                if pattern == 'checkerboard':
                    # Alternating hydrophobic/polar pattern
                    if (x + y + z) % 2 == 0:
                        prop_type = 'hydrophobic'
                    else:
                        prop_type = 'polar'
                        
                elif pattern == 'layered':
                    # Hydrophobic core, polar surface
                    if z < size[2] // 3 or z > 2 * size[2] // 3:
                        prop_type = 'polar'
                    else:
                        prop_type = 'hydrophobic'
                        
                elif pattern == 'random':
                    # Random distribution
                    np.random.seed(x * 1000 + y * 100 + z)  # Reproducible randomness
                    prop_type = np.random.choice(['hydrophobic', 'polar', 'charged_pos', 'charged_neg'])
                    
                else:
                    raise ValueError(f"Unknown pattern: {pattern}")
                
                cavity[pos] = Voxel(
                    position=pos,
                    property_type=prop_type,
                    accessibility=0.5 + 0.5 * np.random.random(),
                    flexibility=0.3 + 0.4 * np.random.random()
                )
    
    return cavity


def create_test_molecule(length: int = 5) -> Dict[Tuple[float, float, float], Voxel]:
    """
    Create a linear test molecule.
    
    Args:
        length: Number of atoms in the molecule
        
    Returns:
        Dictionary mapping positions to Voxel objects
    """
    
    molecule = {}
    
    for i in range(length):
        pos = (float(i), 0.0, 0.0)
        
        # Alternating polar/hydrophobic pattern
        if i % 2 == 0:
            prop_type = 'polar'
        else:
            prop_type = 'hydrophobic'
        
        molecule[pos] = Voxel(
            position=pos,
            property_type=prop_type,
            accessibility=1.0,
            flexibility=0.8
        )
    
    return molecule