"""
Test suite for Spatial Needleman-Wunsch docking framework.

This module contains unit tests for all core functionality including:
- Voxel creation and validation
- Spatial alignment algorithms
- Scoring matrix operations
- Visualization components
- Result analysis
"""

import pytest
import numpy as np
from spatial_docking import (
    SpatialDocking, Voxel, DockingResult,
    create_synthetic_cavity, create_test_molecule,
    spatial_alignment, calculate_placement_score
)
from spatial_docking.scoring import default_compatibility_matrix


class TestVoxel:
    """Test the Voxel dataclass."""
    
    def test_voxel_creation(self):
        """Test basic voxel creation."""
        voxel = Voxel(
            position=(1.0, 2.0, 3.0),
            property_type='hydrophobic',
            accessibility=0.8,
            flexibility=0.6
        )
        
        assert voxel.position == (1.0, 2.0, 3.0)
        assert voxel.property_type == 'hydrophobic'
        assert voxel.accessibility == 0.8
        assert voxel.flexibility == 0.6
        assert voxel.nearby_residue is None
    
    def test_voxel_validation(self):
        """Test voxel parameter validation."""
        # Test invalid accessibility
        with pytest.raises(ValueError, match="Accessibility must be between 0.0 and 1.0"):
            Voxel((0, 0, 0), 'polar', accessibility=1.5)
        
        # Test invalid flexibility
        with pytest.raises(ValueError, match="Flexibility must be between 0.0 and 1.0"):
            Voxel((0, 0, 0), 'polar', flexibility=-0.1)


class TestSyntheticData:
    """Test synthetic cavity and molecule generation."""
    
    def test_synthetic_cavity_creation(self):
        """Test creation of synthetic cavities."""
        cavity = create_synthetic_cavity(size=(5, 5, 5), pattern='checkerboard')
        
        assert len(cavity) == 125  # 5*5*5
        assert all(isinstance(v, Voxel) for v in cavity.values())
        
        # Check property types
        properties = [v.property_type for v in cavity.values()]
        assert 'hydrophobic' in properties
        assert 'polar' in properties
    
    def test_different_cavity_patterns(self):
        """Test different cavity patterns."""
        patterns = ['checkerboard', 'layered', 'random']
        
        for pattern in patterns:
            cavity = create_synthetic_cavity(size=(3, 3, 3), pattern=pattern)
            assert len(cavity) == 27
            assert all(isinstance(v, Voxel) for v in cavity.values())
    
    def test_test_molecule_creation(self):
        """Test creation of test molecules."""
        molecule = create_test_molecule(length=4)
        
        assert len(molecule) == 4
        assert all(isinstance(v, Voxel) for v in molecule.values())
        
        # Check linear arrangement
        positions = sorted(molecule.keys())
        for i, pos in enumerate(positions):
            assert pos == (float(i), 0.0, 0.0)


class TestSpatialDocking:
    """Test the main SpatialDocking class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.docker = SpatialDocking(grid_spacing=0.5, verbose=False)
        self.cavity = create_synthetic_cavity(size=(5, 5, 5), pattern='checkerboard')
        self.molecule = create_test_molecule(length=3)
        self.compatibility_matrix = default_compatibility_matrix()
    
    def test_docking_initialization(self):
        """Test SpatialDocking initialization."""
        docker = SpatialDocking(
            grid_spacing=1.0,
            max_translation=8,
            gap_penalty=-2.0,
            verbose=True
        )
        
        assert docker.grid_spacing == 1.0
        assert docker.max_translation == 8
        assert docker.gap_penalty == -2.0
        assert docker.verbose is True
    
    def test_basic_alignment(self):
        """Test basic molecular alignment."""
        result = self.docker.align(self.cavity, self.molecule)
        
        assert isinstance(result, DockingResult)
        assert isinstance(result.score, float)
        assert isinstance(result.translation, tuple)
        assert len(result.translation) == 3
        assert result.cavity_grid is not None
        assert result.molecule_grid is not None
    
    def test_reproducibility(self):
        """Test that docking is deterministic."""
        result1 = self.docker.align(self.cavity, self.molecule)
        result2 = self.docker.align(self.cavity, self.molecule)
        
        assert result1.score == result2.score
        assert result1.translation == result2.translation
    
    def test_multiple_molecule_docking(self):
        """Test docking multiple molecules."""
        molecules = [
            create_test_molecule(length=2),
            create_test_molecule(length=3),
            create_test_molecule(length=4)
        ]
        
        results = self.docker.align_multiple(self.cavity, molecules)
        
        assert len(results) == 3
        assert all(isinstance(r, DockingResult) for r in results)
        
        # Results should be sorted by score (best first)
        scores = [r.score for r in results]
        assert scores == sorted(scores, reverse=True)
    
    def test_pose_scoring(self):
        """Test scoring of specific poses."""
        translation = (1.0, 1.0, 1.0)
        score = self.docker.score_pose(
            self.cavity, self.molecule, translation, self.compatibility_matrix
        )
        
        assert isinstance(score, float)
    
    def test_binding_site_analysis(self):
        """Test binding site interaction analysis."""
        result = self.docker.align(self.cavity, self.molecule)
        analysis = self.docker.get_binding_site_analysis(result)
        
        assert 'total_score' in analysis
        assert 'interactions' in analysis
        assert 'interaction_counts' in analysis
        assert 'residue_contributions' in analysis
        
        assert analysis['total_score'] == result.score
        assert isinstance(analysis['interactions'], list)


class TestScoringFunctions:
    """Test scoring and compatibility matrix functions."""
    
    def test_default_compatibility_matrix(self):
        """Test default compatibility matrix."""
        matrix = default_compatibility_matrix()
        
        assert isinstance(matrix, dict)
        
        # Check some expected interactions
        assert ('hydrophobic', 'hydrophobic') in matrix
        assert ('polar', 'polar') in matrix
        assert ('charged_pos', 'charged_neg') in matrix
        
        # Hydrophobic-hydrophobic should be favorable
        assert matrix[('hydrophobic', 'hydrophobic')] > 0
        
        # Hydrophobic-hydrophilic should be unfavorable
        assert matrix[('hydrophobic', 'hydrophilic')] < 0
    
    def test_placement_score_calculation(self):
        """Test calculation of placement scores."""
        cavity = create_synthetic_cavity(size=(3, 3, 3), pattern='checkerboard')
        molecule = create_test_molecule(length=2)
        matrix = default_compatibility_matrix()
        
        score = calculate_placement_score(
            cavity, molecule, offset=(0, 0, 0), compatibility_matrix=matrix
        )
        
        assert isinstance(score, float)


class TestVisualization:
    """Test visualization components."""
    
    def test_visualization_no_errors(self):
        """Test that visualization runs without errors."""
        docker = SpatialDocking(verbose=False)
        cavity = create_synthetic_cavity(size=(3, 3, 3))
        molecule = create_test_molecule(length=2)
        
        result = docker.align(cavity, molecule)
        
        # This should not raise any exceptions
        try:
            docker.visualize(result, show_plot=False)
        except Exception as e:
            pytest.fail(f"Visualization failed with error: {e}")


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_empty_molecule(self):
        """Test behavior with empty molecule."""
        docker = SpatialDocking(verbose=False)
        cavity = create_synthetic_cavity(size=(3, 3, 3))
        empty_molecule = {}
        
        result = docker.align(cavity, empty_molecule)
        
        # Should handle gracefully
        assert isinstance(result, DockingResult)
        assert result.score <= 0  # Empty molecule shouldn't score well
    
    def test_empty_cavity(self):
        """Test behavior with empty cavity."""
        docker = SpatialDocking(verbose=False)
        empty_cavity = {}
        molecule = create_test_molecule(length=2)
        
        result = docker.align(empty_cavity, molecule)
        
        # Should handle gracefully
        assert isinstance(result, DockingResult)
    
    def test_large_translation_range(self):
        """Test with large translation search range."""
        docker = SpatialDocking(max_translation=20, verbose=False)
        cavity = create_synthetic_cavity(size=(3, 3, 3))
        molecule = create_test_molecule(length=1)
        
        result = docker.align(cavity, molecule)
        
        # Should complete without timeout
        assert isinstance(result, DockingResult)
    
    def test_invalid_visualization_input(self):
        """Test visualization with invalid input."""
        docker = SpatialDocking(verbose=False)
        
        # Create result without grids
        invalid_result = DockingResult(
            score=1.0,
            translation=(0, 0, 0),
            cavity_grid=None,
            molecule_grid=None
        )
        
        with pytest.raises(ValueError, match="must contain cavity_grid and molecule_grid"):
            docker.visualize(invalid_result)


class TestPerformance:
    """Test performance characteristics."""
    
    def test_scaling_behavior(self):
        """Test that algorithm scales reasonably with size."""
        docker = SpatialDocking(verbose=False)
        
        sizes = [(3, 3, 3), (5, 5, 5), (7, 7, 7)]
        molecule = create_test_molecule(length=2)
        
        for size in sizes:
            cavity = create_synthetic_cavity(size=size)
            result = docker.align(cavity, molecule)
            
            # Should complete for all sizes
            assert isinstance(result, DockingResult)
    
    def test_molecule_size_scaling(self):
        """Test scaling with different molecule sizes."""
        docker = SpatialDocking(verbose=False)
        cavity = create_synthetic_cavity(size=(5, 5, 5))
        
        for length in [1, 3, 5, 8]:
            molecule = create_test_molecule(length=length)
            result = docker.align(cavity, molecule)
            
            # Should complete for all molecule sizes
            assert isinstance(result, DockingResult)


# Integration test
def test_full_workflow():
    """Test complete docking workflow from start to finish."""
    # Create data
    cavity = create_synthetic_cavity(size=(6, 6, 4), pattern='layered')
    molecules = [
        create_test_molecule(length=2),
        create_test_molecule(length=4)
    ]
    
    # Initialize docking
    docker = SpatialDocking(grid_spacing=0.5, verbose=False)
    
    # Run multiple docking
    results = docker.align_multiple(cavity, molecules)
    
    # Analyze best result
    best_result = results[0]
    analysis = docker.get_binding_site_analysis(best_result)
    
    # Visualize (without showing)
    docker.visualize(best_result, show_plot=False)
    
    # All steps should complete successfully
    assert len(results) == 2
    assert isinstance(analysis, dict)
    assert 'interactions' in analysis


if __name__ == '__main__':
    # Run tests if script is executed directly
    pytest.main([__file__, '-v'])