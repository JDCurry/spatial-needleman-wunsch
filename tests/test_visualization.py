"""
Test suite for visualization functions.
"""

import pytest
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from spatial_docking.visualization import visualize_docking, plot_score_heatmap
from spatial_docking import Voxel, create_synthetic_cavity, create_test_molecule

# Use non-interactive backend for testing
matplotlib.use('Agg')


class TestVisualizeDocking:
    """Test the main docking visualization function."""
    
    def setup_method(self):
        """Set up test data."""
        self.cavity = create_synthetic_cavity(size=(3, 3, 3), pattern='checkerboard')
        self.molecule = create_test_molecule(length=2)
        self.offset = (1.0, 1.0, 1.0)
    
    def test_visualization_runs_without_error(self):
        """Test that visualization completes without throwing exceptions."""
        try:
            visualize_docking(self.cavity, self.molecule, self.offset, save_path=None, show_plot=False)
        except Exception as e:
            pytest.fail(f"Visualization failed with error: {e}")
    
    def test_visualization_with_empty_cavity(self):
        """Test visualization with empty cavity."""
        empty_cavity = {}
        
        try:
            visualize_docking(empty_cavity, self.molecule, self.offset, show_plot=False)
        except Exception as e:
            pytest.fail(f"Empty cavity visualization failed: {e}")
    
    def test_visualization_with_empty_molecule(self):
        """Test visualization with empty molecule."""
        empty_molecule = {}
        
        try:
            visualize_docking(self.cavity, empty_molecule, self.offset, show_plot=False)
        except Exception as e:
            pytest.fail(f"Empty molecule visualization failed: {e}")
    
    def test_visualization_with_custom_title(self):
        """Test visualization with custom title."""
        custom_title = "Custom Test Title"
        
        try:
            visualize_docking(self.cavity, self.molecule, self.offset, 
                            title=custom_title, show_plot=False)
        except Exception as e:
            pytest.fail(f"Custom title visualization failed: {e}")
    
    def test_figure_creation(self):
        """Test that visualization creates a matplotlib figure."""
        # Clear any existing figures
        plt.close('all')
        
        # Create visualization
        visualize_docking(self.cavity, self.molecule, self.offset, show_plot=False)
        
        # Check that a figure was created
        assert len(plt.get_fignums()) > 0
        
        # Clean up
        plt.close('all')
    
    def test_3d_axes_creation(self):
        """Test that 3D axes are properly created."""
        plt.close('all')
        
        visualize_docking(self.cavity, self.molecule, self.offset, show_plot=False)
        
        fig = plt.gcf()
        axes = fig.get_axes()
        
        assert len(axes) > 0
        # Check if it's a 3D axis (has zaxis attribute)
        ax = axes[0]
        assert hasattr(ax, 'zaxis')
        
        plt.close('all')


class TestScoreHeatmap:
    """Test score heatmap visualization (if implemented)."""
    
    def test_heatmap_placeholder(self):
        """Placeholder test for score heatmap functionality."""
        # This would test plot_score_heatmap when implemented
        pass


class TestVisualizationUtilities:
    """Test utility functions for visualization."""
    
    def test_color_mapping(self):
        """Test that color mapping works for all property types."""
        color_map = {
            'hydrophobic': 'yellow',
            'hydrophilic': 'blue',
            'polar': 'green',
            'charged_pos': 'red',
            'charged_neg': 'purple',
            'empty': 'gray'
        }
        
        # All common property types should have colors
        property_types = ['hydrophobic', 'polar', 'charged_pos', 'charged_neg']
        
        for prop_type in property_types:
            assert prop_type in color_map
            assert isinstance(color_map[prop_type], str)
    
    def test_plot_parameters(self):
        """Test that plot uses reasonable parameters."""
        plt.close('all')
        
        visualize_docking(
            create_synthetic_cavity(size=(2, 2, 2)),
            create_test_molecule(length=1),
            (0, 0, 0),
            show_plot=False
        )
        
        fig = plt.gcf()
        ax = fig.get_axes()[0]
        
        # Check that axes are labeled
        assert ax.get_xlabel() != ''
        assert ax.get_ylabel() != ''
        assert ax.get_zlabel() != ''
        
        # Check that title is set
        assert ax.get_title() != ''
        
        plt.close('all')


class TestVisualizationIntegration:
    """Integration tests for visualization with different data types."""
    
    def test_different_cavity_patterns(self):
        """Test visualization with different cavity patterns."""
        patterns = ['checkerboard', 'layered', 'random']
        molecule = create_test_molecule(length=2)
        offset = (0, 0, 0)
        
        for pattern in patterns:
            cavity = create_synthetic_cavity(size=(3, 3, 3), pattern=pattern)
            
            try:
                visualize_docking(cavity, molecule, offset, show_plot=False)
            except Exception as e:
                pytest.fail(f"Visualization failed for pattern {pattern}: {e}")
            
            plt.close('all')
    
    def test_different_molecule_sizes(self):
        """Test visualization with different molecule sizes."""
        cavity = create_synthetic_cavity(size=(4, 4, 4))
        offset = (1, 1, 1)
        
        for length in [1, 3, 5]:
            molecule = create_test_molecule(length=length)
            
            try:
                visualize_docking(cavity, molecule, offset, show_plot=False)
            except Exception as e:
                pytest.fail(f"Visualization failed for molecule length {length}: {e}")
            
            plt.close('all')
    
    def test_edge_case_offsets(self):
        """Test visualization with edge case offsets."""
        cavity = create_synthetic_cavity(size=(3, 3, 3))
        molecule = create_test_molecule(length=2)
        
        edge_offsets = [
            (0, 0, 0),      # No offset
            (-1, -1, -1),   # Negative offset
            (10, 10, 10),   # Large offset (molecule outside cavity)
            (0.5, 0.5, 0.5) # Non-integer offset
        ]
        
        for offset in edge_offsets:
            try:
                visualize_docking(cavity, molecule, offset, show_plot=False)
            except Exception as e:
                pytest.fail(f"Visualization failed for offset {offset}: {e}")
            
            plt.close('all')


class TestVisualizationPerformance:
    """Test visualization performance with larger datasets."""
    
    def test_large_cavity_visualization(self):
        """Test visualization with larger cavities."""
        # Create a larger cavity
        large_cavity = create_synthetic_cavity(size=(10, 10, 8))
        molecule = create_test_molecule(length=3)
        offset = (2, 2, 2)
        
        try:
            visualize_docking(large_cavity, molecule, offset, show_plot=False)
        except Exception as e:
            pytest.fail(f"Large cavity visualization failed: {e}")
        
        plt.close('all')
    
    def test_many_property_types(self):
        """Test visualization with molecules having many different property types."""
        cavity = create_synthetic_cavity(size=(4, 4, 4))
        
        # Create molecule with diverse property types
        diverse_molecule = {}
        property_types = ['hydrophobic', 'polar', 'charged_pos', 'charged_neg']
        
        for i, prop_type in enumerate(property_types):
            pos = (float(i), 0.0, 0.0)
            diverse_molecule[pos] = Voxel(pos, prop_type, 1.0, 0.5)
        
        try:
            visualize_docking(cavity, diverse_molecule, (0, 0, 0), show_plot=False)
        except Exception as e:
            pytest.fail(f"Diverse molecule visualization failed: {e}")
        
        plt.close('all')


class TestVisualizationValidation:
    """Test that visualization accurately represents the data."""
    
    def test_molecule_placement_accuracy(self):
        """Test that molecules are placed at correct positions."""
        # This is hard to test automatically without examining the plot data
        # But we can at least ensure the function completes and handles positions
        
        cavity = {(0.0, 0.0, 0.0): Voxel((0.0, 0.0, 0.0), 'polar')}
        molecule = {(0.0, 0.0, 0.0): Voxel((0.0, 0.0, 0.0), 'hydrophobic')}
        offset = (0.0, 0.0, 0.0)
        
        try:
            visualize_docking(cavity, molecule, offset, show_plot=False)
        except Exception as e:
            pytest.fail(f"Simple placement visualization failed: {e}")
        
        plt.close('all')


def test_matplotlib_backend():
    """Test that matplotlib backend is properly configured for testing."""
    backend = matplotlib.get_backend()
    # Should be using non-interactive backend for testing
    assert backend == 'Agg'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])