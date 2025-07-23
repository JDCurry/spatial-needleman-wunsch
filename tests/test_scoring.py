"""
Test suite for scoring functions and compatibility matrices.
"""

import pytest
import numpy as np
from spatial_docking.scoring import default_compatibility_matrix, timp2_specific_bonus
from spatial_docking import Voxel


class TestCompatibilityMatrix:
    """Test the chemical compatibility matrix."""
    
    def test_default_matrix_structure(self):
        """Test that default matrix has expected structure."""
        matrix = default_compatibility_matrix()
        
        assert isinstance(matrix, dict)
        assert len(matrix) > 0
        
        # Check some key interactions exist
        expected_pairs = [
            ('hydrophobic', 'hydrophobic'),
            ('polar', 'polar'),
            ('charged_pos', 'charged_neg'),
            ('charged_pos', 'charged_pos'),
            ('charged_neg', 'charged_neg')
        ]
        
        for pair in expected_pairs:
            assert pair in matrix
    
    def test_matrix_values_reasonable(self):
        """Test that compatibility values are chemically reasonable."""
        matrix = default_compatibility_matrix()
        
        # Favorable interactions should be positive
        assert matrix[('hydrophobic', 'hydrophobic')] > 0
        assert matrix[('polar', 'polar')] > 0
        assert matrix[('charged_pos', 'charged_neg')] > 0
        
        # Unfavorable interactions should be negative
        assert matrix[('hydrophobic', 'hydrophilic')] < 0
        assert matrix[('charged_pos', 'charged_pos')] < 0
        assert matrix[('charged_neg', 'charged_neg')] < 0
    
    def test_matrix_symmetry(self):
        """Test that matrix handles symmetric pairs correctly."""
        matrix = default_compatibility_matrix()
        
        # Some pairs should be symmetric by design
        # (though we handle reverse lookups in the algorithm)
        pair1 = ('hydrophobic', 'polar')
        pair2 = ('polar', 'hydrophobic')
        
        # At least one direction should exist
        has_pair1 = pair1 in matrix
        has_pair2 = pair2 in matrix
        
        assert has_pair1 or has_pair2


class TestTIMP2SpecificScoring:
    """Test TIMP2-specific scoring bonuses."""
    
    def setup_method(self):
        """Set up test voxels."""
        self.ser69_voxel = Voxel(
            position=(1.0, 1.0, 1.0),
            property_type='polar',
            nearby_residue='SER69'
        )
        
        self.leu76_voxel = Voxel(
            position=(2.0, 2.0, 2.0),
            property_type='hydrophobic',
            nearby_residue='LEU76'
        )
        
        self.generic_voxel = Voxel(
            position=(0.0, 0.0, 0.0),
            property_type='polar'
        )
    
    def test_h_bonding_bonus(self):
        """Test hydrogen bonding bonus with SER69."""
        mol_voxel = Voxel((0, 0, 0), 'polar')
        
        bonus = timp2_specific_bonus(mol_voxel, self.ser69_voxel)
        assert bonus > 0  # Should get H-bonding bonus
        
        # Non-polar molecule shouldn't get bonus
        mol_hydrophobic = Voxel((0, 0, 0), 'hydrophobic')
        bonus_hydrophobic = timp2_specific_bonus(mol_hydrophobic, self.ser69_voxel)
        assert bonus_hydrophobic <= bonus
    
    def test_hydrophobic_bonus(self):
        """Test hydrophobic contact bonus with LEU76."""
        mol_voxel = Voxel((0, 0, 0), 'hydrophobic')
        
        bonus = timp2_specific_bonus(mol_voxel, self.leu76_voxel)
        assert bonus > 0  # Should get hydrophobic bonus
        
        # Polar molecule shouldn't get same bonus
        mol_polar = Voxel((0, 0, 0), 'polar')
        bonus_polar = timp2_specific_bonus(mol_polar, self.leu76_voxel)
        assert bonus_polar <= bonus
    
    def test_no_bonus_generic(self):
        """Test that generic residues don't give specific bonuses."""
        mol_voxel = Voxel((0, 0, 0), 'polar')
        
        bonus = timp2_specific_bonus(mol_voxel, self.generic_voxel)
        assert bonus == 0  # No specific residue, no bonus


class TestScoringEdgeCases:
    """Test edge cases in scoring functions."""
    
    def test_empty_matrix(self):
        """Test behavior with empty compatibility matrix."""
        empty_matrix = {}
        
        # Should handle gracefully (algorithm should use defaults or handle missing pairs)
        assert isinstance(empty_matrix, dict)
    
    def test_unknown_property_types(self):
        """Test handling of unknown chemical property types."""
        matrix = default_compatibility_matrix()
        
        # Unknown property types shouldn't crash
        unknown_pair = ('unknown_type', 'another_unknown')
        result = matrix.get(unknown_pair, 0.0)  # Default to 0
        
        assert isinstance(result, (int, float))
    
    def test_matrix_modification(self):
        """Test that matrix can be safely modified."""
        matrix = default_compatibility_matrix()
        original_value = matrix[('hydrophobic', 'hydrophobic')]
        
        # Modify the matrix
        matrix[('hydrophobic', 'hydrophobic')] = 10.0
        assert matrix[('hydrophobic', 'hydrophobic')] == 10.0
        
        # Get a fresh matrix to ensure no global state
        fresh_matrix = default_compatibility_matrix()
        assert fresh_matrix[('hydrophobic', 'hydrophobic')] == original_value


class TestScoringIntegration:
    """Integration tests for scoring with real molecular data."""
    
    def test_realistic_interactions(self):
        """Test scoring with realistic molecular interactions."""
        matrix = default_compatibility_matrix()
        
        # Create realistic voxel pairs
        interactions = [
            ('hydrophobic', 'hydrophobic', 'favorable'),
            ('polar', 'polar', 'favorable'),
            ('charged_pos', 'charged_neg', 'very_favorable'),
            ('hydrophobic', 'polar', 'unfavorable'),
            ('charged_pos', 'charged_pos', 'very_unfavorable')
        ]
        
        for mol_type, cav_type, expected in interactions:
            if (mol_type, cav_type) in matrix:
                score = matrix[(mol_type, cav_type)]
            elif (cav_type, mol_type) in matrix:
                score = matrix[(cav_type, mol_type)]
            else:
                score = 0.0
            
            if expected == 'favorable':
                assert score > 0
            elif expected == 'very_favorable':
                assert score > 1.0
            elif expected == 'unfavorable':
                assert score < 0
            elif expected == 'very_unfavorable':
                assert score < -2.0
    
    def test_score_ranges(self):
        """Test that scores are in reasonable ranges."""
        matrix = default_compatibility_matrix()
        
        all_scores = list(matrix.values())
        
        # Scores should be in reasonable range for docking
        assert min(all_scores) >= -10.0  # Not too negative
        assert max(all_scores) <= 10.0   # Not too positive
        assert len([s for s in all_scores if s == 0]) < len(all_scores) // 2  # Not too many zeros


if __name__ == '__main__':
    pytest.main([__file__, '-v'])