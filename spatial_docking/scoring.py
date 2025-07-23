# scoring.py
"""Chemical compatibility matrices and TIMP2-specific scoring functions."""

def default_compatibility_matrix():
    return {
        ('hydrophobic', 'hydrophobic'): 2.0,
        ('hydrophobic', 'hydrophilic'): -3.0,
        ('hydrophobic', 'polar'): -1.0,
        ('hydrophobic', 'empty'): -0.5,
        ('hydrophilic', 'hydrophilic'): 1.5,
        ('hydrophilic', 'polar'): 1.0,
        ('hydrophilic', 'charged_pos'): 0.5,
        ('hydrophilic', 'charged_neg'): 0.5,
        ('polar', 'polar'): 1.5,
        ('polar', 'charged_pos'): 1.0,
        ('polar', 'charged_neg'): 1.0,
        ('charged_pos', 'charged_neg'): 3.0,
        ('charged_pos', 'charged_pos'): -4.0,
        ('charged_neg', 'charged_neg'): -4.0,
        ('empty', 'empty'): 0.0,
    }

def timp2_specific_bonus(mol_voxel, cavity_voxel):
    bonus = 0.0
    if cavity_voxel.nearby_residue in ['SER69', 'TYR73'] and mol_voxel.property_type in ['polar', 'charged_pos', 'charged_neg']:
        bonus += 2.5
    if cavity_voxel.nearby_residue == 'LEU76' and mol_voxel.property_type == 'hydrophobic':
        bonus += 2.0
    if cavity_voxel.property_type == 'hydrophobic' and mol_voxel.property_type in ['charged_pos', 'charged_neg']:
        bonus -= 3.0
    return bonus
