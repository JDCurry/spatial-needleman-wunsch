# alignment.py
"""Core dynamic programming logic for spatial Needleman-Wunsch docking."""

import numpy as np

def calculate_placement_score(cavity_grid, molecule_grid, offset, compatibility_matrix):
    score = 0.0
    overlap_count = 0

    for mol_pos, mol_voxel in molecule_grid.items():
        cavity_pos = tuple(mol_pos[i] + offset[i] for i in range(3))
        if cavity_pos in cavity_grid:
            cavity_voxel = cavity_grid[cavity_pos]
            pair = (mol_voxel.property_type, cavity_voxel.property_type)
            reverse_pair = (cavity_voxel.property_type, mol_voxel.property_type)
            if pair in compatibility_matrix:
                score += compatibility_matrix[pair]
            elif reverse_pair in compatibility_matrix:
                score += compatibility_matrix[reverse_pair]
            overlap_count += 1

    if overlap_count > 0:
        score = score / overlap_count * np.sqrt(overlap_count)
    return score

def spatial_alignment(cavity_grid, molecule_grid, compatibility_matrix, gap_penalty=-1.0):
    cavity_positions = list(cavity_grid.keys())
    best_score = float('-inf')
    best_offset = (0, 0, 0)

    for cavity_pos in cavity_positions:
        placement_score = calculate_placement_score(cavity_grid, molecule_grid, cavity_pos, compatibility_matrix)
        score = placement_score + gap_penalty
        if score > best_score:
            best_score = score
            best_offset = cavity_pos

    return best_score, best_offset
