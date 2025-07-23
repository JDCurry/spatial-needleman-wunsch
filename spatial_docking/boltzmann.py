# boltzmann.py
"""Boltzmann-weighted probabilistic scoring of spatial docking ensembles."""

import numpy as np

def boltzmann_ensemble_docking(spatial_alignment_all_paths, cavity_grid, molecule_grid, temperature=300):
    kT = 0.0019872041 * temperature
    all_scores = spatial_alignment_all_paths(cavity_grid, molecule_grid)
    boltzmann_weights = {}
    Z = 0.0

    for position, score in all_scores.items():
        energy = -score
        weight = np.exp(-energy / kT)
        boltzmann_weights[position] = weight
        Z += weight

    probabilities = {pos: w / Z for pos, w in boltzmann_weights.items()}
    return probabilities