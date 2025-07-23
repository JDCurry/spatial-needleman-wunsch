# Spatial Needleman-Wunsch

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**A deterministic dynamic programming framework for 3D molecular docking**

Spatial Needleman-Wunsch extends the mathematical guarantees of classical sequence alignment to three-dimensional molecular space, providing **deterministic, interpretable, and optimal** molecular docking solutions.

## Why Deterministic Docking?

Current molecular docking methods rely on stochastic sampling and heuristic optimization, leading to:
- **Inconsistent results** across multiple runs
- **No guarantee of finding optimal solutions**
- **Limited interpretability** of scoring decisions

Spatial Needleman-Wunsch solves these problems by:
- ✅ **Guaranteed optimal alignment** given the scoring function
- ✅ **Reproducible results** every time
- ✅ **Complete interpretability** of every score component
- ✅ **Chemical intuition preservation** through explicit compatibility matrices

## Quick Start

```bash
# Install dependencies
pip install spatial-needleman-wunsch

# Or clone and install from source
git clone https://github.com/JDCurry/spatial-needleman-wunsch.git
cd spatial-needleman-wunsch
pip install -e .
```

```python
from spatial_docking import SpatialDocking
from spatial_docking.scoring import default_compatibility_matrix
from spatial_docking.utils import load_molecule, load_cavity

# Load your protein cavity and ligand molecule
cavity = load_cavity("examples/data/protein_cavity.json")
molecule = load_molecule("examples/data/ligand.sdf")

# Initialize the docking algorithm
docker = SpatialDocking(grid_spacing=0.5)

# Run deterministic alignment
score, pose = docker.align(cavity, molecule, 
                          compatibility_matrix=default_compatibility_matrix())

print(f"Optimal docking score: {score:.3f}")
print(f"Best translation: {pose['translation']}")

# Visualize results
docker.visualize(cavity, molecule, pose)
```

## Core Algorithm

The framework discretizes both protein cavities and ligand molecules into 3D voxel grids, then uses dynamic programming to systematically explore all possible molecular placements:

```python
def spatial_alignment(cavity_grid, molecule_grid, compatibility_matrix):
    """
    Find optimal molecular alignment using 3D dynamic programming.
    
    Returns:
        score: Optimal compatibility score
        translation: Best (x, y, z) offset for molecule placement
    """
    best_score = -∞
    best_translation = (0, 0, 0)
    
    for translation in enumerate_translations(cavity_grid):
        score = calculate_placement_score(cavity_grid, molecule_grid, 
                                        translation, compatibility_matrix)
        if score > best_score:
            best_score = score
            best_translation = translation
    
    return best_score, best_translation
```

## Key Features

### **Deterministic Core**
- **3D voxel representation** of cavities and molecules
- **Chemical compatibility matrices** (analogous to BLOSUM for sequences)
- **Dynamic programming** ensures optimal solutions
- **Gap penalties** for cavity filling and molecular clashes

### **Advanced Extensions**
- **Boltzmann ensemble docking**: Thermodynamically weighted pose distributions
- **Pareto frontier analysis**: Multi-objective optimization (shape, electrostatics, entropy)
- **Conformational sampling**: Flexible ligand docking with torsional freedom
- **Adaptive scoring**: Machine learning-based parameter optimization

### **Interpretability**
- Every score component traceable to specific voxel interactions
- Chemical reasoning preserved throughout alignment process
- Visual debugging of molecular compatibility decisions

## Examples

### Basic Docking
```python
# Simple protein-ligand docking
score, pose = docker.align(cavity, molecule)
```

### Ensemble Analysis
```python
# Boltzmann-weighted ensemble of poses
from spatial_docking.boltzmann import ensemble_docking

probabilities = ensemble_docking(cavity, molecule, temperature=300)
expected_binding_energy = calculate_ensemble_average(probabilities)
```

### Multi-Objective Optimization
```python
# Pareto frontier analysis
from spatial_docking.pareto import pareto_optimal_docking

objectives = {
    'shape_fit': calculate_shape_complementarity,
    'electrostatics': calculate_electrostatic_energy,
    'entropy': calculate_conformational_entropy
}

pareto_poses = pareto_optimal_docking(cavity, molecule, objectives)
```

### Flexible Docking
```python
# Include rotatable bond sampling
from spatial_docking.conformers import flexible_docking

best_score, best_pose, best_conformation = flexible_docking(
    cavity, molecule, rotatable_bonds=['C-C', 'C-N', 'C-O']
)
```

## Contributing

We welcome contributions!

**Areas of interest:**
- GPU acceleration for larger molecular systems
- Integration with existing docking pipelines
- Novel scoring function development
- Protein flexibility incorporation

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Inspired by the foundational work of Needleman & Wunsch on sequence alignment
- Built using NumPy, SciPy, and Matplotlib

## Roadmap

- **v1.1**: GPU acceleration via CUDA/OpenCL
- **v1.2**: Protein side-chain flexibility
- **v1.3**: Machine learning-optimized scoring functions
- **v2.0**: Integration with molecular dynamics workflows

---
