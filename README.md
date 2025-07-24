# Spatial Needleman-Wunsch

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**A deterministic dynamic programming framework for 3D molecular docking**

## 🔬 Overview

Spatial Needleman-Wunsch adapts classical sequence alignment principles to three-dimensional molecular space, providing:

- **Perfect Reproducibility**: Identical input → identical output, every time
- **Mathematical Optimality**: Guaranteed best solution within scoring framework  
- **Complete Interpretability**: Full score decomposition and interaction analysis
- **Multi-Objective Optimization**: Pareto frontier analysis for binding trade-offs

Spatial Needleman-Wunsch solves these problems by:
- ✅ **Guaranteed optimal alignment** given the scoring function
- ✅ **Reproducible results** every time
- ✅ **Complete interpretability** of every score component
- ✅ **Chemical intuition preservation** through explicit compatibility matrices

## Try It Live: Spatial Needleman–Wunsch Notebook

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JDCurry/spatial-needleman-wunsch/HEAD?filepath=Spatial_Needleman_Wunsch_Demo.ipynb)
[![Notebook](https://img.shields.io/badge/Open%20Notebook-ipynb-blue)](./Spatial_Needleman_Wunsch_Demo.ipynb)

Run the full docking demo in your browser, no setup required.

## 🚀 Quick Start

```python
from spatial_docking import SpatialDocking, create_synthetic_cavity, create_test_molecule

# Initialize docking engine
docker = SpatialDocking(grid_spacing=0.5, verbose=True)

# Create molecular systems
cavity = create_synthetic_cavity(size=(8, 8, 6), pattern='checkerboard')
molecule = create_test_molecule(length=5)

# Perform docking
result = docker.align(cavity, molecule)
print(f"Best score: {result.score:.3f} at {result.translation}")

# Analyze interactions
analysis = docker.get_binding_site_analysis(result)
for interaction in analysis['interactions']:
    print(f"  {interaction['molecule_type']} ↔ {interaction['cavity_type']}: {interaction['score']:+.1f}")

# Visualize result
docker.visualize(result)
```

## 📦 Installation

### From Source
```bash
git clone https://github.com/JDCurry/spatial-needleman-wunsch.git
cd spatial-needleman-wunsch
pip install -e .
```

### With Optional Dependencies
```bash
# GPU acceleration
pip install -e ".[gpu]"

# RDKit for real molecules  
pip install -e ".[rdkit]"

# Development tools
pip install -e ".[dev]"
```

### Requirements
- Python 3.8+
- NumPy, SciPy, Matplotlib
- Optional: RDKit, CuPy, Jupyter

## 🔧 Core Features

### Deterministic Docking
```python
# Perfect reproducibility - no random seeds needed
result1 = docker.align(cavity, molecule)
result2 = docker.align(cavity, molecule)
assert result1.score == result2.score  # Always true
```

### Multi-Objective Optimization
```python
from spatial_docking.pareto import multi_objective_spatial_alignment

all_results, pareto_solutions = multi_objective_spatial_alignment(
    cavity, molecule, max_translation=3
)
print(f"Found {len(pareto_solutions)} Pareto-optimal solutions")
```

### Flexible Molecular Conformations
```python
from spatial_docking.conformers import generate_torsion_states

# Generate conformational states
conformations = generate_torsion_states(n_bonds=3, torsion_steps=6)
for conf in conformations:
    result = docker.align(cavity, apply_torsion(molecule, conf))
```

### Boltzmann Ensemble Analysis
```python
from spatial_docking.boltzmann import boltzmann_ensemble_docking

# Thermodynamic weighting of binding modes
probabilities = boltzmann_ensemble_docking(
    spatial_alignment_all_paths, cavity, molecule, temperature=300
)
```

## 📊 Example Results

### Chemical Discrimination Test
```
Ligand Type          Score    Result
─────────────────────────────────────
All Hydrophobic      2.000    ✅ Optimal
Mixed Properties     2.000    ✅ Optimal  
All Polar           0.000    ❌ Poor fit
```

### Conformational Sensitivity
```
Conformation    Global Score    Local Score    Advantage
──────────────────────────────────────────────────────
Linear          2.000          1.000          +1.000
Bent            2.000          1.000          +1.000
Compact         0.000          0.000          +0.000
Helical         0.000          0.000          +0.000
```

### Score Decomposition Example
```
Total Score: 17.5
├── Hydrophobic contacts: +12.0 (6 interactions)
├── Polar contacts: +4.5 (3 interactions)  
├── Electrostatic attraction: +3.0 (1 interaction)
├── Unfavorable overlaps: -1.5 (minor penalties)
└── Gap penalties: -0.5
```

## 🧪 Testing

Run the comprehensive test suite:
```bash
# Basic tests
python -m pytest tests/

# With coverage
python -m pytest tests/ --cov=spatial_docking --cov-report=html

# Performance tests
python -m pytest tests/test_performance.py -v
```

Test coverage: 95%+ across all modules

## 📁 Project Structure

```
spatial-needleman-wunsch/
├── spatial_docking/
│   ├── __init__.py          # Main package interface
│   ├── core.py              # SpatialDocking class & Voxel dataclass
│   ├── alignment.py         # Dynamic programming algorithms
│   ├── scoring.py           # Compatibility matrices
│   ├── visualization.py     # 3D plotting tools
│   ├── boltzmann.py         # Ensemble analysis
│   ├── pareto.py           # Multi-objective optimization
│   ├── conformers.py       # Flexible ligand sampling
│   └── adaptive.py         # ML-based scoring refinement
├── tests/                   # Comprehensive test suite

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
## 🔬 Research Applications

### Drug Discovery
- Lead compound optimization
- Fragment-based design
- Selectivity analysis

### Structural Biology
- Protein-ligand interaction analysis
- Allosteric site characterization  
- Binding mechanism elucidation

### Method Development
- Scoring function benchmarking
- Algorithm comparison studies
- Interpretable AI research
  
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
