# visualization.py
"""3D plotting of cavity and molecular docking results."""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch

def visualize_docking(cavity_grid, molecule_grid, offset):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    color_map = {
        'hydrophobic': 'yellow',
        'hydrophilic': 'blue',
        'polar': 'green',
        'charged_pos': 'red',
        'charged_neg': 'purple',
        'empty': 'gray'
    }

    for pos, voxel in cavity_grid.items():
        if voxel.property_type != 'empty':
            ax.scatter(*pos, color=color_map[voxel.property_type], alpha=0.4, s=100, marker='o')

    for mol_pos, mol_voxel in molecule_grid.items():
        placed_pos = tuple(mol_pos[i] + offset[i] for i in range(3))
        ax.scatter(*placed_pos, color=color_map[mol_voxel.property_type], alpha=0.9, s=200, marker='^')

    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title('Spatial Docking Visualization')

    legend_elements = [Patch(facecolor=color, label=prop)
                       for prop, color in color_map.items() if prop != 'empty']
    ax.legend(handles=legend_elements, loc='upper right')
    plt.show()
