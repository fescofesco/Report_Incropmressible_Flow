"""
Create a staggered grid structure diagram for the finite volume method
Shows where scalars (p, T) and velocity components (u, w) are stored
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Circle, FancyArrowPatch

# Create figure
fig, ax = plt.subplots(figsize=(14, 10))

# Grid parameters
n_cells_r = 4  # Number of cells to show in radial direction
n_cells_z = 5  # Number of cells to show in axial direction
cell_size = 1.0

# Draw grid cells
for i in range(n_cells_r + 1):
    for j in range(n_cells_z + 1):
        # Vertical lines
        ax.plot([j * cell_size, j * cell_size], [0, n_cells_r * cell_size],
                'k-', linewidth=1.5)
        # Horizontal lines
        ax.plot([0, n_cells_z * cell_size], [i * cell_size, i * cell_size],
                'k-', linewidth=1.5)

# Highlight one control volume
highlight_i = 2
highlight_j = 2
highlight = patches.Rectangle(
    (highlight_j * cell_size, highlight_i * cell_size),
    cell_size, cell_size,
    linewidth=3, edgecolor='red', facecolor='yellow', alpha=0.2
)
ax.add_patch(highlight)

# Add label for highlighted cell
ax.text((highlight_j + 0.5) * cell_size, (highlight_i + 0.5) * cell_size,
        r'$CV_{i,j}$', fontsize=14, ha='center', va='center',
        fontweight='bold', color='red')

# Plot scalar variables (p, T) at cell centers
for i in range(n_cells_r):
    for j in range(n_cells_z):
        x_center = (j + 0.5) * cell_size
        y_center = (i + 0.5) * cell_size

        # Draw circle for scalar location
        circle = Circle((x_center, y_center), 0.08, color='blue',
                       fill=True, zorder=10)
        ax.add_patch(circle)

        # Label only a few cells to avoid clutter
        if i == highlight_i and j == highlight_j:
            ax.text(x_center, y_center - 0.35, r'$p_{i,j}, T_{i,j}$',
                   fontsize=11, ha='center', color='blue', fontweight='bold')

# Plot axial velocity (w) at vertical faces (between cells in z-direction)
for i in range(n_cells_r):
    for j in range(n_cells_z + 1):
        x_face = j * cell_size
        y_center = (i + 0.5) * cell_size

        # Draw square for w location
        square_w = patches.Rectangle((x_face - 0.06, y_center - 0.06),
                                     0.12, 0.12, color='green',
                                     fill=True, zorder=10)
        ax.add_patch(square_w)

        # Label specific locations
        if i == highlight_i and j == highlight_j:
            ax.text(x_face - 0.3, y_center, r'$w_{i,j-\frac{1}{2}}$',
                   fontsize=10, ha='right', va='center', color='green',
                   fontweight='bold')
        elif i == highlight_i and j == highlight_j + 1:
            ax.text(x_face + 0.3, y_center, r'$w_{i,j+\frac{1}{2}}$',
                   fontsize=10, ha='left', va='center', color='green',
                   fontweight='bold')

# Plot radial velocity (u) at horizontal faces (between cells in r-direction)
for i in range(n_cells_r + 1):
    for j in range(n_cells_z):
        x_center = (j + 0.5) * cell_size
        y_face = i * cell_size

        # Draw triangle for u location
        triangle = plt.Polygon([
            [x_center - 0.08, y_face - 0.08],
            [x_center + 0.08, y_face - 0.08],
            [x_center, y_face + 0.08]
        ], color='red', fill=True, zorder=10)
        ax.add_patch(triangle)

        # Label specific locations
        if i == highlight_i and j == highlight_j:
            ax.text(x_center, y_face - 0.3, r'$u_{i-\frac{1}{2},j}$',
                   fontsize=10, ha='center', va='top', color='red',
                   fontweight='bold')
        elif i == highlight_i + 1 and j == highlight_j:
            ax.text(x_center, y_face + 0.3, r'$u_{i+\frac{1}{2},j}$',
                   fontsize=10, ha='center', va='bottom', color='red',
                   fontweight='bold')

# Add coordinate axes labels
ax.text(n_cells_z * cell_size + 0.5, -0.5, r'$z$ (axial)',
        fontsize=14, fontweight='bold')
ax.text(-0.7, n_cells_r * cell_size + 0.3, r'$r$ (radial)',
        fontsize=14, fontweight='bold')

# Add arrows showing coordinate directions
ax.annotate('', xy=(n_cells_z * cell_size + 0.3, 0),
            xytext=(n_cells_z * cell_size - 0.5, 0),
            arrowprops=dict(arrowstyle='->', lw=2, color='black'))
ax.annotate('', xy=(0, n_cells_r * cell_size + 0.3),
            xytext=(0, n_cells_r * cell_size - 0.5),
            arrowprops=dict(arrowstyle='->', lw=2, color='black'))

# Add grid spacing labels
ax.annotate('', xy=((highlight_j + 1) * cell_size, -0.8),
            xytext=(highlight_j * cell_size, -0.8),
            arrowprops=dict(arrowstyle='<->', lw=1.5, color='black'))
ax.text((highlight_j + 0.5) * cell_size, -1.1, r'$\Delta z$',
        fontsize=12, ha='center')

ax.annotate('', xy=(-0.8, (highlight_i + 1) * cell_size),
            xytext=(-0.8, highlight_i * cell_size),
            arrowprops=dict(arrowstyle='<->', lw=1.5, color='black'))
ax.text(-1.3, (highlight_i + 0.5) * cell_size, r'$\Delta r$',
        fontsize=12, ha='center', rotation=90, va='center')

# Add index labels for one row and column
for j in range(n_cells_z):
    ax.text((j + 0.5) * cell_size, -1.8, f'$j={j}$',
            fontsize=10, ha='center', color='gray')
for i in range(n_cells_r):
    ax.text(-1.8, (i + 0.5) * cell_size, f'$i={i}$',
            fontsize=10, ha='center', va='center', color='gray')

# Add legend
legend_x = n_cells_z * cell_size + 1.5
legend_y_start = n_cells_r * cell_size - 0.5

# Scalars (p, T)
ax.add_patch(Circle((legend_x, legend_y_start), 0.08, color='blue', fill=True))
ax.text(legend_x + 0.3, legend_y_start,
        r'Scalars: $p, T$ (cell centers)', fontsize=11, va='center')

# Axial velocity w
ax.add_patch(patches.Rectangle((legend_x - 0.06, legend_y_start - 0.8 - 0.06),
                               0.12, 0.12, color='green', fill=True))
ax.text(legend_x + 0.3, legend_y_start - 0.8,
        r'Axial velocity: $w$ (vertical faces)', fontsize=11, va='center')

# Radial velocity u
triangle_legend = plt.Polygon([
    [legend_x - 0.08, legend_y_start - 1.6 - 0.08],
    [legend_x + 0.08, legend_y_start - 1.6 - 0.08],
    [legend_x, legend_y_start - 1.6 + 0.08]
], color='red', fill=True)
ax.add_patch(triangle_legend)
ax.text(legend_x + 0.3, legend_y_start - 1.6,
        r'Radial velocity: $u$ (horizontal faces)', fontsize=11, va='center')

# Add control volume description
ax.text(legend_x, legend_y_start - 2.5,
        r'$\mathbf{Control\ Volume:}$', fontsize=12, fontweight='bold',
        ha='left')
ax.text(legend_x, legend_y_start - 2.9,
        r'$\Delta V = 2\pi r \Delta r \Delta z$', fontsize=10, ha='left')

# Add title
ax.text(n_cells_z * cell_size / 2, n_cells_r * cell_size + 1.2,
        'Staggered Grid Structure for Cylindrical Coordinates',
        fontsize=16, fontweight='bold', ha='center')

ax.text(n_cells_z * cell_size / 2, n_cells_r * cell_size + 0.8,
        'Finite Volume Method - Variable Storage Locations',
        fontsize=12, ha='center', style='italic')

# Add note
note_text = (
    r'$\mathbf{Note:}$ Staggered arrangement prevents pressure-velocity coupling issues\n'
    r'Velocity components stored at control volume faces where fluxes are computed'
)
ax.text(n_cells_z * cell_size / 2, -2.5, note_text,
        fontsize=10, ha='center', va='top',
        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

# Set axis properties
ax.set_xlim(-2.5, n_cells_z * cell_size + 5)
ax.set_ylim(-3, n_cells_r * cell_size + 1.5)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.savefig('../Plots_python/grid_structure_diagram.png', dpi=300, bbox_inches='tight')
print("Grid structure diagram saved to ../Plots_python/grid_structure_diagram.png")
plt.show()
