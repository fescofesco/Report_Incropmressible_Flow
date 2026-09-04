"""
Create detailed sketch of a single staggered grid cell
showing where u_ij, w_ij, p_ij are located
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(12, 10))

# Cell dimensions
cell_width = 4.0
cell_height = 3.0
center_x = 5.0
center_y = 4.0

# Draw the central cell (i,j)
cell = patches.Rectangle((center_x - cell_width/2, center_y - cell_height/2),
                         cell_width, cell_height,
                         linewidth=3, edgecolor='black',
                         facecolor='lightblue', alpha=0.3)
ax.add_patch(cell)

# Draw neighboring cells (faint)
# Left cell (i-1, j)
left_cell = patches.Rectangle((center_x - 1.5*cell_width, center_y - cell_height/2),
                              cell_width, cell_height,
                              linewidth=2, edgecolor='gray',
                              facecolor='white', alpha=0.2, linestyle='--')
ax.add_patch(left_cell)

# Right cell (i+1, j)
right_cell = patches.Rectangle((center_x + cell_width/2, center_y - cell_height/2),
                               cell_width, cell_height,
                               linewidth=2, edgecolor='gray',
                               facecolor='white', alpha=0.2, linestyle='--')
ax.add_patch(right_cell)

# Bottom cell (i, j-1)
bottom_cell = patches.Rectangle((center_x - cell_width/2, center_y - 1.5*cell_height),
                                cell_width, cell_height,
                                linewidth=2, edgecolor='gray',
                                facecolor='white', alpha=0.2, linestyle='--')
ax.add_patch(bottom_cell)

# Top cell (i, j+1)
top_cell = patches.Rectangle((center_x - cell_width/2, center_y + cell_height/2),
                             cell_width, cell_height,
                             linewidth=2, edgecolor='gray',
                             facecolor='white', alpha=0.2, linestyle='--')
ax.add_patch(top_cell)

# Add pressure p_ij at cell center
p_circle = Circle((center_x, center_y), 0.25, color='blue', zorder=10)
ax.add_patch(p_circle)
ax.text(center_x, center_y - 0.7, r'$p_{i,j}$', fontsize=16,
        ha='center', fontweight='bold', color='blue')

# Add temperature θ_ij at cell center (slightly offset)
theta_circle = Circle((center_x + 0.4, center_y + 0.4), 0.15, color='purple', zorder=10)
ax.add_patch(theta_circle)
ax.text(center_x + 0.4, center_y + 1.0, r'$\theta_{i,j}$', fontsize=14,
        ha='center', fontweight='bold', color='purple')

# Add radial velocities u at horizontal faces
# u at bottom face (i-1/2, j)
u_bottom = patches.Rectangle((center_x - 0.15, center_y - cell_height/2 - 0.15),
                             0.3, 0.3, color='red', zorder=10)
ax.add_patch(u_bottom)
ax.text(center_x - 1.2, center_y - cell_height/2, r'$u_{i-\frac{1}{2},j}$',
        fontsize=14, ha='right', va='center', fontweight='bold', color='red')

# u at top face (i+1/2, j)
u_top = patches.Rectangle((center_x - 0.15, center_y + cell_height/2 - 0.15),
                          0.3, 0.3, color='red', zorder=10)
ax.add_patch(u_top)
ax.text(center_x + 1.2, center_y + cell_height/2, r'$u_{i+\frac{1}{2},j}$',
        fontsize=14, ha='left', va='center', fontweight='bold', color='red')

# Add axial velocities w at vertical faces
# w at left face (i, j-1/2)
w_left = patches.Rectangle((center_x - cell_width/2 - 0.15, center_y - 0.15),
                           0.3, 0.3, color='green', zorder=10)
ax.add_patch(w_left)
ax.text(center_x - cell_width/2, center_y - 0.8, r'$w_{i,j-\frac{1}{2}}$',
        fontsize=14, ha='center', fontweight='bold', color='green')

# w at right face (i, j+1/2)
w_right = patches.Rectangle((center_x + cell_width/2 - 0.15, center_y - 0.15),
                            0.3, 0.3, color='green', zorder=10)
ax.add_patch(w_right)
ax.text(center_x + cell_width/2, center_y + 0.8, r'$w_{i,j+\frac{1}{2}}$',
        fontsize=14, ha='center', fontweight='bold', color='green')

# Add dimension labels
# Δr (radial direction - vertical)
ax.annotate('', xy=(center_x - cell_width/2 - 1.0, center_y + cell_height/2),
            xytext=(center_x - cell_width/2 - 1.0, center_y - cell_height/2),
            arrowprops=dict(arrowstyle='<->', lw=2, color='black'))
ax.text(center_x - cell_width/2 - 1.5, center_y, r'$\Delta r$',
        fontsize=16, ha='center', va='center', rotation=90, fontweight='bold')

# Δz (axial direction - horizontal)
ax.annotate('', xy=(center_x + cell_width/2, center_y - cell_height/2 - 0.8),
            xytext=(center_x - cell_width/2, center_y - cell_height/2 - 0.8),
            arrowprops=dict(arrowstyle='<->', lw=2, color='black'))
ax.text(center_x, center_y - cell_height/2 - 1.3, r'$\Delta z$',
        fontsize=16, ha='center', fontweight='bold')

# Add coordinate axes
ax.arrow(1, 1, 0, 1.5, head_width=0.15, head_length=0.1,
         fc='black', ec='black', linewidth=2)
ax.text(0.7, 2.8, r'$r$', fontsize=16, fontweight='bold')

ax.arrow(1, 1, 1.5, 0, head_width=0.1, head_length=0.15,
         fc='black', ec='black', linewidth=2)
ax.text(2.7, 0.7, r'$z$', fontsize=16, fontweight='bold')

# Add cell label
ax.text(center_x, center_y + cell_height/2 + 0.8,
        r'Cell $(i,j)$', fontsize=18, ha='center', fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

# Add legend
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', label=r'Scalars: $p, \theta$ (cell center)',
               markerfacecolor='blue', markersize=12),
    plt.Line2D([0], [0], marker='s', color='w', label=r'Radial velocity: $u$ (horizontal faces)',
               markerfacecolor='red', markersize=12),
    plt.Line2D([0], [0], marker='s', color='w', label=r'Axial velocity: $w$ (vertical faces)',
               markerfacecolor='green', markersize=12)
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=12,
          framealpha=0.9, edgecolor='black', fancybox=True)

# Add volume element formula
volume_text = (r'$\mathbf{Volume\ Element:}$' + '\n'
               r'$dV = 2\pi r \Delta r \Delta z$' + '\n\n'
               r'$\mathbf{Staggered\ Grid\ Benefits:}$' + '\n'
               r'$\bullet$ No pressure-velocity decoupling' + '\n'
               r'$\bullet$ No checkerboard modes' + '\n'
               r'$\bullet$ Natural flux calculations at faces')
ax.text(0.5, 6.5, volume_text, fontsize=11,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8, pad=0.8))

# Add title
ax.text(center_x, 9.5, 'Staggered Grid: Variable Storage Locations',
        fontsize=18, fontweight='bold', ha='center')

ax.text(center_x, 9.0, 'Finite Volume Method - Cylindrical Coordinates',
        fontsize=14, ha='center', style='italic')

# Set limits and aspect
ax.set_xlim(0, 11)
ax.set_ylim(0, 10)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.savefig('../Plots_python/staggered_cell_detailed.png', dpi=300, bbox_inches='tight')
print("Detailed staggered cell sketch saved!")
plt.show()
