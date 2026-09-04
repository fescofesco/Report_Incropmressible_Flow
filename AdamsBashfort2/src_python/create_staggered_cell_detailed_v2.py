"""
Create detailed sketch of a staggered grid showing neighboring cells
with ARROWS for velocity vectors (not dots)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Circle, FancyArrowPatch

fig, ax = plt.subplots(figsize=(14, 11))

# Cell dimensions
cell_width = 3.0
cell_height = 2.5

# Grid of cells (3x3 to show neighbors)
for i in range(-1, 2):
    for j in range(-1, 2):
        center_x = 7.0 + j * cell_width
        center_y = 6.0 + i * cell_height

        # Highlight central cell (i,j)
        if i == 0 and j == 0:
            cell = patches.Rectangle((center_x - cell_width/2, center_y - cell_height/2),
                                    cell_width, cell_height,
                                    linewidth=3, edgecolor='black',
                                    facecolor='lightblue', alpha=0.4)
            ax.add_patch(cell)
            ax.text(center_x, center_y + cell_height/2 + 0.6,
                   r'$(i,j)$', fontsize=16, ha='center',
                   fontweight='bold',
                   bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))
        else:
            # Neighboring cells
            cell = patches.Rectangle((center_x - cell_width/2, center_y - cell_height/2),
                                    cell_width, cell_height,
                                    linewidth=2, edgecolor='gray',
                                    facecolor='white', alpha=0.3, linestyle='--')
            ax.add_patch(cell)
            # Label neighboring cells
            if i == 0 and j == -1:
                ax.text(center_x, center_y + cell_height/2 + 0.4,
                       r'$(i,j-1)$', fontsize=12, ha='center', color='gray')
            elif i == 0 and j == 1:
                ax.text(center_x, center_y + cell_height/2 + 0.4,
                       r'$(i,j+1)$', fontsize=12, ha='center', color='gray')
            elif i == -1 and j == 0:
                ax.text(center_x, center_y + cell_height/2 + 0.4,
                       r'$(i-1,j)$', fontsize=12, ha='center', color='gray')
            elif i == 1 and j == 0:
                ax.text(center_x, center_y + cell_height/2 + 0.4,
                       r'$(i+1,j)$', fontsize=12, ha='center', color='gray')

# Central cell coordinates
center_x = 7.0
center_y = 6.0

# Add pressure and temperature at cell centers for all cells
for i in range(-1, 2):
    for j in range(-1, 2):
        cx = 7.0 + j * cell_width
        cy = 6.0 + i * cell_height

        # Pressure (blue circle)
        p_size = 0.18 if (i==0 and j==0) else 0.12
        p_circle = Circle((cx, cy), p_size, color='blue', zorder=10,
                         alpha=1.0 if (i==0 and j==0) else 0.5)
        ax.add_patch(p_circle)

        if i == 0 and j == 0:
            ax.text(cx - 0.5, cy, r'$p_{i,j}$', fontsize=14,
                   ha='right', fontweight='bold', color='blue')

# Add ARROWS for radial velocities u at horizontal faces (NORMAL to faces = vertical arrows)
arrow_scale = 0.7
arrow_width = 0.15

# u at bottom face of central cell (i-1/2, j) - pointing UP (positive r direction)
u_bottom_x = center_x
u_bottom_y = center_y - cell_height/2
arrow_u_bottom = FancyArrowPatch((u_bottom_x, u_bottom_y - arrow_scale/3),
                                (u_bottom_x, u_bottom_y + arrow_scale/2),
                                arrowstyle='->', mutation_scale=25, linewidth=3,
                                color='red', zorder=15)
ax.add_artist(arrow_u_bottom)
ax.text(u_bottom_x - 1.2, u_bottom_y, r'$u_{i-\frac{1}{2},j}$',
       fontsize=13, ha='right', va='center', fontweight='bold', color='red')

# u at top face of central cell (i+1/2, j) - pointing UP (positive r direction)
u_top_x = center_x
u_top_y = center_y + cell_height/2
arrow_u_top = FancyArrowPatch((u_top_x, u_top_y - arrow_scale/3),
                             (u_top_x, u_top_y + arrow_scale/2),
                             arrowstyle='->', mutation_scale=25, linewidth=3,
                             color='red', zorder=15)
ax.add_artist(arrow_u_top)
ax.text(u_top_x + 1.2, u_top_y, r'$u_{i+\frac{1}{2},j}$',
       fontsize=13, ha='left', va='center', fontweight='bold', color='red')

# Add ARROWS for axial velocities w at vertical faces (NORMAL to faces = horizontal arrows)
# w at left face of central cell (i, j-1/2) - pointing RIGHT (positive z direction)
w_left_x = center_x - cell_width/2
w_left_y = center_y
arrow_w_left = FancyArrowPatch((w_left_x - arrow_scale/3, w_left_y),
                              (w_left_x + arrow_scale/2, w_left_y),
                              arrowstyle='->', mutation_scale=25, linewidth=3,
                              color='green', zorder=15)
ax.add_artist(arrow_w_left)
ax.text(w_left_x, w_left_y - 0.9, r'$w_{i,j-\frac{1}{2}}$',
       fontsize=13, ha='center', fontweight='bold', color='green')

# w at right face of central cell (i, j+1/2) - pointing RIGHT (positive z direction)
w_right_x = center_x + cell_width/2
w_right_y = center_y
arrow_w_right = FancyArrowPatch((w_right_x - arrow_scale/3, w_right_y),
                               (w_right_x + arrow_scale/2, w_right_y),
                               arrowstyle='->', mutation_scale=25, linewidth=3,
                               color='green', zorder=15)
ax.add_artist(arrow_w_right)
ax.text(w_right_x, w_right_y + 0.9, r'$w_{i,j+\frac{1}{2}}$',
       fontsize=13, ha='center', fontweight='bold', color='green')

# Add dimension labels
# Δr (radial direction - vertical)
ax.annotate('', xy=(center_x - cell_width/2 - 1.5, center_y + cell_height/2),
           xytext=(center_x - cell_width/2 - 1.5, center_y - cell_height/2),
           arrowprops=dict(arrowstyle='<->', lw=2.5, color='black'))
ax.text(center_x - cell_width/2 - 2.0, center_y, r'$\Delta r$',
       fontsize=18, ha='center', va='center', rotation=90, fontweight='bold')

# Δz (axial direction - horizontal)
ax.annotate('', xy=(center_x + cell_width/2, center_y - cell_height/2 - 1.2),
           xytext=(center_x - cell_width/2, center_y - cell_height/2 - 1.2),
           arrowprops=dict(arrowstyle='<->', lw=2.5, color='black'))
ax.text(center_x, center_y - cell_height/2 - 1.7, r'$\Delta z$',
       fontsize=18, ha='center', fontweight='bold')

# Add coordinate axes
ax.arrow(1.5, 1.5, 0, 2, head_width=0.2, head_length=0.15,
        fc='black', ec='black', linewidth=2.5)
ax.text(1.2, 3.8, r'$r$', fontsize=18, fontweight='bold')

ax.arrow(1.5, 1.5, 2, 0, head_width=0.15, head_length=0.2,
        fc='black', ec='black', linewidth=2.5)
ax.text(3.8, 1.2, r'$z$', fontsize=18, fontweight='bold')

# Add legend
legend_elements = [
    patches.Patch(facecolor='blue', edgecolor='blue',
                 label=r'Scalars: $p, \theta$ at cell centers'),
    FancyArrowPatch((0, 0), (0.5, 0), arrowstyle='->', mutation_scale=20,
                   linewidth=3, color='red',
                   label=r'Radial velocity: $u$ at horizontal faces'),
    FancyArrowPatch((0, 0), (0, 0.5), arrowstyle='->', mutation_scale=20,
                   linewidth=3, color='green',
                   label=r'Axial velocity: $w$ at vertical faces')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=13,
         framealpha=0.95, edgecolor='black', fancybox=True)

# Add information box
info_text = (r'$\mathbf{Staggered\ Grid\ Benefits:}$' + '\n\n'
            r'$\bullet$ Scalars at cell centers' + '\n'
            r'$\bullet$ Velocities at cell faces' + '\n'
            r'$\bullet$ No pressure-velocity decoupling' + '\n'
            r'$\bullet$ No checkerboard modes' + '\n'
            r'$\bullet$ Natural flux calculations' + '\n'
            r'$\bullet$ Exact conservation' + '\n\n'
            r'$\mathbf{Volume:}\ dV = 2\pi r \Delta r \Delta z$')
ax.text(1.2, 10.0, info_text, fontsize=11,
       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9, pad=0.9))

# Add title
ax.text(7.0, 11.5, 'Staggered Grid: Variable Storage with Neighboring Cells',
       fontsize=18, fontweight='bold', ha='center')
ax.text(7.0, 11.0, 'Arrows show velocity directions at cell faces',
       fontsize=13, ha='center', style='italic', color='navy')

# Set limits and aspect
ax.set_xlim(0, 14)
ax.set_ylim(0, 12)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.savefig('../Plots_python/staggered_cell_detailed.png', dpi=300, bbox_inches='tight')
print("Updated staggered cell diagram with arrows and neighboring cells saved!")
plt.show()
