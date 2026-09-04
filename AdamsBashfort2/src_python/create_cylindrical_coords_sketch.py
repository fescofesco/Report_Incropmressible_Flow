"""
Create a sketch showing cylindrical coordinate system with velocity components
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform

class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)

    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)

fig = plt.figure(figsize=(14, 10))
ax = fig.add_subplot(111, projection='3d')

# Draw cylindrical coordinate axes
origin = np.array([0, 0, 0])

# r-axis (radial)
r_arrow = Arrow3D(0, 0, 0, 2, 0, 0,
                  mutation_scale=20, lw=3, arrowstyle='-|>', color='red')
ax.add_artist(r_arrow)
ax.text(2.3, 0, 0, r'$r$', fontsize=20, fontweight='bold', color='red')

# φ-axis (azimuthal)
phi_points = np.linspace(0, np.pi/3, 20)
x_phi = 1.5 * np.cos(phi_points)
y_phi = 1.5 * np.sin(phi_points)
z_phi = np.zeros_like(phi_points)
ax.plot(x_phi, y_phi, z_phi, 'g-', linewidth=3)
ax.text(1.0, 1.3, 0, r'$\phi$', fontsize=20, fontweight='bold', color='green')

# z-axis (axial)
z_arrow = Arrow3D(0, 0, 0, 0, 0, 2.5,
                  mutation_scale=20, lw=3, arrowstyle='-|>', color='blue')
ax.add_artist(z_arrow)
ax.text(0, 0, 2.7, r'$z$', fontsize=20, fontweight='bold', color='blue')

# Draw a point P at (r, φ, z)
r_p = 1.5
phi_p = np.pi/6
z_p = 1.5

x_p = r_p * np.cos(phi_p)
y_p = r_p * np.sin(phi_p)

ax.scatter([x_p], [y_p], [z_p], color='black', s=100, zorder=10)
ax.text(x_p + 0.2, y_p + 0.2, z_p, 'P', fontsize=16, fontweight='bold')

# Draw velocity components at point P
scale = 0.8

# u (radial velocity)
u_dir = np.array([np.cos(phi_p), np.sin(phi_p), 0])
u_arrow = Arrow3D(x_p, y_p, z_p, scale*u_dir[0], scale*u_dir[1], scale*u_dir[2],
                  mutation_scale=20, lw=2.5, arrowstyle='-|>', color='red')
ax.add_artist(u_arrow)
ax.text(x_p + scale*u_dir[0] + 0.3, y_p + scale*u_dir[1], z_p,
        r'$u$ (radial)', fontsize=14, color='red', fontweight='bold')

# v (azimuthal velocity)
v_dir = np.array([-np.sin(phi_p), np.cos(phi_p), 0])
v_arrow = Arrow3D(x_p, y_p, z_p, scale*v_dir[0], scale*v_dir[1], scale*v_dir[2],
                  mutation_scale=20, lw=2.5, arrowstyle='-|>', color='green')
ax.add_artist(v_arrow)
ax.text(x_p + scale*v_dir[0] - 0.3, y_p + scale*v_dir[1] + 0.3, z_p,
        r'$v$ (azimuthal)', fontsize=14, color='green', fontweight='bold')

# w (axial velocity)
w_arrow = Arrow3D(x_p, y_p, z_p, 0, 0, scale,
                  mutation_scale=20, lw=2.5, arrowstyle='-|>', color='blue')
ax.add_artist(w_arrow)
ax.text(x_p, y_p + 0.3, z_p + scale + 0.2,
        r'$w$ (axial)', fontsize=14, color='blue', fontweight='bold')

# Draw projection lines
ax.plot([0, x_p], [0, y_p], [0, 0], 'k--', alpha=0.3, linewidth=1)
ax.plot([x_p, x_p], [y_p, y_p], [0, z_p], 'k--', alpha=0.3, linewidth=1)

# Draw r dimension
ax.plot([0, x_p], [0, y_p], [0, 0], 'r-', linewidth=2, alpha=0.5)
ax.text(x_p/2, y_p/2, -0.3, r'$r$', fontsize=12, color='red')

# Draw z dimension
ax.plot([0, 0], [0, 0], [0, z_p], 'b-', linewidth=2, alpha=0.5)
ax.text(-0.3, 0, z_p/2, r'$z$', fontsize=12, color='blue')

# Draw circular arc for phi
phi_arc = np.linspace(0, phi_p, 20)
x_arc = 0.5 * np.cos(phi_arc)
y_arc = 0.5 * np.sin(phi_arc)
z_arc = np.zeros_like(phi_arc)
ax.plot(x_arc, y_arc, z_arc, 'g-', linewidth=2, alpha=0.7)

# Add title and labels
ax.set_title('Cylindrical Coordinate System\nVelocity Components: $\\mathbf{u} = u\\,\\mathbf{e}_r + v\\,\\mathbf{e}_\\phi + w\\,\\mathbf{e}_z$',
             fontsize=16, fontweight='bold', pad=20)

# Add text box with explanation
textstr = ('$u$ = radial velocity (r-direction)\n'
           '$v$ = azimuthal velocity ($\\phi$-direction)\n'
           '$w$ = axial velocity (z-direction)\n\n'
           'For axisymmetric flow: $v = 0$, $\\partial/\\partial\\phi = 0$')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax.text2D(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=12,
          verticalalignment='top', bbox=props)

# Set axis limits
ax.set_xlim([-0.5, 2.5])
ax.set_ylim([-0.5, 2])
ax.set_zlim([0, 3])

ax.set_xlabel('X', fontsize=12)
ax.set_ylabel('Y', fontsize=12)
ax.set_zlabel('Z', fontsize=12)

# Set viewing angle
ax.view_init(elev=20, azim=45)

plt.tight_layout()
plt.savefig('../Plots_python/cylindrical_coordinates_sketch.png', dpi=300, bbox_inches='tight')
print("Cylindrical coordinates sketch saved!")
plt.show()
