"""
Create the geometry diagram showing the cylindrical pipe setup
This recreates the figure from the assignment with all labeled parameters
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Wedge
from matplotlib.patches import Arc

# Create figure
fig, ax = plt.subplots(figsize=(14, 6))

# Pipe dimensions (in arbitrary units for visualization)
L = 10.0   # Length
D = 2.0    # Diameter (height)
z_start = 1.0

# Draw the pipe walls
pipe_top = patches.Rectangle((z_start, D/2), L, 0.15,
                             linewidth=1.5, edgecolor='black',
                             facecolor='gray', hatch='////')
pipe_bottom = patches.Rectangle((z_start, -D/2-0.15), L, 0.15,
                                linewidth=1.5, edgecolor='black',
                                facecolor='gray', hatch='////')
ax.add_patch(pipe_top)
ax.add_patch(pipe_bottom)

# Draw inlet plane
inlet_line = patches.Rectangle((z_start-0.05, -D/2), 0.05, D,
                               linewidth=1.5, edgecolor='black',
                               facecolor='lightblue', alpha=0.3)
ax.add_patch(inlet_line)

# Draw outlet plane
outlet_line = patches.Rectangle((z_start+L, -D/2), 0.05, D,
                                linewidth=1.5, edgecolor='black',
                                facecolor='white', alpha=0.3)
ax.add_patch(outlet_line)

# Draw coordinate system
# r-axis (vertical) - from center
ax.arrow(-0.3, 0, 0, D/2+0.3, head_width=0.15, head_length=0.1,
         fc='black', ec='black', linewidth=2)
ax.text(-0.5, D/2+0.5, r'$\mathbf{r}$', fontsize=18, fontweight='bold')

# z-axis (horizontal) - at centerline (r=0)
ax.arrow(z_start-0.5, 0, 1.5, 0, head_width=0.1, head_length=0.15,
         fc='black', ec='black', linewidth=2)
ax.text(z_start+1.2, -0.25, r'$\mathbf{z}$', fontsize=18, fontweight='bold')

# Draw velocity profiles at inlet (uniform)
n_arrows_inlet = 7
r_positions = np.linspace(-D/2 + 0.2, D/2 - 0.2, n_arrows_inlet)
for r in r_positions:
    ax.arrow(z_start - 0.7, r, 0.5, 0, head_width=0.12, head_length=0.08,
             fc='blue', ec='blue', linewidth=1.5, alpha=0.7)

# Label inlet velocity
ax.text(z_start - 1.5, 0, r'$W_{\mathrm{in}}$', fontsize=16,
        color='blue', fontweight='bold', ha='center', va='center')

# Label inlet temperature
ax.text(z_start - 0.4, -D/2 - 0.35, r'$T_{\mathrm{in}}$', fontsize=16,
        color='red', fontweight='bold')

# Draw developing velocity profile (parabolic)
z_profile = z_start + L/2
n_points = 50
r_profile = np.linspace(-D/2 + 0.1, D/2 - 0.1, n_points)
# Parabolic profile: w = 2*W_in*(1 - (2r/D)^2)
w_profile = 2.0 * 0.8 * (1 - (2*r_profile/D)**2)

# Draw parabola
profile_z = z_profile + w_profile * 0.0
for i in range(len(r_profile)-1):
    ax.plot([z_profile, z_profile + w_profile[i]*0.5],
            [r_profile[i], r_profile[i]], 'b-', linewidth=1, alpha=0.5)

# Draw parabolic envelope
r_envelope = np.linspace(-D/2 + 0.1, D/2 - 0.1, 100)
w_envelope = 2.0 * 0.8 * (1 - (2*r_envelope/D)**2)
ax.plot(z_profile + w_envelope * 0.5, r_envelope, 'b--', linewidth=2, alpha=0.8)

# Label developing velocity
ax.text(z_profile, D/2 + 0.6, r'$w = w(z,r)$', fontsize=14,
        ha='center', style='italic')

# Draw outlet velocity profile (fully developed parabolic)
z_outlet = z_start + L
r_outlet = np.linspace(-D/2 + 0.1, D/2 - 0.1, n_points)
w_outlet = 2.0 * 0.8 * (1 - (2*r_outlet/D)**2)

# Draw arrows for outlet
n_arrows_outlet = 9
r_arrows = np.linspace(-D/2 + 0.15, D/2 - 0.15, n_arrows_outlet)
for r in r_arrows:
    w_local = 2.0 * 0.8 * (1 - (2*r/D)**2)
    if w_local > 0.05:
        ax.arrow(z_outlet - 0.3, r, w_local * 0.4, 0,
                head_width=0.08, head_length=0.06,
                fc='blue', ec='blue', linewidth=1.2, alpha=0.7)

# Label outlet velocity
ax.text(z_outlet + 0.8, 0, r'$w$', fontsize=16,
        color='blue', fontweight='bold')

# Draw temperature indicators
ax.text(z_outlet + 0.8, D/2 + 0.35, r'$T$', fontsize=16,
        color='red', fontweight='bold')

# Draw heat flux arrows (top wall)
n_heat = 8
z_heat = np.linspace(z_start + 0.5, z_start + L - 0.5, n_heat)
for z in z_heat:
    ax.arrow(z, D/2 + 0.5, 0, -0.25, head_width=0.15, head_length=0.08,
             fc='red', ec='red', linewidth=1.5, alpha=0.8)

# Label top heat flux
ax.text(z_start + L/2, D/2 + 0.9, r'$q_w = \mathrm{const.}$',
        fontsize=14, color='red', fontweight='bold', ha='center')

# Draw heat flux arrows (bottom wall)
for z in z_heat:
    ax.arrow(z, -D/2 - 0.65, 0, 0.25, head_width=0.15, head_length=0.08,
             fc='red', ec='red', linewidth=1.5, alpha=0.8)

# Label bottom heat flux
ax.text(z_start + L/2, -D/2 - 1.05, r'$q_w = \mathrm{const.}$',
        fontsize=14, color='red', fontweight='bold', ha='center')

# Dimension lines
# Diameter D
ax.plot([z_start + L + 0.6, z_start + L + 0.6], [-D/2, D/2],
        'k-', linewidth=1.5)
ax.plot([z_start + L + 0.5, z_start + L + 0.7], [-D/2, -D/2],
        'k-', linewidth=1.5)
ax.plot([z_start + L + 0.5, z_start + L + 0.7], [D/2, D/2],
        'k-', linewidth=1.5)
ax.text(z_start + L + 1.1, 0, r'$D$', fontsize=16,
        fontweight='bold', va='center')

# Length L
ax.plot([z_start, z_start], [-D/2 - 1.4, -D/2 - 1.6],
        'k-', linewidth=1.5)
ax.plot([z_start + L, z_start + L], [-D/2 - 1.4, -D/2 - 1.6],
        'k-', linewidth=1.5)
ax.annotate('', xy=(z_start + L, -D/2 - 1.5), xytext=(z_start, -D/2 - 1.5),
            arrowprops=dict(arrowstyle='<->', linewidth=1.5, color='black'))
ax.text(z_start + L/2, -D/2 - 1.9, r'$L$', fontsize=16,
        fontweight='bold', ha='center')

# Add title
ax.text(z_start + L/2, D/2 + 1.5,
        'Cylindrical Pipe Flow with Constant Wall Heat Flux',
        fontsize=16, fontweight='bold', ha='center')

# Add parameter box
param_text = (r'$\mathbf{Given\ Parameters:}$' + '\n'
              r'$\bullet\ D$ = Pipe diameter' + '\n'
              r'$\bullet\ L = 50D$ = Pipe length' + '\n'
              r'$\bullet\ W_{\mathrm{in}}$ = Inlet velocity (uniform)' + '\n'
              r'$\bullet\ T_{\mathrm{in}}$ = Inlet temperature (uniform)' + '\n'
              r'$\bullet\ q_w = \mathrm{const.}$ = Wall heat flux' + '\n'
              r'$\bullet\ \rho, \mu, \lambda, c_v = \mathrm{const.}$ = Fluid properties' + '\n'
              r'$\bullet\ Re = \rho W_{\mathrm{in}} D / \mu = 100$' + '\n'
              r'$\bullet\ Pr = \mu c_v / \lambda = 5$')

ax.text(z_start + L + 2.5, 0, param_text, fontsize=11,
        verticalalignment='center', bbox=dict(boxstyle='round',
        facecolor='wheat', alpha=0.8, pad=0.8))

# Set axis properties
ax.set_xlim(-2, z_start + L + 6)
ax.set_ylim(-D/2 - 2.5, D/2 + 2)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.savefig('../Plots_python/geometry_diagram.png', dpi=300, bbox_inches='tight')
print("Geometry diagram saved to ../Plots_python/geometry_diagram.png")
plt.show()
