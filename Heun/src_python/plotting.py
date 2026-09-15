"""
Plotting functions for flow simulation results
All functions return matplotlib axes and do NOT perform calculations
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def plot_velocity_contour_axial(w_star, r_grid, z_grid, title="Axial Velocity w/W_in"):
    """
    Plot contour of non-dimensional axial velocity

    Parameters:
    -----------
    w_star : ndarray (n_r, n_z)
        Non-dimensional axial velocity field
    r_grid : ndarray (n_r,)
        Radial grid coordinates
    z_grid : ndarray (n_z,)
        Axial grid coordinates
    title : str
        Plot title

    Returns:
    --------
    ax : matplotlib axes
    """
    fig, ax = plt.subplots(figsize=(12, 4))

    Z, R = np.meshgrid(z_grid, r_grid)
    contour = ax.contourf(Z, R, w_star, levels=50, cmap='jet')
    ax.contour(Z, R, w_star, levels=10, colors='black', linewidths=0.5, alpha=0.3)

    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('w/W_in', rotation=270, labelpad=20)

    ax.set_xlabel('z/D')
    ax.set_ylabel('r/D')
    ax.set_title(title)
    ax.set_aspect('auto')

    return ax


def plot_velocity_contour_radial(u_star, r_grid, z_grid, title="Radial Velocity u/W_in"):
    """
    Plot contour of non-dimensional radial velocity

    Parameters:
    -----------
    u_star : ndarray (n_r, n_z)
        Non-dimensional radial velocity field
    r_grid : ndarray (n_r,)
        Radial grid coordinates
    z_grid : ndarray (n_z,)
        Axial grid coordinates
    title : str
        Plot title

    Returns:
    --------
    ax : matplotlib axes
    """
    fig, ax = plt.subplots(figsize=(12, 4))

    Z, R = np.meshgrid(z_grid, r_grid)
    contour = ax.contourf(Z, R, u_star, levels=50, cmap='RdBu_r')
    ax.contour(Z, R, u_star, levels=10, colors='black', linewidths=0.5, alpha=0.3)

    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('u/W_in', rotation=270, labelpad=20)

    ax.set_xlabel('z/D')
    ax.set_ylabel('r/D')
    ax.set_title(title)
    ax.set_aspect('auto')

    return ax


def plot_temperature_contour(theta, r_grid, z_grid, title="Non-dimensional Temperature θ"):
    """
    Plot contour of non-dimensional temperature

    Parameters:
    -----------
    theta : ndarray (n_r, n_z)
        Non-dimensional temperature field (T-T_in)·ρ·W_in·c_v/q_w
    r_grid : ndarray (n_r,)
        Radial grid coordinates
    z_grid : ndarray (n_z,)
        Axial grid coordinates
    title : str
        Plot title

    Returns:
    --------
    ax : matplotlib axes
    """
    fig, ax = plt.subplots(figsize=(12, 4))

    Z, R = np.meshgrid(z_grid, r_grid)
    contour = ax.contourf(Z, R, theta, levels=50, cmap='hot')
    ax.contour(Z, R, theta, levels=10, colors='black', linewidths=0.5, alpha=0.3)

    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('(T-T_in)·ρ·W_in·c_v/q_w', rotation=270, labelpad=20)

    ax.set_xlabel('z/D')
    ax.set_ylabel('r/D')
    ax.set_title(title)
    ax.set_aspect('auto')

    return ax


def plot_velocity_profile(w_star, r_grid, z_positions, z_grid, w_analytical=None):
    """
    Plot velocity profiles at selected axial positions

    Parameters:
    -----------
    w_star : ndarray (n_r, n_z)
        Non-dimensional axial velocity field
    r_grid : ndarray (n_r,)
        Radial grid coordinates
    z_positions : list of float
        Axial positions where to extract profiles (in units of D)
    z_grid : ndarray (n_z,)
        Axial grid coordinates
    w_analytical : callable, optional
        Function w_analytical(r) for fully developed profile

    Returns:
    --------
    ax : matplotlib axes
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    colors = cm.viridis(np.linspace(0, 1, len(z_positions)))

    for i, z_pos in enumerate(z_positions):
        # Find closest grid index
        z_idx = np.argmin(np.abs(z_grid - z_pos))
        w_profile = w_star[:, z_idx]
        ax.plot(r_grid, w_profile, 'o-', color=colors[i],
                label=f'z/D = {z_grid[z_idx]:.1f}', markersize=4)

    if w_analytical is not None:
        w_analytical_values = w_analytical(r_grid)
        ax.plot(r_grid, w_analytical_values, 'k--', linewidth=2,
                label='Analytical (fully developed)')

    ax.set_xlabel('r/D')
    ax.set_ylabel('w/W_in')
    ax.set_title('Axial Velocity Profiles at Different z-positions')
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def plot_temperature_profile(theta, r_grid, z_positions, z_grid, theta_analytical=None):
    """
    Plot temperature profiles at selected axial positions

    Parameters:
    -----------
    theta : ndarray (n_r, n_z)
        Non-dimensional temperature field
    r_grid : ndarray (n_r,)
        Radial grid coordinates
    z_positions : list of float
        Axial positions where to extract profiles (in units of D)
    z_grid : ndarray (n_z,)
        Axial grid coordinates
    theta_analytical : callable, optional
        Function theta_analytical(r, z) for the approved fully-developed profile
        (see analytical.py::calculate_analytical_temperature). Only plotted ONCE, at
        the last (largest) z_position -- that formula is the fully-developed
        asymptote and is not a valid prediction at the smaller z_positions (it can
        extrapolate to negative theta there; see report Appendix B).

    Returns:
    --------
    ax : matplotlib axes
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    colors = cm.plasma(np.linspace(0, 1, len(z_positions)))

    for i, z_pos in enumerate(z_positions):
        # Find closest grid index
        z_idx = np.argmin(np.abs(z_grid - z_pos))
        theta_profile = theta[:, z_idx]
        ax.plot(r_grid, theta_profile, 'o-', color=colors[i],
                label=f'z/D = {z_grid[z_idx]:.1f}', markersize=4)

    z_fd_idx = np.argmin(np.abs(z_grid - z_positions[-1]))
    z_fd = z_grid[z_fd_idx]

    if theta_analytical is not None:
        theta_analytical_values = theta_analytical(r_grid, z_fd)
        ax.plot(r_grid, theta_analytical_values, 'k--', linewidth=2,
                label=f'z/D = {z_fd:.1f}')

    ax.set_xlabel('r/D')
    ax.set_ylabel('(T-T_in)·ρ·W_in·c_v/q_w')
    ax.set_title('Temperature Profiles at Different z-positions')
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def plot_convergence_history(residuals_dict, title="Convergence History"):
    """
    Plot convergence history of residuals

    Parameters:
    -----------
    residuals_dict : dict
        Dictionary with keys like 'continuity', 'momentum', 'energy'
        and values as lists of residual values over iterations
    title : str
        Plot title

    Returns:
    --------
    ax : matplotlib axes
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    for key, residuals in residuals_dict.items():
        ax.semilogy(residuals, label=key, linewidth=2)

    ax.set_xlabel('Iteration / Time Step')
    ax.set_ylabel('Residual')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3, which='both')

    return ax


def plot_comparison_with_analytical(numerical, analytical, r_grid,
                                    ylabel="w/W_in", title="Comparison"):
    """
    Compare numerical and analytical profiles

    Parameters:
    -----------
    numerical : ndarray (n_r,)
        Numerical solution profile
    analytical : ndarray (n_r,)
        Analytical solution profile
    r_grid : ndarray (n_r,)
        Radial grid coordinates
    ylabel : str
        Y-axis label
    title : str
        Plot title

    Returns:
    --------
    ax : matplotlib axes
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8),
                                    gridspec_kw={'height_ratios': [3, 1]})

    # Main comparison plot
    ax1.plot(r_grid, numerical, 'o-', label='Numerical', markersize=6)
    ax1.plot(r_grid, analytical, 'k--', linewidth=2, label='Analytical')
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Error plot
    error = np.abs(numerical - analytical)
    ax2.plot(r_grid, error, 'r-', linewidth=2)
    ax2.set_xlabel('r/D')
    ax2.set_ylabel('|Error|')
    ax2.set_title('Absolute error')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    return ax1


def plot_streamlines(u_star, w_star, r_grid, z_grid, title="Streamlines"):
    """
    Plot streamlines of the flow

    Parameters:
    -----------
    u_star : ndarray (n_r, n_z)
        Non-dimensional radial velocity
    w_star : ndarray (n_r, n_z)
        Non-dimensional axial velocity
    r_grid : ndarray (n_r,)
        Radial grid coordinates
    z_grid : ndarray (n_z,)
        Axial grid coordinates
    title : str
        Plot title

    Returns:
    --------
    ax : matplotlib axes
    """
    fig, ax = plt.subplots(figsize=(12, 4))

    Z, R = np.meshgrid(z_grid, r_grid)

    # Plot streamlines
    ax.streamplot(Z, R, w_star, u_star, density=1.5, color='b', linewidth=1)

    # Add velocity magnitude as background
    velocity_mag = np.sqrt(u_star**2 + w_star**2)
    contour = ax.contourf(Z, R, velocity_mag, levels=20, cmap='YlOrRd', alpha=0.6)

    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Velocity magnitude', rotation=270, labelpad=20)

    ax.set_xlabel('z/D')
    ax.set_ylabel('r/D')
    ax.set_title(title)
    ax.set_aspect('auto')

    return ax
