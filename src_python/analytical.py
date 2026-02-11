"""
Analytical solutions for fully developed flow in a cylindrical pipe
These are calculation functions that return results (no plotting)
"""

import numpy as np


def calculate_analytical_velocity(r_grid, W_in=1.0, D=1.0):
    """
    Calculate analytical velocity profile for fully developed laminar pipe flow

    Profile: w = 2*W_in*(1 - (2r/D)²)

    Parameters:
    -----------
    r_grid : ndarray
        Radial coordinates (can be dimensional or non-dimensional)
    W_in : float, optional
        Inlet velocity (default 1.0 for non-dimensional)
    D : float, optional
        Pipe diameter (default 1.0 for non-dimensional)

    Returns:
    --------
    w_analytical : ndarray
        Analytical axial velocity profile
    """
    r_ratio = 2.0 * r_grid / D
    w_analytical = 2.0 * W_in * (1.0 - r_ratio**2)

    return w_analytical


def calculate_analytical_temperature(r_grid, z, Re, Pr, D=1.0):
    """
    Calculate analytical temperature profile for fully developed thermal flow

    Profile: θ = 4*(z/D) + Re*Pr*[(2r/D)² - (1/4)*(2r/D)⁴ - 3/4]

    Where θ = (T - T_in)*ρ*W_in*c_v/q_w

    Parameters:
    -----------
    r_grid : ndarray
        Radial coordinates (non-dimensional)
    z : float
        Axial position (non-dimensional)
    Re : float
        Reynolds number
    Pr : float
        Prandtl number
    D : float, optional
        Pipe diameter (default 1.0 for non-dimensional)

    Returns:
    --------
    theta_analytical : ndarray
        Analytical non-dimensional temperature profile
    """
    r_ratio = 2.0 * r_grid / D

    # Linear term in z
    term1 = 4.0 * z / D

    # Radial profile term
    term2 = Re * Pr * (r_ratio**2 - 0.25 * r_ratio**4 - 0.75)

    theta_analytical = term1 + term2

    return theta_analytical


def calculate_analytical_velocity_nondim(r_star):
    """
    Calculate analytical velocity profile (non-dimensional form)

    Profile: w* = 2*(1 - (2*r*)²)  where r* = r/D

    Parameters:
    -----------
    r_star : ndarray
        Non-dimensional radial coordinates (r/D)

    Returns:
    --------
    w_star_analytical : ndarray
        Analytical non-dimensional axial velocity
    """
    return calculate_analytical_velocity(r_star, W_in=1.0, D=1.0)


def calculate_analytical_temperature_nondim(r_star, z_star, Re, Pr):
    """
    Calculate analytical temperature profile (non-dimensional form)

    Profile: θ = 4*z* + Re*Pr*[(2*r*)² - (1/4)*(2*r*)⁴ - 3/4]

    Parameters:
    -----------
    r_star : ndarray
        Non-dimensional radial coordinates (r/D)
    z_star : float
        Non-dimensional axial position (z/D)
    Re : float
        Reynolds number
    Pr : float
        Prandtl number

    Returns:
    --------
    theta_analytical : ndarray
        Analytical non-dimensional temperature
    """
    return calculate_analytical_temperature(r_star, z_star, Re, Pr, D=1.0)


def calculate_centerline_velocity():
    """
    Calculate centerline velocity for fully developed flow

    Returns:
    --------
    w_centerline : float
        Non-dimensional centerline velocity (w_centerline/W_in = 2.0)
    """
    return 2.0


def calculate_volumetric_flow_rate(r_grid, w_profile):
    """
    Calculate volumetric flow rate from velocity profile

    Q = ∫∫ w · r dr dφ = 2π ∫ w·r dr

    Parameters:
    -----------
    r_grid : ndarray
        Radial grid coordinates
    w_profile : ndarray
        Axial velocity profile at a given z-position

    Returns:
    --------
    Q : float
        Volumetric flow rate
    """
    # Trapezoidal integration: Q = 2π ∫ w·r dr
    integrand = w_profile * r_grid
    Q = 2.0 * np.pi * np.trapz(integrand, r_grid)

    return Q


def calculate_mean_velocity(r_grid, w_profile, R):
    """
    Calculate area-averaged mean velocity

    W_mean = Q / A = Q / (π R²)

    Parameters:
    -----------
    r_grid : ndarray
        Radial grid coordinates
    w_profile : ndarray
        Axial velocity profile
    R : float
        Pipe radius

    Returns:
    --------
    W_mean : float
        Mean velocity
    """
    Q = calculate_volumetric_flow_rate(r_grid, w_profile)
    A = np.pi * R**2
    W_mean = Q / A

    return W_mean


def calculate_analytical_flow_rate_nondim(R_star=0.5):
    """
    Calculate analytical volumetric flow rate for fully developed flow

    For parabolic profile w* = 2(1-(2r*)²), we have:
    Q* = π/2 * R*²  (non-dimensional)

    For R* = 0.5 (half diameter), Q* = π/8

    Parameters:
    -----------
    R_star : float, optional
        Non-dimensional pipe radius (default 0.5 for D=1)

    Returns:
    --------
    Q_star : float
        Non-dimensional flow rate
    """
    Q_star = 0.5 * np.pi * R_star**2

    return Q_star


def calculate_nusselt_number_fully_developed():
    """
    Calculate Nusselt number for fully developed thermal flow
    with constant wall heat flux (H1 boundary condition)

    For laminar pipe flow: Nu = 48/11 ≈ 4.364

    Returns:
    --------
    Nu : float
        Nusselt number
    """
    Nu = 48.0 / 11.0  # Exact analytical value

    return Nu


def calculate_friction_factor(Re):
    """
    Calculate Darcy friction factor for fully developed laminar pipe flow

    f = 64/Re

    Parameters:
    -----------
    Re : float
        Reynolds number

    Returns:
    --------
    f : float
        Darcy friction factor
    """
    f = 64.0 / Re

    return f


def calculate_entry_length_hydrodynamic(Re, D=1.0):
    """
    Calculate hydrodynamic entry length (approximation)

    L_h ≈ 0.05 * Re * D  (for laminar flow)

    Parameters:
    -----------
    Re : float
        Reynolds number
    D : float, optional
        Pipe diameter (default 1.0 for non-dimensional)

    Returns:
    --------
    L_h : float
        Hydrodynamic entry length
    """
    L_h = 0.05 * Re * D

    return L_h


def calculate_entry_length_thermal(Re, Pr, D=1.0):
    """
    Calculate thermal entry length (approximation)

    L_t ≈ 0.05 * Re * Pr * D  (for laminar flow)

    Parameters:
    -----------
    Re : float
        Reynolds number
    Pr : float
        Prandtl number
    D : float, optional
        Pipe diameter (default 1.0 for non-dimensional)

    Returns:
    --------
    L_t : float
        Thermal entry length
    """
    L_t = 0.05 * Re * Pr * D

    return L_t
