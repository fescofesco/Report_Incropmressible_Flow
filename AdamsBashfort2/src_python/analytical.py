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

    Profile: θ = 4*(z/D) + Re*Pr*[(1/2)*(2r/D)² - (1/8)*(2r/D)⁴ - 7/48]

    Where θ = (T - T_in)*ρ*W_in*c_v/q_w

    This is the approved fully-developed profile. It is obtained by direct
    integration of the non-dimensional energy equation with the wall condition
    ∂θ/∂r* = Re*Pr; the constant -7/48 follows from the exact global energy balance
    (the flow-weighted/bulk average of the radial bracket must vanish, so
    θ_bulk = 4z* for all z >= 0). It gives Nu = 48/11 (see
    calculate_nusselt_number_fully_developed). An earlier transcription of the
    assignment printed the radial terms without the factor 1/2 and with -3/4 instead
    of -7/48 (that form gives Nu = 24/11 and violates θ_bulk = 4z*); it is not used
    anywhere in this codebase. Full derivation and sources:
    Report/nusselt_number_analysis.md.

    Caveat (inherent, not a bug): the radial *shape* of this profile is only valid
    once the flow is thermally fully developed (z* beyond ~25); the *level*
    θ_bulk = 4z* is correct from the inlet on.

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
    term2 = Re * Pr * (0.5 * r_ratio**2 - 0.125 * r_ratio**4 - 7.0 / 48.0)

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

    Profile: θ = 4*z* + Re*Pr*[(1/2)*(2*r*)² - (1/8)*(2*r*)⁴ - 7/48]
    (see calculate_analytical_temperature for the derivation and
    Report/nusselt_number_analysis.md for the full analysis)

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

    For parabolic profile w* = 2(1-(2r*)²) with mean velocity 1 (matching
    the uniform inlet), Q* = mean_velocity * area = 1 * pi*R*^2:
    Q* = pi * R*^2  (non-dimensional)

    For R* = 0.5 (half diameter), Q* = pi/4

    Cross-checked against numerical integration of the velocity profile
    (calculate_volumetric_flow_rate): matches to 1e-6. The previous formula
    (pi/2 * R*^2 = pi/8 at R*=0.5) was off by a factor of 2.

    Parameters:
    -----------
    R_star : float, optional
        Non-dimensional pipe radius (default 0.5 for D=1)

    Returns:
    --------
    Q_star : float
        Non-dimensional flow rate
    """
    Q_star = np.pi * R_star**2

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
