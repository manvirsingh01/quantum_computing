import numpy as np
import pandas as pd
from scipy.constants import pi
from scipy.optimize import minimize_scalar

# --- Variables from the user's image ---
# Crystal length (m) - Not used in phase-matching calculation
L = 2e-3
# Pump wavelength (m)
lambda_p_meters = 405e-9
# Poling period (m)
poling_period_meters = 10e-6
# Angles (rad) - Not used in this Type-0 calculation
theta_s = np.deg2rad(0)
theta_i = np.deg2rad(0)
# Speed of light (m/s) - Not used in calculation
c = 299792458

# DataFrames from the image - Note: These are not used in the n_z/n_y functions below.
# The functions use their own internal thermo-optic coefficients.
temp_fit_coeff_n1 = pd.DataFrame({
    'z': [9.9587e-6, 9.9228e-6, -8.9603e-6, 4.1010e-6],
    'y': [6.2897e-6, 6.3061e-6, -6.0629e-6, 2.6486e-6]
})
temp_fit_coeff_n2 = pd.DataFrame({
    'z': [-1.1882e-8, 10.459e-8, -9.8136e-8, 3.1481e-8],
    'y': [-0.14445e-8, 2.2244e-8, -3.5770e-8, 1.3470e-8]
})


# --- End of variables from image ---


def calculate_poly_coeffs(wl_um, coeffs):
    """Calculates polynomial expression: n(λ) = a0 + a1/λ + a2/λ² + a3/λ³."""
    a0, a1, a2, a3 = coeffs
    return a0 + a1 / wl_um + a2 / wl_um ** 2 + a3 / wl_um ** 3


def n_y(wl_um, T_C):
    """Calculates the temperature-dependent refractive index for y-polarization in PPKTP."""
    l2 = wl_um ** 2
    n0 = np.sqrt(2.29743 + 0.038146 / (l2 - 0.04355) + 23.9573 / (l2 - 30.134))
    # Thermo-optic coefficients for y-polarization
    n1_coeffs = (0.1513e-4, -0.2152e-4, 0.1982e-4, -0.0635e-4)
    n2_coeffs = (0.1560e-7, -0.4497e-7, 0.5199e-7, -0.2003e-7)
    delta_T = T_C - 25.0
    n1 = calculate_poly_coeffs(wl_um, n1_coeffs)
    n2 = calculate_poly_coeffs(wl_um, n2_coeffs)
    delta_n = n1 * delta_T + n2 * delta_T ** 2
    return n0 + delta_n


def n_z(wl_um, T_C):
    """Calculates the temperature-dependent refractive index for z-polarization in PPKTP."""
    l2 = wl_um ** 2
    A = 4.59423
    B = 0.06206 / (l2 - 0.04763)
    C = 110.80672 / (l2 - 86.12171)
    # The sign before C was corrected from + to - for the correct physical model
    n0 = np.sqrt(A + B - C)
    # Thermo-optic coefficients for z-polarization
    n1_coeffs = (0.4997e-4, -1.325e-4, 1.233e-4, -0.505e-4)
    n2_coeffs = (0.3423e-7, -1.464e-7, 1.656e-7, -0.739e-7)
    delta_T = T_C - 25.0
    n1 = calculate_poly_coeffs(wl_um, n1_coeffs)
    n2 = calculate_poly_coeffs(wl_um, n2_coeffs)
    delta_n = n1 * delta_T + n2 * delta_T ** 2
    return n0 + delta_n


def k_vector(wl_m, n):
    """Calculates wave vector magnitude k = 2πn/λ."""
    return 2 * pi * n / wl_m


def delta_k(T_C, wl_p_m, poling_period_m):
    """Calculates the phase mismatch Δk = k_p - k_s - k_i - K for Type-0 (z->zz)."""
    # Assuming degenerate down-conversion for signal and idler wavelengths
    wl_s_m = 2 * wl_p_m
    wl_i_m = wl_s_m

    # Convert wavelengths from meters to micrometers for refractive index functions
    wl_p_um = wl_p_m * 1e6
    wl_s_um = wl_s_m * 1e6
    wl_i_um = wl_i_m * 1e6

    # Get refractive indices for pump, signal, and idler (all z-polarized)
    n_p = n_y(wl_p_um, T_C)
    n_s = n_z(wl_s_um, T_C)
    n_i = n_z(wl_i_um, T_C)

    # Calculate wave vectors
    k_p = k_vector(wl_p_m, n_p)
    k_s = k_vector(wl_s_m, n_s)
    k_i = k_vector(wl_i_m, n_i)

    # Calculate the grating vector from the poling period
    K = 2 * pi / poling_period_m

    # Return the absolute phase mismatch
    return abs(k_p - k_s - k_i - K)


def find_optimal_temperature(wl_p_m, poling_period_m):
    """Finds the temperature that minimizes phase mismatch Δk."""
    poling_period_um = poling_period_m * 1e6
    print(f"--- Calculating for Poling Period: {poling_period_um:.2f} µm ---")

    # Create a lambda function to pass to the optimizer
    objective_function = lambda T_C: delta_k(T_C, wl_p_m, poling_period_m)

    # Run the bounded scalar optimization
    result = minimize_scalar(objective_function, bounds=(20, 1000000000), method='bounded')

    if result.success:
        print(f"Optimum temperature for phase matching: {result.x:.2f} °C")
        print(f"Phase mismatch Δk at optimum T: {result.fun:.4e} m⁻¹\n")
    else:
        print("Optimization failed: Could not find a phase-matching temperature.\n")


# --- Main execution block ---
if __name__ == "__main__":
    # Run the optimization using the variables defined at the top of the script
    find_optimal_temperature(
        wl_p_m=lambda_p_meters,
        poling_period_m=poling_period_meters
    )