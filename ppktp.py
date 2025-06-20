import numpy as np
import pandas as pd
import scipy
from scipy.constants import pi
from scipy.optimize import minimize_scalar

# --- Print Version Information ---
print(f"SciPy Version: {scipy.__version__}")
print(f"NumPy Version: {np.__version__}\n")

# --- Input Parameters ---
# (The rest of the code is exactly the same as before)
# Crystal length (m)
L = 2e-3
# Pump wavelength (m)
lambda_p_meters = 405e-9
# Poling period (m)
poling_period_meters = 10e-6
# Angles (rad)
theta_s = np.deg2rad(0)
theta_i = np.deg2rad(0)
# Speed of light (m/s)
c = 299792458

# --- Coefficients ---
temp_fit_coeff_n1 = pd.DataFrame({
    'z': [9.9587e-6, 9.9228e-6, -8.9603e-6, 4.1010e-6],
    'y': [6.2897e-6, 6.3061e-6, -6.0629e-6, 2.6486e-6]
})
temp_fit_coeff_n2 = pd.DataFrame({
    'z': [-1.1882e-8, 10.459e-8, -9.8136e-8, 3.1481e-8],
    'y': [-0.14445e-8, 2.2244e-8, -3.5770e-8, 1.3470e-8]
})

def calculate_poly_coeffs(wl_um, coeffs):
    a0, a1, a2, a3 = coeffs
    return a0 + a1 / wl_um + a2 / wl_um ** 2 + a3 / wl_um ** 3

def n_y(wl_um, T_C):
    l2 = wl_um ** 2
    n0_squared = 3.45018 + (0.04341 / (l2 - 0.04597)) + (16.98825 / (l2 - 39.43799))
    n0 = np.sqrt(n0_squared)
    n1_coeffs = temp_fit_coeff_n1['y'].values
    n2_coeffs = temp_fit_coeff_n2['y'].values
    delta_T = T_C - 20.0
    n1 = calculate_poly_coeffs(wl_um, n1_coeffs)
    n2 = calculate_poly_coeffs(wl_um, n2_coeffs)
    delta_n = n1 * delta_T + n2 * delta_T ** 2
    return n0 + delta_n

def n_z(wl_um, T_C):
    l2 = wl_um ** 2
    n0_squared = 4.59423 + (0.06206 / (l2 - 0.04763)) + (110.80672 / (l2 - 86.12171))
    n0 = np.sqrt(n0_squared)
    n1_coeffs = temp_fit_coeff_n1['z'].values
    n2_coeffs = temp_fit_coeff_n2['z'].values
    delta_T = T_C - 20.0
    n1 = calculate_poly_coeffs(wl_um, n1_coeffs)
    n2 = calculate_poly_coeffs(wl_um, n2_coeffs)
    delta_n = n1 * delta_T + n2 * delta_T ** 2
    return n0 + delta_n

def k_vector(wl_m, n):
    return 2 * pi * n / wl_m

def delta_k(T_C, wl_p_m, poling_period_m):
    wl_s_m = 2 * wl_p_m
    wl_i_m = wl_s_m
    wl_p_um = wl_p_m * 1e6
    wl_s_um = wl_s_m * 1e6
    wl_i_um = wl_i_m * 1e6
    n_p = n_y(wl_p_um, T_C)
    n_s = n_z(wl_s_um, T_C)
    n_i = n_y(wl_i_um, T_C)
    k_p = k_vector(wl_p_m, n_p)
    k_s = k_vector(wl_s_m, n_s)
    k_i = k_vector(wl_i_m, n_i)
    K = 2 * pi / poling_period_m
    return abs(k_p - k_s - k_i - K)

def find_optimal_temperature(wl_p_m, poling_period_m):
    poling_period_um = poling_period_m * 1e6
    print(f"--- Calculating for Poling Period: {poling_period_um:.2f} µm ---")
    objective_function = lambda T_C: delta_k(T_C, wl_p_m, poling_period_m)
    result = minimize_scalar(objective_function, bounds=(0, 100), method='bounded')
    if result.success:
        print(f"Optimum temperature in range [0, 100]°C: {result.x:.2f} °C")
        print(f"Phase mismatch Δk at this temperature: {result.fun:.4e} m⁻¹\n")
    else:
        print("Optimization failed: Could not find a phase-matching temperature in the specified range.\n")

if __name__ == "__main__":
    find_optimal_temperature(
        wl_p_m=lambda_p_meters,
        poling_period_m=poling_period_meters
    )