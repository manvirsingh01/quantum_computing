import numpy as np
from scipy.constants import pi
from scipy.optimize import minimize_scalar

# ================================
# Sellmeier and Temperature Model
# ================================

def calculate_poly_coeffs(wl_um, coeffs):
    """
    Calculates temperature-dependent coefficients n1 or n2 for given wavelength.
    Uses: n(λ) = a0 + a1/λ + a2/λ² + a3/λ³

    Args:
        wl_um (float): Wavelength in micrometers.
        coeffs (tuple): Coefficients (a0, a1, a2, a3).

    Returns:
        float: Resulting coefficient value.
    """
    a0, a1, a2, a3 = coeffs
    return a0 + a1 / wl_um + a2 / wl_um**2 + a3 / wl_um**3


def n_z(wl_um, T_C):
    """
    Calculates the refractive index for z-polarized light in PPKTP crystal.

    Args:
        wl_um (float): Wavelength in micrometers.
        T_C (float): Temperature in Celsius.

    Returns:
        float: Temperature-dependent refractive index.
    """
    # Sellmeier base index (static terms)
    l2 = wl_um**2
    A = 4.59423
    B = 0.06206 / (l2 - 0.04763)
    C = 110.80672 / (l2 - 86.12171)
    n0 = np.sqrt(A + B + C)

    # Temperature-dependent polynomial coefficients from literature
    n1_coeffs = (0.4997e-4, -1.325e-4, 1.233e-4, -0.505e-4)
    n2_coeffs = (0.3423e-7, -1.464e-7, 1.656e-7, -0.739e-7)

    # Temperature correction: Δn(T) = n1·ΔT + n2·ΔT²
    delta_T = T_C - 25.0
    n1 = calculate_poly_coeffs(wl_um, n1_coeffs)
    n2 = calculate_poly_coeffs(wl_um, n2_coeffs)
    delta_n = n1 * delta_T + n2 * delta_T**2

    return n0 + delta_n


# ================================
# Phase Mismatch Calculation
# ================================

def k_vector(wl_um, n):
    """Computes the wave vector magnitude (k = 2πn/λ)."""
    return 2 * pi * n / (wl_um * 1e-6)


def delta_k(T_C, wl_p=0.405, wl_s=0.810, poling_period_um=9.0):
    """
    Computes phase mismatch Δk for given temperature.

    Args:
        T_C (float): Temperature in Celsius.
        wl_p (float): Pump wavelength in microns.
        wl_s (float): Signal wavelength in microns.
        poling_period_um (float): Quasi-phase-matching period in microns.

    Returns:
        float: Absolute value of phase mismatch Δk (in 1/m).
    """
    wl_i = 1 / (1 / wl_p - 1 / wl_s)  # Energy conservation

    n_p = n_z(wl_p, T_C)
    n_s = n_z(wl_s, T_C)
    n_i = n_z(wl_i, T_C)

    k_p = k_vector(wl_p, n_p)
    k_s = k_vector(wl_s, n_s)
    k_i = k_vector(wl_i, n_i)
    K = 2 * pi / (poling_period_um * 1e-6)  # Grating vector

    return abs(k_p - k_s - k_i - K)


# ================================
# Optimization for Phase Matching
# ================================

def find_optimal_temperature():
    result = minimize_scalar(delta_k, bounds=(20, 200), method='bounded')
    if result.success:
        print(f"✅ Optimum temperature for phase matching: {result.x:.2f} °C")
        print(f"🔬 Minimum phase mismatch Δk: {result.fun:.4e} m⁻¹")
    else:
        print("❌ Optimization failed.")
        print(f"Reason: {result.message}")


# ================================
# Run the Script
# ================================

if __name__ == "__main__":
    find_optimal_temperature()
