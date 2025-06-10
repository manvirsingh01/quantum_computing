import numpy as np
from scipy.optimize import minimize_scalar

# Constants
c = 299_792_458  # Speed of light in m/s

# Wavelengths (in micrometers)
wl_p = 0.405  # Pump
wl_s = 0.810  # Signal
wl_i = 1 / (1 / wl_p - 1 / wl_s)  # Idler (energy conservation)

# Grating period (in micrometers)
Lambda = 9.0

# Sellmeier equations with temperature dependence for PPKTP
def n_y(wl, T):
    l2 = wl ** 2
    return np.sqrt(
        3.29100 + 0.04140 / (l2 - 0.03978) + 9.35522 / (l2 - 31.45571)
        + 1e-6 * (0.00094186 * l2 - 0.0000002937 * l2**2) * T
    )

def n_z(wl, T):
    l2 = wl ** 2
    return np.sqrt(
        4.59423 - 0.06206 / (l2 - 0.04763) - 110.80672 / (l2 - 86.12171)
        + 1e-6 * (0.0011 * l2 - 0.00000035 * l2**2) * T
    )

# Phase mismatch function Δk(T)
def phase_mismatch(T):
    n_p = n_y(wl_p, T)
    n_s = n_y(wl_s, T)
    n_i = n_z(wl_i, T)

    # Wavevectors (in m⁻¹)
    k_p = 2 * np.pi * n_p / (wl_p * 1e-6)
    k_s = 2 * np.pi * n_s / (wl_s * 1e-6)
    k_i = 2 * np.pi * n_i / (wl_i * 1e-6)

    # QPM grating vector
    K_qpm = 2 * np.pi / (Lambda * 1e-6)

    # Phase mismatch Δk
    delta_k = k_p - k_s - k_i - K_qpm
    return abs(delta_k)

# Minimize |Δk| to find phase matching temperature
result = minimize_scalar(phase_mismatch, bounds=(20, 100), method='bounded')

# Extract optimal temperature and mismatch
optimal_T = result.x
min_delta_k = result.fun

# Threshold to consider Δk ≈ 0 (phase matching)
tolerance = 1e-5

# Variable a will hold the optimal temperature if condition is satisfied
a = None
if min_delta_k < tolerance:
    a = optimal_T

# Output
print("\n=== Type-I SPDC Phase Matching in PPKTP ===")
print(f"Pump wavelength (λₚ):   {wl_p:.4f} µm")
print(f"Signal wavelength (λₛ): {wl_s:.4f} µm")
print(f"Idler wavelength (λᵢ):  {wl_i:.4f} µm")
print(f"Grating period (Λ):     {Lambda:.4f} µm")
print(f"Optimal temperature:    {optimal_T:.2f} °C")
print(f"Minimum |Δk|:           {min_delta_k:.4e} m⁻¹")

if a is not None:
    print(f"✔ Phase matching condition satisfied. Value of a = {a:.2f} °C")
else:
    print("✘ Phase matching not satisfied within tolerance. a is not assigned.")
