import numpy as np
from scipy.optimize import minimize_scalar

# Speed of light in vacuum (m/s)
c = 299792458

# Sellmeier equations for BBO crystal
# These functions return the refractive indices for ordinary and extraordinary rays
# Wavelength (wl_um) is in micrometers (um)
def sellmeier_o(wl_um):
    l2 = wl_um ** 2
    return np.sqrt(2.7359 + 0.01878 / (l2 - 0.01822) - 0.01354 * l2)

def sellmeier_e(wl_um):
    l2 = wl_um ** 2
    return np.sqrt(2.3753 + 0.01224 / (l2 - 0.01667) - 0.01516 * l2)

# Extraordinary refractive index depends on angle theta (rad) with respect to the optical axis
def ne_theta(no, ne, theta_rad):
    sin2 = np.sin(theta_rad) ** 2
    cos2 = np.cos(theta_rad) ** 2
    return (ne * no) / np.sqrt(ne**2 * cos2 + no**2 * sin2)

# Define wavelengths (in micrometers, um)
wl_p = 0.4   # Pump wavelength (extraordinary polarized)
wl_s = 0.8   # Signal wavelength (ordinary polarized)
wl_i = 1 / (1 / wl_p - 1 / wl_s)  # Idler wavelength from energy conservation (extraordinary polarized)

# Refractive indices for the signal, idler, and pump waves
no_s = sellmeier_o(wl_s)             # Signal is ordinary polarized
ne_i = sellmeier_e(wl_i)             # Idler is extraordinary polarized
no_p = sellmeier_o(wl_p)             # Used for computing angular dependence

# Phase mismatch function: |Δk| = |k_p - k_s - k_i|
def phase_mismatch(theta_deg):
    theta_rad = np.radians(theta_deg)
    ne_p_theta = ne_theta(no_p, sellmeier_e(wl_p), theta_rad)  # Pump extraordinary index as function of angle

    # Wavevector magnitudes: k = 2π * n / λ
    k_p = ne_p_theta * 2 * np.pi / (wl_p * 1e-6)
    k_s = no_s * 2 * np.pi / (wl_s * 1e-6)
    k_i = ne_i * 2 * np.pi / (wl_i * 1e-6)

    delta_k = k_p - k_s - k_i  # Vector momentum conservation condition
    return abs(delta_k)

# Minimize |Δk| to find the phase-matching angle theta (in degrees)
result = minimize_scalar(phase_mismatch, bounds=(20, 40), method='bounded')

# Final results
optimal_theta = result.x        # Optimal pump angle for phase matching
min_delta_k = result.fun        # Corresponding minimum phase mismatch Δk

print(f"Optimal phase-matching angle (theta): {optimal_theta:.4f} degrees")
print(f"Minimum phase mismatch (Delta k): {min_delta_k:.2f} m^-1")
