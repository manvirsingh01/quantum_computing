import numpy as np
from scipy.constants import c, pi
from scipy.optimize import minimize_scalar


# Convert nm to micrometers
def nm_to_um(nm):
    return nm / 1000


# Temperature-dependent Sellmeier equation for z-polarization in PPKTP
def n_z(wl_um, T):
    l2 = wl_um ** 2
    # Static Sellmeier terms
    A = 4.59423
    B = 0.06206 / (l2 - 0.04763)
    C = 110.80672 / (l2 - 86.12171)
    n0 = np.sqrt(A + B + C)

    # Temperature-dependent correction (coefficients from literature)
    # These are example coefficients (adjust as per experimental data if needed)
    a0, a1, a2, a3 = 1.0e-5, -2.5e-6, 1.0e-7, -1.0e-9
    delta_n = (a0 + a1 / wl_um + a2 / wl_um ** 2 + a3 / wl_um ** 3) * (T - 25)
    return n0 + delta_n


# Calculate wave vector magnitude
def k(wavelength_um, n):
    return 2 * pi * n / (wavelength_um * 1e-6)


# Δk (phase mismatch) function
def delta_k(T, wl_p=0.405, wl_s=0.810, Lambda=9.0):  # Lambda in μm
    wl_i = 1 / (1 / wl_p - 1 / wl_s)
    n_p = n_z(wl_p, T)
    n_s = n_z(wl_s, T)
    n_i = n_z(wl_i, T)

    k_p = k(wl_p, n_p)
    k_s = k(wl_s, n_s)
    k_i = k(wl_i, n_i)

    K = 2 * pi / (Lambda * 1e-6)  # grating vector in m⁻¹

    return abs(k_p - k_s - k_i - K)


# Minimize phase mismatch over temperature range
result = minimize_scalar(delta_k, bounds=(20, 1000000), method='bounded')

if result.success:
    optimal_temp = result.x
    print(f"Optimum temperature for phase matching: {optimal_temp:.2f} °C")
    print(f"Phase mismatch Δk at optimum T: {result.fun:.4e} m⁻¹")
else:
    print("Failed to find optimum temperature.")
