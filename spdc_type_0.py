import numpy as np
from scipy.constants import c

# Sellmeier equation for KTP (n_z, assuming all polarizations are z-polarized for Type-0)
def n_z(wl_um):
    l2 = wl_um ** 2
    return np.sqrt(4.59423 + (0.06206 / (l2 - 0.04763)) + (110.80672 / (l2 - 86.12171)))

# Wavelengths (in micrometers)
wl_p = 0.405  # Pump wavelength
wl_s = 0.810  # Signal wavelength
wl_i = 0.810  # From energy conservation

# Refractive indices (all same polarization)
n_p = n_z(wl_p)
n_s = n_z(wl_s)
n_i = n_z(wl_i)

# Wavevectors
k_p = 2 * np.pi * n_p / wl_p
k_s = 2 * np.pi * n_s / wl_s
k_i = 2 * np.pi * n_i / wl_i

# Quasi-Phase Matching condition: k_p = k_s + k_i + 2π / Λ
# Solve for Λ
delta_k = k_p - k_s - k_i
Lambda = 2 * np.pi / delta_k

print("Pump Wavelength (um):", wl_p)
print("Signal Wavelength (um):", wl_s)
print("Idler Wavelength (um):", wl_i)
print("Refractive Indices (n_p, n_s, n_i):", n_p, n_s, n_i)
print("Required Grating Period Λ (um):", Lambda)
