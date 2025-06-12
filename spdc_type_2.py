import numpy as np
from scipy.optimize import fsolve


def n_o_bbo(wavelength):
    """
    Calculates the ordinary refractive index (no) of BBO.
    Wavelength must be in micrometers (μm).
    """
    wavelength_sq = wavelength ** 2
    return np.sqrt(2.7359 + 0.01878 / (wavelength_sq - 0.01822) - 0.01354 * wavelength_sq)


def n_e_bbo(wavelength):
    """
    Calculates the principal extraordinary refractive index (ne) of BBO.
    Wavelength must be in micrometers (μm).
    """
    wavelength_sq = wavelength ** 2
    return np.sqrt(2.3753 + 0.01224 / (wavelength_sq - 0.01667) - 0.01516 * wavelength_sq)


def n_e_theta(n_o, n_e, theta_rad):
    """
    Calculates the effective extraordinary refractive index for a given angle (theta)
    to the optical axis. Angle must be in radians.
    """
    cos2_theta = np.cos(theta_rad) ** 2
    sin2_theta = np.sin(theta_rad) ** 2
    n_e_theta_sq = (n_e ** 2 * n_o ** 2) / (n_e ** 2 * cos2_theta + n_o ** 2 * sin2_theta)
    return np.sqrt(n_e_theta_sq)


def phase_mismatch(theta_p_deg, lambda_p, lambda_s, lambda_i):
    """
    Calculates the phase mismatch (delta_k) for collinear Type-II SPDC (e -> o + e).
    The function is zero when phase matching is achieved.
    - theta_p_deg: Pump angle in degrees.
    - Wavelengths are in micrometers (μm).
    """
    theta_p_rad = np.radians(theta_p_deg)

    # Get principal refractive indices at each wavelength
    n_o_p = n_o_bbo(lambda_p)
    n_e_p = n_e_bbo(lambda_p)
    n_o_s = n_o_bbo(lambda_s)
    # n_e_s is not needed for the 'o' polarized signal
    n_o_i = n_o_bbo(lambda_i)
    n_e_i = n_e_bbo(lambda_i)

    # Calculate effective refractive index for 'e' polarized photons at angle theta
    n_p_eff = n_e_theta(n_o_p, n_e_p, theta_p_rad)
    n_i_eff = n_e_theta(n_o_i, n_e_i, theta_p_rad)
    # Signal is 'o' polarized, so its index is just n_o_s
    n_s_eff = n_o_s

    # Momentum conservation for the collinear case: k_p = k_s + k_i
    # This simplifies to: n_p/lambda_p - n_s/lambda_s - n_i/lambda_i = 0
    delta_k = (n_p_eff / lambda_p) - (n_s_eff / lambda_s) - (n_i_eff / lambda_i)

    return delta_k


# --- Main Calculation ---
# Use the examples provided in the text to find the phase-matching angle.

print("--- BBO Type-II SPDC Phase-Matching Angle Calculation (Collinear e -> o + e) ---")
# --- Example 1: Based on the 405 nm pump example ---
lambda_p2 = 0.405  # 405 nm pump
lambda_s2 = 2 * lambda_p2  # 810 nm signal (degenerate)
lambda_i2 = 2 * lambda_p2  # 810 nm idler (degenerate)
initial_angle_guess2 = 30.0  # An initial guess

solution1 = fsolve(phase_mismatch, initial_angle_guess2, args=(lambda_p2, lambda_s2, lambda_i2))
phase_matching_angle2 = solution1[0]

print(f"\n--- Example 1: Pump @ {lambda_p2 * 1000:.0f} nm ---")
print(f"Pump Wavelength: {lambda_p2 * 1000:.0f} nm")
print(f"Degenerate Signal/Idler Wavelength: {lambda_s2 * 1000:.0f} nm")
print(f"Calculated Collinear Phase-Matching Angle (θp): {phase_matching_angle2:.2f} degrees")
print(
    "(Note: The text's 3° separation between signal/idler refers to a non-collinear setup, which requires a more complex vector-based calculation.)")