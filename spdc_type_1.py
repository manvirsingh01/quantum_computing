import numpy as np
from scipy.optimize import fsolve

# --- Refractive Index Model (KDP) ---

def n_o_kdp(wavelength):
    lambda_sq = wavelength**2
    n_o_sq = 2.259276 + (0.0184 / (lambda_sq - 0.0179)) - 0.0155 * lambda_sq
    return np.sqrt(n_o_sq)

def n_e_kdp(wavelength):
    lambda_sq = wavelength**2
    n_e_sq = 2.132668 + (0.0128 / (lambda_sq - 0.0156)) - 0.0044 * lambda_sq
    return np.sqrt(n_e_sq)

def n_e_theta(n_o, n_e, theta_rad):
    cos2 = np.cos(theta_rad)**2
    sin2 = np.sin(theta_rad)**2
    n_e_theta_sq = (n_e**2 * n_o**2) / (n_e**2 * cos2 + n_o**2 * sin2)
    return np.sqrt(n_e_theta_sq)

# --- Phase Matching Condition ---

def phase_mismatch(theta_deg, lambda_p, lambda_s, lambda_i):
    theta_rad = np.radians(theta_deg)
    n_o_p = n_o_kdp(lambda_p)
    n_e_p = n_e_kdp(lambda_p)
    n_p_eff = n_e_theta(n_o_p, n_e_p, theta_rad)
    n_s = n_o_kdp(lambda_s)
    n_i = n_o_kdp(lambda_i)
    return (n_p_eff / lambda_p) - (n_s / lambda_s) - (n_i / lambda_i)

# --- Effective Nonlinearity ---

def calculate_deff_kdp(theta_deg, phi_deg):
    d36 = 0.39  # pm/V
    theta_rad = np.radians(theta_deg)
    phi_rad = np.radians(phi_deg)
    return d36 * np.sin(2 * phi_rad) * np.sin(theta_rad)

# --- Execution ---

lambda_pump = 0.405  # μm
lambda_signal = 2 * lambda_pump
lambda_idler = 2 * lambda_pump

phi_angle_deg = 45.0
initial_theta = 40.0

print("KDP Type-I SPDC: Phase Matching and d_eff")
print("-" * 50)
print(f"Pump: {lambda_pump * 1000:.0f} nm → Signal/Idler: {lambda_signal * 1000:.0f} nm")
print("-" * 50)

solution = fsolve(phase_mismatch, initial_theta, args=(lambda_pump, lambda_signal, lambda_idler))
theta_pm = solution[0]

deff = calculate_deff_kdp(theta_pm, phi_angle_deg)

print(f"Phase-Matching Angle (θ): {theta_pm:.2f}°")
print(f"Effective Nonlinearity (d_eff): {deff:.4f} pm/V")
