import numpy as np


# --- Sellmeier Equations for KDP ---

def n_o_kdp(wavelength):
    lambda_sq = wavelength ** 2
    n_o_sq = 2.259276 + (0.0184 / (lambda_sq - 0.0179)) - 0.0155 * lambda_sq
    return np.sqrt(n_o_sq)


def n_e_kdp(wavelength):
    lambda_sq = wavelength ** 2
    n_e_sq = 2.132668 + (0.0128 / (lambda_sq - 0.0156)) - 0.0044 * lambda_sq
    return np.sqrt(n_e_sq)


def calculate_qpm_period(lambda_p, lambda_s, lambda_i):
    n_p = n_e_kdp(lambda_p)
    n_s = n_e_kdp(lambda_s)
    n_i = n_e_kdp(lambda_i)

    if n_p / lambda_p <= n_s / lambda_s + n_i / lambda_i:
        print("Warning: Dispersion prevents QPM.")
        return None

    lambda_inv = (n_p / lambda_p) - (n_s / lambda_s) - (n_i / lambda_i)
    if lambda_inv <= 0:
        return None

    return 1 / lambda_inv


# --- Execution ---

pump_wavelength = 0.405  # μm
signal_wavelength = 2 * pump_wavelength
idler_wavelength = 2 * pump_wavelength

print("Type-0 SPDC in KDP: QPM Grating Period Calculation")
print("-" * 50)
print(f"Pump λ: {pump_wavelength * 1000:.0f} nm")
print(f"Signal/Idler λ: {signal_wavelength * 1000:.0f} nm")
print("-" * 50)

grating_period = calculate_qpm_period(pump_wavelength, signal_wavelength, idler_wavelength)

if grating_period is not None:
    print(f"Required grating period (Λ): {grating_period:.2f} μm")
else:
    print("No valid grating period found.")
