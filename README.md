# SPDC Phase-Matching Simulation Tools

This repository contains a collection of Python scripts designed to calculate phase-matching conditions for Spontaneous Parametric Down-Conversion (SPDC) in various nonlinear optical crystals. These tools are essential for designing experiments that generate entangled photon pairs, a fundamental resource in quantum computing and quantum information science.

## Overview.        

Spontaneous Parametric Down-Conversion is a nonlinear optical process where a high-energy pump photon spontaneously splits into two lower-energy photons, conventionally named the "signal" and "idler," inside a nonlinear crystal. For this process to occur efficiently, phase-matching conditions (which ensure momentum conservation) must be met.

This project provides scripts to calculate the required parameters—such as crystal angle, temperature, or poling period—to achieve phase-matching for different SPDC configurations.

## Features

- **Multiple SPDC Types:** Simulates conditions for Type-0, Type-I, and Type-II SPDC.
- **Various Phase-Matching Schemes:** Supports both Birefringent Phase-Matching (BPM) and Quasi-Phase-Matching (QPM), including temperature tuning.
- **Common Nonlinear Crystals:** Includes models for KDP, BBO, and periodically-poled KTP (PPKTP).
- **Key Parameter Calculation:**
  - Phase-matching angles ($\theta$).
  - Quasi-Phase-Matching grating periods ($\Lambda$).
  - Optimal crystal temperature for PPKTP.
  - Effective nonlinear coefficient ($d_{eff}$).

## Scripts Overview.         

| Script Name       | Crystal | SPDC Type         | Phase-Matching | Goal                                                                 |
| ----------------- | ------- | ----------------- | -------------- | -------------------------------------------------------------------- |
| `spdc_type_0.py`  | KDP     | Type-0 (e → e + e)  | QPM            | Calculate the required poling period ($\Lambda$).                      |
| `spdc_type_1.py`  | KDP     | Type-I (e → o + o)  | BPM            | Find the phase-matching angle ($\theta$) and $d_{eff}$.                |
| `spdc_type_2.py`  | BBO     | Type-II (e → o + e) | BPM            | Find the collinear phase-matching angle ($\theta$).                  |
| `ppktp.py`        | PPKTP   | Type-II (y→z+y)     | QPM (Temp)     | Find the optimal crystal temperature for a given poling period.      |

## Theoretical Background

### 1. Spontaneous Parametric Down-Conversion (SPDC)

SPDC is governed by the laws of energy and momentum conservation:

-   **Energy Conservation:** $\omega_p = \omega_s + \omega_i$
    (The pump photon frequency equals the sum of signal and idler frequencies).
-   **Momentum Conservation (Phase-Matching):** $\vec{k_p} = \vec{k_s} + \vec{k_i}$
    (The pump wavevector must equal the sum of the signal and idler wavevectors).

The phase-matching condition is crucial for constructive interference of the down-converted fields, leading to efficient photon pair generation.

### 2. Phase-Matching Techniques

-   **Birefringent Phase-Matching (BPM):** Achieved in birefringent crystals by carefully choosing the angle of the pump beam relative to the crystal's optical axis. This angle is used to tune the refractive indices "seen" by the differently polarized photons. Scripts `spdc_type_1.py` and `spdc_type_2.py` use this method.
-   **Quasi-Phase-Matching (QPM):** Achieved by periodically reversing the crystal's nonlinear coefficient, creating a grating structure. This introduces a grating vector, $\vec{K_g}$, to the momentum equation: $\vec{k_p} = \vec{k_s} + \vec{k_i} + \vec{K_g}$. Script `spdc_type_0.py` uses this method.

### 3. Temperature Tuning in QPM

In materials like Periodically-Poled KTP (PPKTP), the refractive indices are not only dependent on wavelength but also highly sensitive to temperature. This property can be exploited for phase-matching. Instead of changing the crystal angle, the temperature is varied to fine-tune the refractive indices. For a fixed pump wavelength and poling period ($\Lambda$), there is an optimal temperature at which the phase mismatch ($\Delta k$) is minimized. The `ppktp.py` script automates finding this temperature by minimizing the phase-mismatch equation using temperature-dependent Sellmeier equations.

### 4. SPDC Polarization Configurations

-   **Type-0:** The pump, signal, and idler photons all have the same polarization (e.g., e → e + e).
-   **Type-I:** The signal and idler photons have the same polarization, which is orthogonal to the pump (e.g., e → o + o).
-   **Type-II:** The signal and idler photons have orthogonal polarizations (e.g., e → o + e).

*(Note: 'e' denotes extraordinary polarization and 'o' denotes ordinary polarization.)*

## Prerequisites

The scripts require Python 3 and the following packages:

-   NumPy
-   SciPy
-   Pandas (only for `ppktp.py`)

## Installation

1.  Clone the repository:
    ```bash
    git clone [https://github.com/your-username/quantum_computing.git](https://github.com/your-username/quantum_computing.git)
    cd quantum_computing
    ```

2.  Install the required packages using pip:
    ```bash
    pip install numpy scipy pandas
    ```

## Usage

Each script can be run directly from the command line. The parameters (like pump wavelength) are hard-coded in the files but can be easily modified.

### `spdc_type_0.py`

Calculates the QPM grating period for Type-0 SPDC in KDP.

**To Run:**
```bash
python spdc_type_0.py
```

**Sample Output:**
```
Type-0 SPDC in KDP: QPM Grating Period Calculation
--------------------------------------------------
Pump λ: 405 nm
Signal/Idler λ: 810 nm
--------------------------------------------------
Required grating period (Λ): 7.78 μm
```

### `spdc_type_1.py`

Calculates the phase-matching angle and effective nonlinearity for Type-I SPDC in KDP.

**To Run:**
```bash
python spdc_type_1.py
```

**Sample Output:**
```
KDP Type-I SPDC: Phase Matching and d_eff
--------------------------------------------------
Pump: 405 nm → Signal/Idler: 810 nm
--------------------------------------------------
Phase-Matching Angle (θ): 50.25°
Effective Nonlinearity (d_eff): 0.2976 pm/V
```

### `spdc_type_2.py`

Calculates the collinear phase-matching angle for Type-II SPDC in BBO.

**To Run:**
```bash
python spdc_type_2.py
```

**Sample Output:**
```
--- BBO Type-II SPDC Phase-Matching Angle Calculation (Collinear e -> o + e) ---

--- Example 1: Pump @ 405 nm ---
Pump Wavelength: 405 nm
Degenerate Signal/Idler Wavelength: 810 nm
Calculated Collinear Phase-Matching Angle (θp): 42.06 degrees
```

### `ppktp.py`

This script models a Type-II process in Periodically-Poled KTP (PPKTP), where phase-matching is achieved by tuning the crystal's temperature. It uses temperature-dependent Sellmeier equations to find the optimal operating temperature that minimizes the phase mismatch for a given pump wavelength and poling period.

**To Run:**
```bash
python ppktp.py
```

**Sample Output:**
```
--- Calculating for Poling Period: 10.00 µm ---
Optimum temperature in range [0, 100]°C: 60.75 °C
Phase mismatch Δk at this temperature: 1.4870e-09 m⁻¹
```

## Contributing

Contributions are welcome! If you would like to improve a script, add a new material model, or implement a different configuration, please feel free to fork the repository and submit a pull request.

1.  Fork the Project.
2.  Create your Feature Branch (`git checkout -b feature/AmazingFeature`).
3.  Commit your Changes (`git commit -m 'Add some AmazingFeature'`).
4.  Push to the Branch (`git push origin feature/AmazingFeature`).
5.  Open a Pull Request.

## License

This project is distributed under the MIT License. See `LICENSE` for more information.
