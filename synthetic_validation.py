"""
Synthetic Validation for Waller Decomposition
Tests if the method can recover known EUR from synthetic data.
"""

import numpy as np
import pywt
from scipy.optimize import curve_fit
from scipy.integrate import quad

# --- Core Functions (same as sensitivity test) ---

def power_law_decline(t, qi, D, n):
    """Generalized power-law decline."""
    return qi * (1 + D * t) ** (-n)

def true_eur_analytical(qi, D, n, t_end):
    """Calculate true EUR analytically for power-law decline."""
    if 0 < n < 1:
        return qi / (D * (1 - n)) * ((1 + D * t_end) ** (1 - n) - 1)
    else:
        result, _ = quad(lambda x: power_law_decline(x, qi, D, n), 0, t_end)
        return result

def reconstruct_dwt_components(signal, wavelet='db4', level=4):
    """Decompose and reconstruct each DWT level independently."""
    coeffs = pywt.wavedec(signal, wavelet, level=level)
    components = []
    for i in range(len(coeffs)):
        coeffs_zeroed = [np.zeros_like(c) for c in coeffs]
        coeffs_zeroed[i] = coeffs[i].copy()
        component = pywt.waverec(coeffs_zeroed, wavelet)
        components.append(component[:len(signal)])
    return components

def separate_rate_components(t, q, p_wf, p_i, wavelet='db4', level=4):
    """Separate production rate into multi-scale components."""
    delta_p = p_i - p_wf
    rnp = delta_p / q
    rnp_log = np.log10(np.maximum(rnp, 1e-6))

    rnp_components_log = reconstruct_dwt_components(rnp_log, wavelet, level)
    rnp_components = [10 ** comp for comp in rnp_components_log]

    inv_rnp = [1.0 / np.maximum(comp, 1e-6) for comp in rnp_components]
    inv_rnp_sum = np.sum(inv_rnp, axis=0)

    q_components = [q * (inv_comp / inv_rnp_sum) for inv_comp in inv_rnp]
    return q_components

def fit_and_compute_eur(t, q_components, t_abandon=10950, q_min=10):
    """Fit decline to each component and sum EUR."""
    eur_total = 0.0

    for q_k in q_components:
        if np.mean(q_k) < 0.01 * np.mean(np.sum(q_components, axis=0)):
            continue
        try:
            popt, _ = curve_fit(
                power_law_decline, t, q_k,
                p0=[q_k[0], 0.01, 0.5],
                bounds=([1e-6, 1e-6, 0.01], [q_k[0] * 10, 1.0, 2.0]),
                maxfev=10000
            )
            qi, D, n = popt

            if n > 0 and D > 0:
                t_econ = ((qi / q_min) ** (1/n) - 1) / D
                t_end = min(t_abandon, t_econ)
            else:
                t_end = t_abandon

            if 0 < n < 1:
                eur_k = qi / (D * (1 - n)) * ((1 + D * t_end) ** (1 - n) - 1)
            else:
                eur_k, _ = quad(lambda x: power_law_decline(x, qi, D, n), 0, t_end)

            eur_total += eur_k
        except:
            continue

    return eur_total

# --- Define Synthetic Well with KNOWN EUR ---

print("=" * 65)
print("WALLER DECOMPOSITION: SYNTHETIC VALIDATION")
print("=" * 65)

# True decline parameters for 3 stages
stages = [
    {"qi": 2000, "D": 0.015, "n": 0.7, "name": "Stage 1 (near-wellbore)"},
    {"qi": 1000, "D": 0.008, "n": 0.5, "name": "Stage 2 (mid-field)"},
    {"qi": 500,  "D": 0.004, "n": 0.4, "name": "Stage 3 (far-field)"},
]

t_abandon = 10950  # 30 years in days
q_min = 10  # Economic limit MSCF/d

# Calculate TRUE EUR for each stage
print("\n--- TRUE WELL PARAMETERS ---\n")
true_eur_total = 0
for stage in stages:
    qi, D, n = stage["qi"], stage["D"], stage["n"]

    # Economic limit time
    if qi > q_min:
        t_econ = ((qi / q_min) ** (1/n) - 1) / D
    else:
        t_econ = t_abandon
    t_end = min(t_abandon, t_econ)

    eur_stage = true_eur_analytical(qi, D, n, t_end) / 1000  # MMSCF
    true_eur_total += eur_stage
    print(f"{stage['name']}: qi={qi}, D={D}, n={n}")
    print(f"  → EUR = {eur_stage:.2f} MMSCF (t_econ = {t_end:.0f} days)\n")

print(f"TRUE TOTAL EUR: {true_eur_total:.2f} MMSCF")
print("-" * 65)

# --- Generate Synthetic Production Data ---

np.random.seed(42)
t = np.linspace(1, 180, 180)  # 180 days observed

# Sum of all stages
q_true = sum(
    stage["qi"] * (1 + stage["D"] * t) ** (-stage["n"])
    for stage in stages
)

# Add realistic noise (3%)
noise_level = 0.03
q_observed = q_true * (1 + noise_level * np.random.randn(len(t)))
q_observed = np.maximum(q_observed, 10)

# Synthetic BHP (declining)
p_i = 5500
p_wf = 2500 - 3 * t + 15 * np.random.randn(len(t))
p_wf = np.clip(p_wf, 500, p_i - 100)

print(f"\n--- OBSERVED DATA ---\n")
print(f"Days observed: {len(t)}")
print(f"Initial rate: {q_observed[0]:.1f} MSCF/d")
print(f"Final rate: {q_observed[-1]:.1f} MSCF/d")
print(f"Noise level: {noise_level*100:.0f}%")
print(f"Initial pressure: {p_i} psi")
print("-" * 65)

# --- Run Waller Decomposition ---

print(f"\n--- WALLER DECOMPOSITION RESULTS ---\n")

q_components = separate_rate_components(t, q_observed, p_wf, p_i, wavelet='db4', level=4)
estimated_eur = fit_and_compute_eur(t, q_components, t_abandon, q_min) / 1000  # MMSCF

error_pct = (estimated_eur - true_eur_total) / true_eur_total * 100

print(f"True EUR:      {true_eur_total:.2f} MMSCF")
print(f"Estimated EUR: {estimated_eur:.2f} MMSCF")
print(f"Error:         {error_pct:+.2f}%")
print("-" * 65)

# --- Test Multiple Noise Realizations ---

print(f"\n--- MONTE CARLO: 100 NOISE REALIZATIONS ---\n")

errors = []
estimates = []

for seed in range(100):
    np.random.seed(seed)
    q_noisy = q_true * (1 + noise_level * np.random.randn(len(t)))
    q_noisy = np.maximum(q_noisy, 10)

    p_wf_noisy = 2500 - 3 * t + 15 * np.random.randn(len(t))
    p_wf_noisy = np.clip(p_wf_noisy, 500, p_i - 100)

    try:
        q_comp = separate_rate_components(t, q_noisy, p_wf_noisy, p_i, wavelet='db4', level=4)
        eur_est = fit_and_compute_eur(t, q_comp, t_abandon, q_min) / 1000
        estimates.append(eur_est)
        errors.append((eur_est - true_eur_total) / true_eur_total * 100)
    except:
        continue

estimates = np.array(estimates)
errors = np.array(errors)

print(f"Realizations: {len(estimates)}")
print(f"True EUR:     {true_eur_total:.2f} MMSCF")
print(f"Mean Est:     {np.mean(estimates):.2f} MMSCF")
print(f"Std Dev:      {np.std(estimates):.2f} MMSCF")
print(f"P10:          {np.percentile(estimates, 10):.2f} MMSCF")
print(f"P50:          {np.percentile(estimates, 50):.2f} MMSCF")
print(f"P90:          {np.percentile(estimates, 90):.2f} MMSCF")
print(f"\nMean Error:   {np.mean(errors):+.2f}%")
print(f"Error Std:    {np.std(errors):.2f}%")
print("-" * 65)

# --- Verdict ---

print(f"\n{'=' * 65}")
if abs(np.mean(errors)) < 5 and np.std(errors) < 10:
    print("✓ VALIDATED: Waller Decomposition recovers true EUR within 5% mean error")
elif abs(np.mean(errors)) < 10:
    print("⚠ PARTIAL: Mean error < 10%, but consider limitations")
else:
    print("✗ FAILED: Method shows significant bias")
print("=" * 65)
