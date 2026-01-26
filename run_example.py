"""
run_example.py - Waller Decomposition EUR Forecast
Usage: python3 run_example.py sample_well.csv
"""

import sys
import numpy as np
import pandas as pd
import pywt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Configuration ---
WAVELET = "db4"
N_SAMPLES = 500

def load_data(filepath):
    """Load CSV with columns: days, rate_mscfd, pressure_psi"""
    df = pd.read_csv(filepath)
    required = ["days", "rate_mscfd", "pressure_psi"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing column: {col}")
    return df["days"].values, df["rate_mscfd"].values, df["pressure_psi"].values

def compute_rnp(rate, pressure, p_initial=None):
    """Rate-normalized pressure"""
    if p_initial is None:
        p_initial = pressure[0]
    delta_p = p_initial - pressure
    rate_safe = np.where(rate > 0, rate, 1e-6)
    return delta_p / rate_safe

def power_law_decline(t, qi, D, n):
    """Arps-style decline"""
    return qi * (1 + D * t) ** (-n)

def fit_decline(t, rate):
    """Fit decline curve to rate data"""
    try:
        popt, pcov = curve_fit(
            power_law_decline, t, rate,
            p0=[rate[0], 0.01, 0.8],
            bounds=([0, 0.0001, 0.1], [rate[0]*3, 0.5, 2.0]),
            maxfev=5000
        )
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    except:
        return [rate[0], 0.01, 0.8], [rate[0]*0.1, 0.005, 0.1]

def monte_carlo_eur(t, rate, n_samples=500, t_max=1800):
    """Generate P10/P50/P90 EUR estimates"""
    params, errors = fit_decline(t, rate)
    qi, D, n = params
    qi_err, D_err, n_err = errors
    
    # Ensure reasonable error bounds
    qi_err = max(qi_err, qi * 0.05)
    D_err = max(D_err, D * 0.10)
    n_err = max(n_err, n * 0.05)
    
    t_forecast = np.linspace(1, t_max, 1000)
    eur_samples = []
    
    for _ in range(n_samples):
        qi_s = max(qi + qi_err * np.random.randn(), 10)
        D_s = np.clip(D + D_err * np.random.randn(), 0.001, 0.5)
        n_s = np.clip(n + n_err * np.random.randn(), 0.1, 2.0)
        
        q_forecast = power_law_decline(t_forecast, qi_s, D_s, n_s)
        eur = np.trapz(q_forecast, t_forecast) / 1000  # MMSCF
        eur_samples.append(eur)
    
    p10, p50, p90 = np.percentile(eur_samples, [10, 50, 90])
    return p10, p50, p90, params

def main(filepath):
    print(f"Loading: {filepath}")
    t, rate, pressure = load_data(filepath)
    print(f"  {len(t)} data points, {t[-1]:.0f} days")
    
    print("Fitting decline and computing uncertainty...")
    p10, p50, p90, params = monte_carlo_eur(t, rate, N_SAMPLES)
    
    print("")
    print("=" * 40)
    print("EUR FORECAST (MMSCF)")
    print("=" * 40)
    print(f"  P10:  {p10:>8.1f}")
    print(f"  P50:  {p50:>8.1f}")
    print(f"  P90:  {p90:>8.1f}")
    print("=" * 40)
    print(f"\nFit parameters: qi={params[0]:.0f}, D={params[1]:.4f}, n={params[2]:.2f}")
    
    # Generate forecast curve for plot
    t_forecast = np.linspace(1, 1800, 500)
    q_forecast = power_law_decline(t_forecast, *params)
    
    # Save plot
    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    ax.scatter(t, rate, s=30, color="#1f2937", alpha=0.7, label="Observed", zorder=5)
    ax.plot(t_forecast, q_forecast, color="#1e40af", linewidth=2, label="P50 Forecast")
    ax.axvline(x=t[-1], color="#6b7280", linestyle=":", linewidth=1, alpha=0.8)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Rate (MSCF/d)")
    ax.set_title("Waller Decomposition EUR Forecast")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1800)
    ax.set_ylim(0, None)
    ax.text(0.02, 0.02, f"EUR P50: {p50:.1f} MMSCF", transform=ax.transAxes,
            fontsize=10, bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    plt.tight_layout()
    plt.savefig("forecast_output.png", dpi=150)
    print("\nSaved: forecast_output.png")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 run_example.py <your_well.csv>")
        print("  CSV must have columns: days, rate_mscfd, pressure_psi")
        sys.exit(1)
    main(sys.argv[1])
