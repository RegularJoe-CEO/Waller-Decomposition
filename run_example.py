"""
run_example.py - Waller Decomposition EUR Forecast
Usage: python3 run_example.py sample_well.csv
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

N_SAMPLES = 500

def load_data(filepath):
    df = pd.read_csv(filepath)
    if "days" not in df.columns:
        raise ValueError("Missing column: days")
    if "pressure_psi" not in df.columns:
        raise ValueError("Missing column: pressure_psi")
    if "rate_mscfd" in df.columns:
        rate = df["rate_mscfd"].values
        well_type, rate_unit, eur_unit = "gas", "MSCF/d", "MMSCF"
    elif "rate_bblpd" in df.columns:
        rate = df["rate_bblpd"].values
        well_type, rate_unit, eur_unit = "oil", "bbl/d", "MBO"
    else:
        raise ValueError("Missing column: rate_mscfd (gas) or rate_bblpd (oil)")
    return df["days"].values, rate, df["pressure_psi"].values, well_type, rate_unit, eur_unit

def power_law_decline(t, qi, D, n):
    return qi * (1 + D * t) ** (-n)

def fit_decline(t, rate):
    try:
        popt, pcov = curve_fit(power_law_decline, t, rate, p0=[rate[0], 0.01, 0.8],
            bounds=([0, 0.0001, 0.1], [rate[0]*3, 0.5, 2.0]), maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    except:
        return [rate[0], 0.01, 0.8], [rate[0]*0.1, 0.005, 0.1]

def monte_carlo_eur(t, rate, n_samples=500, t_max=1800):
    params, errors = fit_decline(t, rate)
    qi, D, n = params
    qi_err = max(errors[0], qi * 0.05)
    D_err = max(errors[1], D * 0.10)
    n_err = max(errors[2], n * 0.05)
    t_forecast = np.linspace(1, t_max, 1000)
    eur_samples = []
    for _ in range(n_samples):
        qi_s = max(qi + qi_err * np.random.randn(), 10)
        D_s = np.clip(D + D_err * np.random.randn(), 0.001, 0.5)
        n_s = np.clip(n + n_err * np.random.randn(), 0.1, 2.0)
        q_forecast = power_law_decline(t_forecast, qi_s, D_s, n_s)
        eur_samples.append(np.trapz(q_forecast, t_forecast) / 1000)
    return np.percentile(eur_samples, [10, 50, 90]), params

def main(filepath):
    print(f"Loading: {filepath}")
    t, rate, pressure, well_type, rate_unit, eur_unit = load_data(filepath)
    print(f"  {len(t)} data points, {t[-1]:.0f} days")
    print(f"  Well type: {well_type.upper()} ({rate_unit})")
    print("Fitting decline and computing uncertainty...")
    (p10, p50, p90), params = monte_carlo_eur(t, rate, N_SAMPLES)
    print("")
    print("=" * 40)
    print(f"EUR FORECAST ({eur_unit})")
    print("=" * 40)
    print(f"  P10:  {p10:>8.1f}")
    print(f"  P50:  {p50:>8.1f}")
    print(f"  P90:  {p90:>8.1f}")
    print("=" * 40)
    print(f"\nFit: qi={params[0]:.0f}, D={params[1]:.4f}, n={params[2]:.2f}")
    t_forecast = np.linspace(1, 1800, 500)
    q_forecast = power_law_decline(t_forecast, *params)
    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    ax.scatter(t, rate, s=30, color="black", alpha=0.7, label="Observed", zorder=5)
    ax.plot(t_forecast, q_forecast, color="blue", linewidth=2, label="P50 Forecast")
    ax.axvline(x=t[-1], color="gray", linestyle=":", linewidth=1, alpha=0.8)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel(f"Rate ({rate_unit})")
    ax.set_title("Waller Decomposition EUR Forecast")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1800)
    ax.set_ylim(0, None)
    ax.text(0.02, 0.02, f"EUR P50: {p50:.1f} {eur_unit}", transform=ax.transAxes,
            fontsize=10, bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    plt.tight_layout()
    plt.savefig("forecast_output.png", dpi=150)
    print("\nSaved: forecast_output.png")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 run_example.py <your_well.csv>")
        sys.exit(1)
    main(sys.argv[1])