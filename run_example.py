"""
run_example.py - Waller Decomposition EUR Forecast
Usage: 
  Single well:  python3 run_example.py sample_well.csv
  Batch mode:   python3 run_example.py --batch ./wells/
"""

import sys
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

N_SAMPLES = 500
ECON_LIMIT_GAS = 50
ECON_LIMIT_OIL = 5
ECON_LIMIT_WATER = 10
T_MAX = 1800
FLOWBACK_DAYS = 60

def load_data(filepath):
    df = pd.read_csv(filepath)
    if "days" not in df.columns:
        raise ValueError("Missing column: days")
    
    data = {"days": df["days"].values}
    well_type = None
    
    if "rate_mscfd" in df.columns:
        data["gas"] = df["rate_mscfd"].values
        well_type = "gas"
    if "rate_bblpd" in df.columns:
        data["oil"] = df["rate_bblpd"].values
        well_type = "oil"
    if "rate_bwpd" in df.columns:
        data["water"] = df["rate_bwpd"].values
    if "pressure_psi" in df.columns:
        data["pressure"] = df["pressure_psi"].values
    
    if "gas" not in data and "oil" not in data:
        raise ValueError("Need rate_mscfd (gas) or rate_bblpd (oil)")
    
    if "gas" in data and "oil" in data:
        well_type = "both"
    
    return data, well_type

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

def forecast_phase(t, rate, econ_limit, n_samples=500, t_max=1800):
    params, errors = fit_decline(t, rate)
    qi, D, n = params
    qi_err = max(errors[0], qi * 0.05)
    D_err = max(errors[1], D * 0.10)
    n_err = max(errors[2], n * 0.05)
    
    t_forecast = np.linspace(1, t_max, 1000)
    eur_samples = []
    q_samples = []
    
    for _ in range(n_samples):
        qi_s = max(qi + qi_err * np.random.randn(), 10)
        D_s = np.clip(D + D_err * np.random.randn(), 0.001, 0.5)
        n_s = np.clip(n + n_err * np.random.randn(), 0.1, 2.0)
        q_forecast = power_law_decline(t_forecast, qi_s, D_s, n_s)
        q_forecast = np.where(q_forecast >= econ_limit, q_forecast, 0)
        q_samples.append(q_forecast)
        eur_samples.append(np.trapz(q_forecast, t_forecast) / 1000)
    
    q_samples = np.array(q_samples)
    q_p10 = np.percentile(q_samples, 10, axis=0)
    q_p50 = np.percentile(q_samples, 50, axis=0)
    q_p90 = np.percentile(q_samples, 90, axis=0)
    
    p10, p50, p90 = np.percentile(eur_samples, [10, 50, 90])
    
    return {
        "p10": p10, "p50": p50, "p90": p90, 
        "params": params, "t": t_forecast,
        "q_p10": q_p10, "q_p50": q_p50, "q_p90": q_p90
    }

def run_single_well(filepath, save_csv=True, save_plot=True, quiet=False):
    if not quiet:
        print(f"Loading: {filepath}")
    data, well_type = load_data(filepath)
    t = data["days"]
    
    flowback_mask = t <= FLOWBACK_DAYS
    n_flowback = np.sum(flowback_mask)
    
    if not quiet:
        print(f"  {len(t)} data points, {t[-1]:.0f} days")
        if n_flowback > 0:
            print(f"  Flowback period: {n_flowback} points ({FLOWBACK_DAYS} days)")
    
    results = {}
    
    if "gas" in data:
        results["gas"] = forecast_phase(t, data["gas"], ECON_LIMIT_GAS, N_SAMPLES, T_MAX)
    if "oil" in data:
        results["oil"] = forecast_phase(t, data["oil"], ECON_LIMIT_OIL, N_SAMPLES, T_MAX)
    if "water" in data:
        results["water"] = forecast_phase(t, data["water"], ECON_LIMIT_WATER, N_SAMPLES, T_MAX)
    
    if not quiet:
        print("")
        print("=" * 50)
        print("EUR FORECAST")
        print("=" * 50)
        if "gas" in results:
            r = results["gas"]
            print(f"  Gas (MMSCF):   P10={r['p10']:>7.1f}  P50={r['p50']:>7.1f}  P90={r['p90']:>7.1f}")
        if "oil" in results:
            r = results["oil"]
            print(f"  Oil (MBO):     P10={r['p10']:>7.1f}  P50={r['p50']:>7.1f}  P90={r['p90']:>7.1f}")
        if "water" in results:
            r = results["water"]
            print(f"  Water (MBW):   P10={r['p10']:>7.1f}  P50={r['p50']:>7.1f}  P90={r['p90']:>7.1f}")
        print("=" * 50)
        
        if "gas" in data and "oil" in data:
            gor = (data["gas"] * 1000) / np.where(data["oil"] > 0, data["oil"], 0.001)
            print(f"\nGOR: Initial={gor[0]:.0f}  Current={gor[-1]:.0f} scf/bbl")
        if "water" in data and "oil" in data:
            wor = data["water"] / np.where(data["oil"] > 0, data["oil"], 0.001)
            print(f"WOR: Initial={wor[0]:.1f}  Current={wor[-1]:.1f} bbl/bbl")
        if "water" in data and "gas" in data and "oil" not in data:
            wgr = data["water"] / np.where(data["gas"] > 0, data["gas"], 0.001) * 1000
            print(f"WGR: Initial={wgr[0]:.1f}  Current={wgr[-1]:.1f} bbl/MMscf")
    
    basename = os.path.splitext(os.path.basename(filepath))[0]
    
    if save_csv:
        csv_rows = []
        t_out = results[list(results.keys())[0]]["t"]
        for i, day in enumerate(t_out):
            row = {"day": int(day)}
            for phase in results:
                row[f"{phase}_p10"] = results[phase]["q_p10"][i]
                row[f"{phase}_p50"] = results[phase]["q_p50"][i]
                row[f"{phase}_p90"] = results[phase]["q_p90"][i]
            csv_rows.append(row)
        
        csv_df = pd.DataFrame(csv_rows)
        csv_path = f"{basename}_forecast.csv"
        csv_df.to_csv(csv_path, index=False)
        if not quiet:
            print(f"\nSaved: {csv_path}")
    
    if save_plot:
        # ARIES/PHDWin style: single combined semi-log plot
        plt.rcParams.update({
            'font.family': 'sans-serif',
            'font.size': 11,
            'axes.linewidth': 0.8,
            'grid.alpha': 0.3
        })
        
        fig, ax = plt.subplots(figsize=(10, 7), dpi=150)
        
        # Industry standard colors
        colors = {"gas": "#DC143C", "oil": "#228B22", "water": "#4169E1"}
        phase_labels = {"gas": "Gas (MSCF/d)", "oil": "Oil (bbl/d)", "water": "Water (bwpd)"}
        eur_units = {"gas": "MMSCF", "oil": "MBO", "water": "MBW"}
        
        # Convert days to months for x-axis
        t_months = t / 30.44
        T_MAX_MONTHS = T_MAX / 30.44
        
        # Plot each phase
        for phase in results.keys():
            r = results[phase]
            t_forecast_months = r["t"] / 30.44
            
            # P50 forecast line (solid)
            ax.semilogy(t_forecast_months, r["q_p50"], color=colors[phase], linewidth=2, 
                       label=f'{phase.capitalize()} P50 - EUR: {r["p50"]:.0f} {eur_units[phase]}')
            
            # P10/P90 lines (dashed, thinner)
            ax.semilogy(t_forecast_months, r["q_p10"], color=colors[phase], linewidth=1, 
                       linestyle='--', alpha=0.6)
            ax.semilogy(t_forecast_months, r["q_p90"], color=colors[phase], linewidth=1, 
                       linestyle='--', alpha=0.6)
            
            # Historical data points (black/gray circles)
            obs_production = data[phase][~flowback_mask]
            t_production_months = t[~flowback_mask] / 30.44
            ax.semilogy(t_production_months, obs_production, 'o', color='black', 
                       markersize=4, alpha=0.7, zorder=5)
            
            # Flowback points (lighter gray)
            if n_flowback > 0:
                obs_flowback = data[phase][flowback_mask]
                t_flowback_months = t[flowback_mask] / 30.44
                ax.semilogy(t_flowback_months, obs_flowback, 'o', color='gray', 
                           markersize=3, alpha=0.5, zorder=4)
        
        # Axis formatting
        ax.set_xlabel('Time (months)', fontsize=12, fontweight='medium')
        ax.set_ylabel('Rate', fontsize=12, fontweight='medium')
        ax.set_xlim(0, T_MAX_MONTHS)
        ax.set_ylim(bottom=1)
        
        # Light gray horizontal grid only (PHDWin style)
        ax.yaxis.grid(True, linestyle='-', alpha=0.3, color='gray')
        ax.xaxis.grid(False)
        
        # Clean legend in upper right
        ax.legend(loc='upper right', framealpha=0.95, fontsize=10)
        
        # EUR callout box
        eur_text = "EUR (P50)\n"
        for phase in results.keys():
            eur_text += f"{phase.capitalize()}: {results[phase]['p50']:.0f} {eur_units[phase]}\n"
        ax.text(0.02, 0.02, eur_text.strip(), transform=ax.transAxes, fontsize=10,
                verticalalignment='bottom', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.9))
        
        # Title
        plot_title = basename.replace('_', ' ').title()
        ax.set_title(plot_title, fontsize=14, fontweight='bold', pad=10)
        
        plt.tight_layout()
        plot_path = f"{basename}_forecast.png"
        plt.savefig(plot_path, dpi=150, facecolor='white', edgecolor='none')
        plt.close()
        if not quiet:
            print(f"Saved: {plot_path}")

    
    return {
        "file": basename,
        "gas_p50": results.get("gas", {}).get("p50"),
        "oil_p50": results.get("oil", {}).get("p50"),
        "water_p50": results.get("water", {}).get("p50")
    }

def run_batch(folder_path):
    csv_files = glob.glob(os.path.join(folder_path, "*.csv"))
    if not csv_files:
        print(f"No CSV files found in {folder_path}")
        return
    
    print(f"Found {len(csv_files)} wells in {folder_path}")
    print("=" * 60)
    
    summary = []
    for filepath in sorted(csv_files):
        try:
            result = run_single_well(filepath, save_csv=True, save_plot=True, quiet=True)
            summary.append(result)
            gas = f"{result['gas_p50']:.1f}" if result['gas_p50'] else "-"
            oil = f"{result['oil_p50']:.1f}" if result['oil_p50'] else "-"
            water = f"{result['water_p50']:.1f}" if result['water_p50'] else "-"
            print(f"  {result['file']:30s}  Gas={gas:>8}  Oil={oil:>8}  Water={water:>8}")
        except Exception as e:
            print(f"  {os.path.basename(filepath):30s}  ERROR: {e}")
    
    print("=" * 60)
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv("batch_summary.csv", index=False)
    print(f"\nSaved: batch_summary.csv")
    
    totals = {}
    if summary_df["gas_p50"].notna().any():
        totals["gas"] = summary_df["gas_p50"].sum()
    if summary_df["oil_p50"].notna().any():
        totals["oil"] = summary_df["oil_p50"].sum()
    if summary_df["water_p50"].notna().any():
        totals["water"] = summary_df["water_p50"].sum()
    
    print("\nPORTFOLIO TOTALS (P50):")
    if "gas" in totals:
        print(f"  Gas:   {totals['gas']:>10.1f} MMSCF")
    if "oil" in totals:
        print(f"  Oil:   {totals['oil']:>10.1f} MBO")
    if "water" in totals:
        print(f"  Water: {totals['water']:>10.1f} MBW")

def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  Single well:  python3 run_example.py <well.csv>")
        print("  Batch mode:   python3 run_example.py --batch <folder>")
        sys.exit(1)
    
    if sys.argv[1] == "--batch":
        if len(sys.argv) < 3:
            print("Error: Specify folder path")
            sys.exit(1)
        run_batch(sys.argv[2])
    else:
        run_single_well(sys.argv[1])

if __name__ == "__main__":
    main()