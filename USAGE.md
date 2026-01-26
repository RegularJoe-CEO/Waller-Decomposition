# Waller Decomposition - Quick Start Guide

## What This Does
Forecasts EUR (Estimated Ultimate Recovery) with P10/P50/P90 uncertainty bounds using wavelet decomposition of rate-normalized pressure data.

## Requirements
```
pip install numpy pandas scipy pywt matplotlib
```

## Step 1: Prepare Your Data

Create a CSV file with these **exact column names**:

| Column | Description | Units |
|--------|-------------|-------|
| `days` | Days since first production | days (integer or decimal) |
| `rate_mscfd` | Gas production rate | MSCF/day |
| `pressure_psi` | Flowing tubing pressure or BHP | psi |

**Example file (`my_well.csv`):**
```
days,rate_mscfd,pressure_psi
1,2450,3850
7,2380,3720
14,2290,3610
30,2050,3440
60,1720,3210
90,1480,3050
120,1290,2920
180,1020,2720
365,610,2370
```

### Data Tips
- **Minimum data:** 30+ days recommended, 90+ days preferred
- **Frequency:** Daily or weekly is fine. Monthly works but reduces accuracy.
- **Gaps:** Small gaps are OK. Large gaps (>30 days) may affect results.
- **Shut-ins:** Remove or interpolate through shut-in periods.

## Step 2: Export from Your System

**From ARIES:**
1. Run production query for your well
2. Export to CSV
3. Rename columns to match: `days`, `rate_mscfd`, `pressure_psi`

**From PHDWin:**
1. Production Data → Export → CSV
2. Rename columns to match

**From Excel/Manual:**
1. Columns A, B, C = days, rate_mscfd, pressure_psi
2. Save As → CSV (Comma delimited)

## Step 3: Run the Forecast

```bash
python3 run_example.py my_well.csv
```

## Step 4: Read the Output

```
========================================
EUR FORECAST (MMSCF)
========================================
  P10:     850.3
  P50:    1070.1
  P90:    1340.7
========================================
```

| Value | Meaning |
|-------|---------|
| P10 | 90% chance EUR exceeds this (conservative) |
| P50 | Best estimate (median) |
| P90 | 10% chance EUR exceeds this (optimistic) |

A plot is saved to `forecast_output.png`.

## Troubleshooting

| Error | Fix |
|-------|-----|
| `Missing column: days` | Check your CSV header row matches exactly |
| `ModuleNotFoundError: pywt` | Run `pip install PyWavelets` |
| Results look wrong | Check units - must be MSCF/day, not MCF or MMSCF |

## Questions?
Open an issue on GitHub or see the full methodology in `waller_decomposition_v3.md`.
