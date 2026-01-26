# Waller Decomposition - Quick Start Guide

## What This Does
Forecasts EUR (Estimated Ultimate Recovery) with P10/P50/P90 uncertainty bounds using wavelet decomposition of rate-normalized pressure data.

**Supports:** Gas wells (MSCF/d) and oil wells (bbl/d)

## Requirements
pip install numpy pandas scipy PyWavelets matplotlib

## Step 1: Prepare Your Data

### For Gas Wells
| Column | Description | Units |
|--------|-------------|-------|
| days | Days since first production | days |
| rate_mscfd | Gas production rate | MSCF/day |
| pressure_psi | Flowing tubing pressure or BHP | psi |

### For Oil Wells
| Column | Description | Units |
| days | Days since first production | days |
| rate_bblpd | Oil production rate | bbl/day |
| pressure_psi | Flowing tubing pressure or BHP | psi |

### Data Tips
- Minimum data: 30+ days recommended, 90+ days preferred
- Shut-ins: Remove or interpolate through shut-in periods

### What About Water?
This method forecasts the primary producing phase (gas or oil). Water is not modeled.

## Step 2: Run the Forecast
python3 run_example.py my_well.csv

## Step 3: Read the Output
| Value | Meaning |
|-------|---------|
| P10 | 90% chance EUR exceeds this (conservative) |
| P50 | Best estimate (median) |
| P90 | 10% chance EUR exceeds this (optimistic) |

Output Units:
- Gas wells: MMSCF
- Oil wells: MBO

## Questions?
Open an issue on GitHub or see waller_decomposition_v3.md
