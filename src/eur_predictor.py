import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class EURPredictor:
    def __init__(self):
        self.economic_limit = 50

    def arps_hyperbolic(self, t, qi, di, b):
        return qi / (1 + b * di * t) ** (1 / b)

    def predict(self, df, early_days=45):
        df = df[df['days'] <= early_days].copy()
        # Simple fit
        t = df['days'].values
        q = df['rate_bblpd'].values  # example oil
        popt, _ = curve_fit(self.arps_hyperbolic, t, q, p0=[2000, 0.001, 1.0], bounds=(0, [10000, 1, 2]))
        # Generate forecast
        t_forecast = np.arange(1, 365*30)
        q_forecast = self.arps_hyperbolic(t_forecast, *popt)
        eur = np.trapezoid(q_forecast, t_forecast) / 365  # rough
        print(f'P50 EUR Oil: {eur:.0f} MBO')
        return self

    def save_outputs(self):
        print('Outputs saved (simulated)')