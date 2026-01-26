# The Waller Decomposition

**Multi-Resolution Wavelet EUR Forecasting for Unconventional Wells**

**Technical White Paper • Version 3.0**

**Author:** Eric Waller  
**Date:** January 2026  
**Status:** Research  
**License:** MIT

---

## Abstract

This paper presents the Waller Decomposition, a multi-resolution wavelet analytical framework for Estimated Ultimate Recovery (EUR) forecasting in multi-stage fractured horizontal wells. The method combines discrete wavelet transforms on Rate-Normalized Pressure (RNP) signals with physics-based drainage distance mapping.

The framework:

1. Decomposes production signals into multi-scale components via discrete wavelet transform
2. Maps wavelet scales to drainage distances using hydraulic diffusivity relationships
3. Fits decline models to each component independently
4. Reconstructs composite EUR from component contributions
5. Quantifies uncertainty via Bayesian parameter estimation
6. Corrects for parent-child well interference using pressure depletion models

This paper presents the theoretical derivation, implementation details, and limitations of the method.

---

## 1. Introduction

### 1.1 The EUR Forecasting Problem

Accurate EUR forecasting in unconventional shale reservoirs remains one of the most consequential challenges in petroleum engineering. Investment decisions worth billions of dollars—drilling programs, reserve bookings, asset acquisitions—depend on EUR estimates made with limited production data.

Traditional Decline Curve Analysis (DCA) methods face fundamental limitations:

| Method | Assumption | Reality in MFHW |
|--------|------------|-----------------|
| Arps (1945) | Single drainage volume | 20-60 independent SRVs |
| Arps | Boundary-dominated flow | Extended linear flow |
| Duong (2011) | Empirical fit | No physical basis |
| Stretched Exponential | Mathematical convenience | No drainage physics |

### 1.2 The Waller Decomposition Hypothesis

The Waller Decomposition is founded on a central hypothesis:

> Wavelet decomposition of RNP signals separates the superposed drainage contributions from individual fracture stages, enabling physics-based extrapolation of each component independently.

Each wavelet scale corresponds to a characteristic drainage distance governed by hydraulic diffusivity. By identifying and extrapolating each scale's contribution, the Waller Decomposition reconstructs the composite EUR from N stage-level forecasts rather than fitting a single empirical curve to the aggregate signal.

---

## 2. Theoretical Foundation

### 2.1 Rate-Normalized Pressure (RNP)

The fundamental diagnostic variable is Rate-Normalized Pressure:

```
RNP(t) = Δp(t) / q(t) = (p_i - p_wf(t)) / q(t)
```

Where:
- p_i = initial reservoir pressure (psi)
- p_wf(t) = flowing bottomhole pressure (psi)
- q(t) = production rate (MSCF/d or bbl/d)

RNP has units of psi/(MSCF/d) and is directly proportional to the pressure transient solution for a producing well.

### 2.2 Discrete Wavelet Transform

The Discrete Wavelet Transform (DWT) decomposes a signal into approximation and detail coefficients at multiple scales:

```
DWT: x[n] → {cA_J, cD_J, cD_{J-1}, ..., cD_1}
```

Where:
- cA_J = approximation coefficients at level J (low-frequency content)
- cD_j = detail coefficients at level j (high-frequency content at scale 2^j)

The Waller Decomposition uses the Daubechies-4 (db4) wavelet as the default basis due to its balance of smoothness and compact support.

### 2.3 Scale-to-Drainage Distance Derivation

This is the key theoretical contribution of the Waller Decomposition. We derive the relationship between wavelet scale and drainage distance from first principles.

**Step 1: Hydraulic Diffusivity**

The pressure diffusivity equation for single-phase flow in a homogeneous reservoir:

```
∂p/∂t = η ∇²p
```

Where hydraulic diffusivity η = k / (φ μ c_t)

**Step 2: Radius of Investigation**

The radius of investigation at time t:

```
r_inv(t) = √(4 η t) = √(4 k t / (φ μ c_t))
```

**Step 3: Wavelet Scale to Time**

For DWT at sampling interval Δt, the characteristic time for level j:

```
τ_j = 2^j × Δt
```

**Step 4: Scale-to-Distance Mapping**

Combining the above:

```
r_j = C_eff × √(2^j × Δt)
```

Where C_eff = √(4η) is the effective drainage constant (ft/√day).

### 2.4 Basin-Specific Calibration Constants

| Basin | Formation | C_eff (ft/√day) | Notes |
|-------|-----------|-----------------|-------|
| Permian | Wolfcamp A | 25-35 | Oil window |
| Permian | Bone Spring | 22-30 | Variable |
| Eagle Ford | Upper | 20-28 | Oil/condensate |
| Eagle Ford | Lower | 25-32 | Dry gas |
| Marcellus | Lower | 35-45 | High GOR |
| Bakken | Middle | 15-22 | Tight matrix |
| Haynesville | — | 40-55 | Overpressured |
| Niobrara | B Bench | 20-28 | Variable quality |

### 2.5 When Scale-to-Stage Mapping Holds

The Waller Decomposition framework assumes wavelet components correspond to drainage from distinct fracture stages or stage groups. This assumption has validity conditions.

**When Assumptions Are Violated**

Waller Decomposition components represent mathematical frequency bands, not physical stages. The method may still provide useful forecasts by capturing multi-timescale behavior, but the physical interpretation (scale = drainage distance) is weakened.

**Recommendation:** For wells with known interference or parent-child effects, interpret Waller Decomposition results as empirical multi-scale fits rather than physics-based stage decomposition.

---

## 3. Implementation

The Waller Decomposition system is implemented as a modular Python framework compatible with standard oilfield data formats and common data providers (Enverus, DrillingInfo, IHS).

### 3.1 Input Requirements

The Waller Decomposition requires three time-series inputs:

| Input | Units | Source |
|-------|-------|--------|
| Production Rate q(t) | MSCF/d or bbl/d | Production allocation |
| Flowing BHP p_wf(t) | psi | Gauges or calculated |
| Initial Pressure p_i | psi | DFIT, MDT, or analog |

### 3.2 Input Validation

```python
import numpy as np

def validate_inputs(t, q, p_wf, p_i):
    """Validate production data for Waller Decomposition analysis."""
    errors = []
    
    if len(t) < 30:
        errors.append(f"Insufficient data: {len(t)} points (minimum 30)")
    
    if np.any(q <= 0):
        errors.append(f"Non-positive rates detected: {np.sum(q <= 0)} points")
    
    if np.any(p_wf <= 0):
        errors.append(f"Non-positive BHP detected: {np.sum(p_wf <= 0)} points")
    
    if p_i <= np.max(p_wf):
        errors.append(f"Initial pressure ({p_i}) must exceed max BHP ({np.max(p_wf)})")
    
    gaps = np.diff(t)
    if np.max(gaps) > 30:
        errors.append(f"Large data gap detected: {np.max(gaps):.0f} days")
    
    return {'valid': len(errors) == 0, 'errors': errors}
```

### 3.3 Data Preprocessing

```python
import numpy as np
from scipy.interpolate import CubicSpline

def preprocess_production(t_days, q, p_wf, is_cumulative=False, min_points=50):
    """Preprocess production data for Waller Decomposition analysis."""
    if is_cumulative:
        q = np.diff(q) / np.diff(t_days)
        t_days = t_days[1:]
        p_wf = p_wf[1:]
    
    q = np.maximum(q, 1.0)
    
    if len(t_days) < min_points:
        t_interp = np.linspace(t_days[0], t_days[-1], min_points)
        q_spline = CubicSpline(t_days, q)
        p_spline = CubicSpline(t_days, p_wf)
        q = np.maximum(q_spline(t_interp), 1.0)
        p_wf = p_spline(t_interp)
        t_days = t_interp
    
    return t_days, q, p_wf
```

### 3.4 Component Separation

```python
import pywt
import numpy as np
from scipy.optimize import curve_fit

def power_law_decline(t, qi, D, n):
    """Generalized power-law decline: q = qi * (1 + D*t)^(-n)"""
    return qi * (1 + D * t) ** (-n)

def reconstruct_dwt_components(signal, wavelet='db4', level=None):
    """Decompose signal via DWT and reconstruct each level independently."""
    coeffs = pywt.wavedec(signal, wavelet, level=level)
    components = []
    
    for i in range(len(coeffs)):
        coeffs_zeroed = [np.zeros_like(c) for c in coeffs]
        coeffs_zeroed[i] = coeffs[i].copy()
        component = pywt.waverec(coeffs_zeroed, wavelet)
        components.append(component[:len(signal)])
    
    return components

def separate_rate_components(t, q, p_wf, p_i, wavelet='db4', level=5):
    """Separate production rate into multi-scale components via RNP decomposition."""
    delta_p = p_i - p_wf
    rnp = delta_p / q
    rnp_log = np.log10(np.maximum(rnp, 1e-6))
    
    rnp_components_log = reconstruct_dwt_components(rnp_log, wavelet, level)
    rnp_components = [10 ** comp for comp in rnp_components_log]
    
    inv_rnp = [1.0 / np.maximum(comp, 1e-6) for comp in rnp_components]
    inv_rnp_sum = np.sum(inv_rnp, axis=0)
    
    q_components = [q * (inv_comp / inv_rnp_sum) for inv_comp in inv_rnp]
    
    q_reconstructed = np.sum(q_components, axis=0)
    reconstruction_error = np.mean(np.abs(q_reconstructed - q) / q) * 100
    
    return q_components, reconstruction_error
```

### 3.5 Decline Fitting with Uncertainty

```python
def fit_components(t, q_components, min_contribution=0.01):
    """Fit power-law decline to each rate component."""
    q_total_mean = np.mean(np.sum(q_components, axis=0))
    stage_params = []
    
    for k, q_k in enumerate(q_components):
        if np.mean(q_k) < min_contribution * q_total_mean:
            continue
        
        try:
            popt, pcov = curve_fit(
                power_law_decline,
                t, q_k,
                p0=[q_k[0], 0.01, 0.5],
                bounds=([1e-6, 1e-6, 0.01], [q_k[0] * 10, 1.0, 2.0]),
                maxfev=10000
            )
            qi, D, n = popt
            perr = np.sqrt(np.diag(pcov))
            
            stage_params.append({
                'level': k,
                'qi': qi, 'D': D, 'n': n,
                'qi_std': perr[0], 'D_std': perr[1], 'n_std': perr[2]
            })
        except:
            continue
    
    return stage_params
```

### 3.6 EUR Calculation

```python
from scipy.integrate import quad

def calculate_eur(stage_params, t_abandon=10950, q_min=10):
    """Calculate EUR for each component and total well."""
    eur_by_stage = []
    
    for stage in stage_params:
        qi, D, n = stage['qi'], stage['D'], stage['n']
        
        if n > 0 and D > 0:
            t_econ = ((qi / q_min) ** (1/n) - 1) / D
            t_end = min(t_abandon, t_econ)
        else:
            t_end = t_abandon
        
        if 0 < n < 1:
            eur_stage = qi / (D * (1 - n)) * ((1 + D * t_end) ** (1 - n) - 1)
        else:
            eur_stage, _ = quad(lambda t: power_law_decline(t, qi, D, n), 0, t_end)
        
        eur_by_stage.append(eur_stage / 1000)
    
    return sum(eur_by_stage), eur_by_stage
```

### 3.7 Parent-Child Correction

For child wells drilled near existing producers, the initial reservoir pressure is reduced by drainage from parent wells.

**Applicability:** Apply when:
- Child well within 1,500 ft of parent well
- Parent produced > 6 months before child spud
- Same target formation

**Depletion Model:**

Parent drainage radius:
```
r_inv = (C/2) × √t_p
```

Depletion factor:
```
f_dep = 1 - (d_pc/r_inv)² if d_pc < r_inv, else 0
```

EUR correction:
```
EUR_child = EUR_Waller × (1 - α × f_dep)
```

Where α = basin interference factor (typically 0.3-0.5).

**Interference Factor by Basin:**

| Basin | Formation | α |
|-------|-----------|---|
| Permian (Midland) | Wolfcamp A | 0.40-0.50 |
| Permian (Delaware) | Wolfcamp | 0.35-0.45 |
| Eagle Ford | Lower | 0.45-0.55 |
| Marcellus | Lower | 0.30-0.40 |
| Bakken | Middle | 0.35-0.45 |
| Haynesville | — | 0.40-0.50 |

### 3.8 Frac Hit Detection

Frac hits occur when hydraulic fracturing in a child well creates pressure communication with a parent well.

**Detection Method:**

Analyze parent well production in a ±30 day window around child frac dates. Flag anomalies exceeding:
- Rate change > 20%
- Pressure change > 100 psi

**Handling:**

| Detection | Action |
|-----------|--------|
| Rate increase >10% | Include data |
| Rate decrease >10% | Exclude 60 days post-hit |
| No significant change | No action |

---

## 4. Uncertainty Quantification

### 4.1 Sources of Uncertainty

| Source | Type | How Addressed |
|--------|------|---------------|
| Measurement noise | Aleatory | Propagated via MCMC |
| Pressure data quality | Epistemic | Quality weighting |
| Wavelet choice | Epistemic | Sensitivity analysis |
| Decline model form | Epistemic | Model comparison |
| Calibration constant C | Epistemic | Basin-specific priors |

### 4.2 Bayesian Parameter Estimation

Instead of point estimates from least-squares fitting, the Waller Decomposition uses Markov Chain Monte Carlo (MCMC) to sample the full posterior distribution of decline parameters.

### 4.3 EUR Probability Distribution

EUR uncertainty is computed by propagating parameter samples through the EUR integral.

```python
import numpy as np
import emcee
from scipy.integrate import quad

def log_likelihood(theta, t, q_obs, q_sigma):
    """Gaussian log-likelihood for decline curve fit."""
    qi, D, n = theta
    if qi <= 0 or D <= 0 or n <= 0:
        return -np.inf
    q_model = qi * (1 + D * t) ** (-n)
    residuals = (q_obs - q_model) / q_sigma
    return -0.5 * np.sum(residuals ** 2)

def log_prior(theta, qi_prior_mean=5000, n_prior_mean=0.8):
    """Weakly informative priors based on basin analogs."""
    qi, D, n = theta
    
    if not (1 < qi < 100000 and 1e-4 < D < 1.0 and 0.1 < n < 2.5):
        return -np.inf
    
    lp_qi = -0.5 * ((np.log(qi) - np.log(qi_prior_mean)) / 1.0) ** 2
    lp_n = -0.5 * ((n - n_prior_mean) / 0.5) ** 2
    lp_D = 0.0
    
    return lp_qi + lp_n + lp_D

def log_probability(theta, t, q_obs, q_sigma):
    """Posterior log-probability."""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, t, q_obs, q_sigma)

def fit_decline_bayesian(t, q, n_samples=5000):
    """Fit decline curve with full Bayesian uncertainty."""
    q_sigma = np.maximum(0.05 * q, 10)
    
    n_walkers = 32
    p0 = np.array([q[0], 0.01, 0.5])
    pos = p0 + 0.1 * p0 * np.random.randn(n_walkers, 3)
    pos = np.abs(pos)
    
    sampler = emcee.EnsembleSampler(
        n_walkers, 3, log_probability,
        args=(t, q, q_sigma)
    )
    sampler.run_mcmc(pos, n_samples, progress=False)
    
    samples = sampler.get_chain(discard=1000, thin=10, flat=True)
    
    results = {}
    for i, name in enumerate(['qi', 'D', 'n']):
        results[name] = {
            'P10': np.percentile(samples[:, i], 10),
            'P50': np.percentile(samples[:, i], 50),
            'P90': np.percentile(samples[:, i], 90),
            'mean': np.mean(samples[:, i]),
            'std': np.std(samples[:, i])
        }
    results['samples'] = samples
    
    return results

def compute_eur_distribution(samples, t_end=10950, q_min=10):
    """Compute EUR distribution from MCMC parameter samples."""
    eur_samples = []
    
    for qi, D, n in samples:
        if D > 0 and n > 0:
            t_econ = ((qi / q_min) ** (1/n) - 1) / D if n != 0 else np.inf
            t_calc = min(t_end, t_econ)
        else:
            t_calc = t_end
        
        if 0 < n < 1:
            eur = qi / (D * (1 - n)) * ((1 + D * t_calc) ** (1 - n) - 1)
        else:
            eur, _ = quad(lambda t: qi * (1 + D * t) ** (-n), 0, t_calc)
        
        eur_samples.append(eur / 1000)
    
    return {
        'P10': np.percentile(eur_samples, 10),
        'P50': np.percentile(eur_samples, 50),
        'P90': np.percentile(eur_samples, 90),
        'mean': np.mean(eur_samples),
        'std': np.std(eur_samples)
    }
```

---

## 5. Limitations

### 5.1 Applicable Conditions

The Waller Decomposition applies to:
- Multi-stage fractured horizontal wells
- Unconventional reservoirs (shale, tight)
- 30-180 days of production data
- Wells with pressure data (BHP)
- Single-phase flow

### 5.2 Inapplicable Conditions

| Condition | Reason |
|-----------|--------|
| <30 days data | Insufficient for wavelet decomposition |
| No pressure data | Cannot compute RNP |
| Conventional wells | No multi-stage superposition |
| Dense simultaneous development | Interference exceeds correction capability |
| Severe liquid loading | Non-monotonic RNP |
| Multi-phase segregation | Violates single-phase assumption |

### 5.3 Data Requirements

| Parameter | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| Production data | 30 days | 90+ days | Daily or monthly |
| Pressure data | Surface tubing | FBHP gauges | Higher quality = better results |
| Initial pressure | Estimated | DFIT/MDT | Critical for RNP calculation |
| Basin constant C | Literature value | Local calibration | From analog wells |

---

## 6. Conclusion

The Waller Decomposition decomposes multi-stage fractured horizontal well production into components at different drainage timescales using discrete wavelet transform. Each component is mapped to a drainage distance via hydraulic diffusivity, fitted with a decline model, and integrated for EUR.

Uncertainty is quantified via Bayesian parameter estimation, yielding P10/P50/P90 EUR distributions. The method requires bottomhole pressure data and is limited to conditions where the underlying assumptions hold, as specified in Section 5.

---

## Appendix A: Dependencies

Python implementation requires:

```
numpy>=1.21, scipy>=1.7, pywt>=1.3, matplotlib>=3.5, emcee>=3.1, pandas>=1.4
```

---

## Appendix B: Basin Calibration Constants

The effective calibration constant C_eff must be calibrated per basin. Starting values based on typical reservoir properties:

| Basin | Formation | C_eff (ft/√day) | Notes |
|-------|-----------|-----------------|-------|
| Permian | Wolfcamp A | 25-35 | Oil window |
| Permian | Bone Spring | 22-30 | Variable |
| Eagle Ford | Upper | 20-28 | Oil/condensate |
| Eagle Ford | Lower | 25-32 | Dry gas |
| Marcellus | Lower | 35-45 | High GOR |
| Bakken | Middle | 15-22 | Tight matrix |
| Haynesville | — | 40-55 | Overpressured |
| Niobrara | B Bench | 20-28 | Variable quality |

**Calibration procedure:** Select 5-10 mature wells (> 5 years) with known EUR, run Waller Decomposition with range of C values, select C that minimizes hindcast error, apply to new wells.

---

## Appendix C: Wavelet Basis Sensitivity Analysis

A common concern with wavelet-based methods is sensitivity to the choice of wavelet basis. To address this, the Waller Decomposition was tested using eight wavelet families on synthetic well data (180 days, 3-stage decline, 5500 psi initial pressure).

### C.1 Results

| Wavelet | EUR (MMSCF) | Δ from Baseline |
|---------|-------------|-----------------|
| Daubechies-4 | 7,100.7 | baseline |
| Daubechies-6 | 6,992.6 | −1.52% |
| Daubechies-8 | 7,213.6 | +1.59% |
| Symlet-4 | 6,949.8 | −2.13% |
| Symlet-6 | 6,966.5 | −1.89% |
| Symlet-8 | 7,011.1 | −1.26% |
| Coiflet-2 | 6,977.4 | −1.74% |
| Coiflet-4 | 7,226.5 | +1.77% |

### C.2 Summary Statistics

- **Mean EUR:** 7,054.8 MMSCF
- **Standard Deviation:** 104.5 MMSCF
- **Coefficient of Variation:** 1.48%

### C.3 Conclusion

EUR estimates vary less than ±2.2% across all tested wavelet families. The coefficient of variation (1.48%) confirms that the Waller Decomposition is **not sensitive to wavelet basis selection**. Daubechies-4 is recommended as the default due to its balance of smoothness and compact support, but results are robust to alternative choices.

---

## Synthetic Validation

Preliminary synthetic validation was performed using a 3-stage power-law decline model (180 days observed, 3% noise). Results showed:

- **True EUR:** 5,547 MMSCF
- **Mean Estimated EUR:** 6,127 MMSCF (100 Monte Carlo realizations)
- **Mean Error:** +10.5%
- **P10–P90 Range:** 5,507–6,767 MMSCF

The true EUR falls within the P10–P90 uncertainty bounds, but the method exhibits a positive bias (~10% overestimation). Further calibration of the wavelet scale-to-drainage distance mapping is needed to reduce systematic error. This remains an active area of development. See Appendix D for full methodology.

---

## Appendix D: Synthetic Validation Methodology

To assess the accuracy of the Waller Decomposition, a synthetic well was constructed with known EUR for comparison.

### D.1 Synthetic Well Design

A 3-stage fractured horizontal well was simulated using superposed power-law decline:

| Stage | Description | q_i (MSCF/d) | D (1/day) | n | EUR (MMSCF) |
|-------|-------------|--------------|-----------|---|-------------|
| 1 | Near-wellbore | 2,000 | 0.015 | 0.7 | 1,612.7 |
| 2 | Mid-field | 1,000 | 0.008 | 0.5 | 2,103.2 |
| 3 | Far-field | 500 | 0.004 | 0.4 | 1,831.2 |
| **Total** | | | | | **5,547.0** |

Each stage follows: q_k(t) = q_{i,k}(1 + D_k t)^{-n_k}

EUR was calculated analytically by integrating to economic limit (10 MSCF/d) or 30 years, whichever came first.

### D.2 Test Conditions

- **Observed period:** 180 days
- **Noise:** 3% Gaussian on rate
- **Initial pressure:** 5,500 psi
- **BHP:** Declining from 2,500 psi with 15 psi noise

### D.3 Monte Carlo Results (100 Realizations)

| Metric | Value |
|--------|-------|
| True EUR | 5,547 MMSCF |
| Mean Estimated EUR | 6,127 MMSCF |
| Standard Deviation | 503 MMSCF |
| P10 | 5,507 MMSCF |
| P50 | 6,150 MMSCF |
| P90 | 6,767 MMSCF |
| Mean Error | +10.5% |
| Error Std Dev | 9.1% |

### D.4 Interpretation

The Waller Decomposition shows a systematic positive bias of approximately 10% on this synthetic case. However:

1. **The true EUR (5,547 MMSCF) falls within the P10–P90 range**, indicating the uncertainty quantification captures the true value.

2. **The bias is consistent**, suggesting a calibration offset rather than random error.

3. **Possible causes:**
   - Wavelet scale-to-distance mapping assumes uniform diffusivity
   - DWT octave spacing may not perfectly match physical drainage scales
   - 180 days may be insufficient to constrain far-field behavior

### D.5 Conclusion

The synthetic validation demonstrates that the Waller Decomposition provides EUR estimates within ~10% of true values, with uncertainty bounds that encompass the true EUR. Reducing the systematic bias through improved scale calibration is identified as a priority for future development.

---

## References

1. Arps, J.J. (1945). Analysis of Decline Curves. *Trans. AIME*, 160(1), 228-247.
2. Lee, W.J. (1982). Well Testing. *SPE Textbook Series*, Vol. 1.
3. Duong, A.N. (2011). Rate-Decline Analysis for Fracture-Dominated Shale Reservoirs. *SPE REE*, 14(3), 377-387.
4. Valko, P.P., & Lee, W.J. (2010). A Better Way to Forecast Production from Unconventional Gas Wells. SPE 134231.
5. Foreman-Mackey, D. (2013). emcee: The MCMC Hammer. *PASP*, 125(925), 306.

---

**© 2026 Eric Waller. Released under MIT License.**

- Website: [ewaller.com](https://ewaller.com)
- Company: [luxiedge.com](https://luxiedge.com)
