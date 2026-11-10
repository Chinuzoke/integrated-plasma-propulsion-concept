"""
uncertainty.py -- sensitivity and illustrative uncertainty analysis.

THREE different things, kept separate:

1. LOCAL ELASTICITY  S_x = d ln(P_del) / d ln(x)  (central finite difference).
   Mathematical sensitivity of the model at the baseline. Not an uncertainty.
2. ONE-AT-A-TIME RANGE SWEEP (tornado) over ILLUSTRATIVE ranges.
3. MONTE CARLO over the same ILLUSTRATIVE ranges (uniform, seeded).

HONESTY NOTE: no measured uncertainty data exist for these inputs in this
project. The ranges below are ILLUSTRATIVE and are NOT derived from
instrument specifications or literature. They show which parameters the
output is most sensitive to under plausible-sized variations; they are not
confidence intervals.

Parameter uncertainty (an input is imperfectly known) is covered here.
MODEL-FORM uncertainty (Gaussian-profile assumption for an M^2 = 11.8 beam,
single-number M^2, fixed-offset pointing, ignored jitter/turbulence/
aberrations/obscuration, incoherent window model) is NOT quantified by this
module; see README 'Limitations'.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr

from config import SystemInputs
from power_budget import compute_budget

PARAM_LABELS = {
    "power_w": "Transmitter power", "wavelength_m": "Wavelength", "waist_m": "Beam waist radius",
    "m2": "M^2", "distance_m": "Propagation distance", "aperture_radius_m": "Receiver aperture radius",
    "pointing_arcsec": "Pointing error", "mirror_R": "Mirror reflectance",
    "cond_T": "Beam-conditioning transmission", "window_T": "Window transmission",
}

# (low, high) in the parameter's own SI/user unit. ILLUSTRATIVE, see module docstring.
ILLUSTRATIVE_RANGES = {
    "power_w": (27_000.0, 33_000.0),          # +/-10 %
    "wavelength_m": (1063e-9, 1065e-9),       # +/-1 nm
    "waist_m": (0.225, 0.275),                # +/-10 %
    "m2": (10.62, 12.98),                     # +/-10 %
    "distance_m": (1.08e6, 1.32e6),           # +/-10 %
    "aperture_radius_m": (4.5, 5.5),          # +/-10 %
    "pointing_arcsec": (0.05, 0.30),
    "mirror_R": (0.998, 0.9998),    # guaranteed spec .. typical measured (manufacturer values)
    "cond_T": (0.93, 0.98),
    "window_T": (0.990, 0.998),
}


def delivered_power(inp: SystemInputs) -> float:
    return compute_budget(inp).delivered_w


def local_elasticities(inp: SystemInputs, rel_step: float = 1e-3) -> dict:
    """S_x = d ln P / d ln x by central differences (keeps fractions <= 1)."""
    out = {}
    for name in PARAM_LABELS:
        x0 = getattr(inp, name)
        hi, lo = x0 * (1 + rel_step), x0 * (1 - rel_step)
        p_hi = delivered_power(inp.with_(**{name: hi}))
        p_lo = delivered_power(inp.with_(**{name: lo}))
        out[name] = float((np.log(p_hi) - np.log(p_lo)) / (np.log(hi) - np.log(lo)))
    return out


def tornado(inp: SystemInputs, ranges: dict | None = None) -> list:
    """Rows (name, low_value, high_value, P_at_low, P_at_high, swing_W), sorted by swing."""
    ranges = ranges or ILLUSTRATIVE_RANGES
    rows = []
    for name, (lo, hi) in ranges.items():
        p_lo = delivered_power(inp.with_(**{name: lo}))
        p_hi = delivered_power(inp.with_(**{name: hi}))
        rows.append((name, lo, hi, p_lo, p_hi, abs(p_hi - p_lo)))
    return sorted(rows, key=lambda r: r[5], reverse=True)


def monte_carlo(inp: SystemInputs, n: int = 2000, seed: int = 12345, ranges: dict | None = None) -> dict:
    """Uniform sampling over the illustrative ranges. Returns stats and Spearman rank correlations."""
    ranges = ranges or ILLUSTRATIVE_RANGES
    rng = np.random.default_rng(seed)
    names = list(ranges)
    samples = {k: rng.uniform(*ranges[k], size=n) for k in names}
    p = np.array([delivered_power(inp.with_(**{k: float(samples[k][i]) for k in names})) for i in range(n)])
    corr = {k: float(spearmanr(samples[k], p).statistic) for k in names}
    return {"n": n, "seed": seed, "mean": float(p.mean()), "std": float(p.std(ddof=1)),
            "p05": float(np.percentile(p, 5)), "p50": float(np.percentile(p, 50)),
            "p95": float(np.percentile(p, 95)), "min": float(p.min()), "max": float(p.max()),
            "spearman": dict(sorted(corr.items(), key=lambda kv: -abs(kv[1])))}
