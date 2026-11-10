"""
scenarios.py -- named alternative cases and sensitivity sweeps.

The reference scenario (baseline_v1) stays reproducible: it uses ONLY the
guaranteed specifications and the unmodified reference inputs. Alternatives
are labelled and never replace it.

The M^2 and pointing sweeps are SENSITIVITY CASES: neither parameter has an
established source, so no sweep value is a prediction or a specification.
"""
from __future__ import annotations

import config as cfg
import pointing_error as pe
from power_budget import compute_budget

# name -> (description, overrides applied to SystemInputs)
SCENARIOS = {
    "baseline_v1": ("Reference scenario: guaranteed manufacturer specifications (mirror 0.998, conditioning 0.960, "
                    "window 0.996), unverified M^2 and pointing as in config.PARAMS", {}),
    "typical_mirror": ("Mirror at TYPICAL measured performance (R = 0.9998, combined abs+scatter <= 150 ppm, "
                       "PROVISIONAL manufacturer measurement); other components unchanged",
                       dict(mirror_R=cfg.value("mirror_typical_reflectance"),
                            mirror_abs_scatter_ppm=cfg.value("mirror_typical_abs_scatter_ppm"))),
    "uncoated_window_comparison": ("COMPARISON ONLY: uncoated fused-silica Fresnel window (the representative "
                                   "Shimadzu window is AR-coated); bulk absorption/scatter not included",
                                   dict(window_mode="uncoated_fresnel")),
}

POINTING_CASES_ARCSEC = (0.0, 0.01, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0)   # UNVERIFIED sensitivity cases
M2_CASES = (1.0, 3.0, 5.0, 11.8, 20.0)                                  # UNVERIFIED sensitivity cases


def run_scenarios(base: cfg.SystemInputs) -> list:
    """[(name, description, budget)]"""
    return [(n, d, compute_budget(base.with_(**ov))) for n, (d, ov) in SCENARIOS.items()]


def pointing_sweep(base: cfg.SystemInputs, cases=POINTING_CASES_ARCSEC) -> list:
    """[(arcsec, displacement_m, capture_efficiency, delivered_W)]"""
    rows = []
    for a in cases:
        b = compute_budget(base.with_(pointing_arcsec=a))
        rows.append((a, pe.displacement_from_arcsec(base.distance_m, a), b.capture.efficiency, b.delivered_w))
    return rows


def m2_sweep(base: cfg.SystemInputs, cases=M2_CASES) -> list:
    """[(M2, beam_radius_m, capture_efficiency, delivered_W)]"""
    rows = []
    for m in cases:
        b = compute_budget(base.with_(m2=m))
        rows.append((m, b.beam.beam_radius_m, b.capture.efficiency, b.delivered_w))
    return rows
