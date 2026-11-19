"""
scan.py -- EXPLORATORY scan of LSP operating inputs at a given laser input.

Purpose: find out whether an ignited steady branch exists at all, and where. NOT a design study and NOT a
performance claim: the operating inputs are not specified by any source in this project, the model is 1-D with
a collimated top-hat beam, and several closures are provisional (see README 'Known limitations').
"""
from __future__ import annotations

import itertools

from .config import LaserInput, LSPConfig
from .run import run_lsp


def exploratory_scan(laser: LaserInput, pressures_pa, mdots_kg_s, radii_m, length_m=0.02, n_cells=100, seeds=("hot",)):
    rows = []
    for R, p, md in itertools.product(radii_m, pressures_pa, mdots_kg_s):
        cfg = LSPConfig(laser=laser, pressure_pa=p, mdot_kg_s=md, plasma_radius_m=R, length_m=length_m, n_cells=n_cells)
        r = run_lsp(cfg, seeds=seeds)
        a = r.accounting
        rows.append({"R_mm": R * 1e3, "p_MPa": p / 1e6, "mdot_mg_s": md * 1e6, "branch": r.classification,
                     "converged": r.primary.converged, "T_peak_K": a["peak_T"], "absorption": a["absorption_efficiency"],
                     "P_abs_W": a["P_abs"], "P_rad_W": a["P_rad"], "P_wall_W": a["P_wall"], "dH_W": a["dH"],
                     "res": r.primary.residual})
    return rows


def format_scan(rows) -> str:
    h = f"{'R[mm]':>6}{'p[MPa]':>8}{'mdot[mg/s]':>11}  {'branch':<26}{'Tpeak[K]':>9}{'abs[%]':>8}{'Prad[W]':>9}{'Pwall[W]':>9}{'dH[W]':>8}"
    L = [h, "-" * len(h)]
    for r in rows:
        L.append(f"{r['R_mm']:6.2f}{r['p_MPa']:8.2f}{r['mdot_mg_s']:11.3f}  {r['branch']:<26}{r['T_peak_K']:9.0f}"
                 f"{100 * r['absorption']:8.2f}{r['P_rad_W']:9.0f}{r['P_wall_W']:9.0f}{r['dH_W']:8.1f}")
    return "\n".join(L)
