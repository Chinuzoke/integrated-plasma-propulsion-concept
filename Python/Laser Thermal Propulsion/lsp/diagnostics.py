"""
diagnostics.py -- absorption and ignition diagnostics computed from the model itself (nothing assumed).

1. absorption_summary(): the absorption coefficient along the solved state, the optical depth, and an INDEPENDENT recomputation of the
   absorbed power (direct Saha + SI inverse-bremsstrahlung, no property splines, no solver fields) compared with the solver's P_abs.
2. uniform_column_balance(): a necessary-condition view of ignition. For a plasma column held uniformly at temperature T over the chamber
   length, the model's deposition is  P_dep(T) = P_in (1 - exp(-alpha(T) L))  and the unavoidable demands are
   P_loss(T) = mdot (h(T) - h_in) + (q_rad(T) + q_wall(T)) A L.  Where P_dep < P_loss for every T no uniform hot state can be sustained
   (a Raizer-type maintenance criterion). It is a diagnostic, not a replacement for the full solution.
"""
from __future__ import annotations

import numpy as np

from . import absorption as ab
from . import argon_thermo as at
from . import energy as en


def absorption_summary(res, cfg) -> dict:
    """Absorption diagnostics for a solved state; includes an independent recomputation of P_abs."""
    T = res.T
    dx = res.geometry.dx
    alpha = res.fields["alpha"]
    comp = at.composition(T, cfg.pressure_pa)
    alpha_direct = ab.alpha_ff(T, comp["n_e"], ab.zsq_n_ion(comp), cfg.laser.wavelength_m, cfg.gaunt_ff)
    tau = float(np.sum(alpha * dx))
    I = cfg.I0_w_m2 * np.exp(-np.cumsum(alpha_direct * dx))
    P_abs_indep = cfg.laser.power_w * float(-np.expm1(-np.sum(alpha_direct * dx)))
    j = int(np.argmax(alpha))
    return {"alpha_max_1_m": float(alpha.max()), "x_at_alpha_max_mm": float(1e3 * res.x[j]), "T_at_alpha_max_K": float(T[j]),
            "n_e_at_alpha_max": float(comp["n_e"][j]), "ion_fraction_at_alpha_max": float(comp["n_e"][j] / comp["n_heavy"][j]),
            "optical_depth_total": tau, "absorbed_fraction_from_tau": float(-np.expm1(-tau)),
            "P_abs_solver_W": float(res.balance["P_abs"]), "P_abs_independent_W": P_abs_indep,
            "independent_rel_diff": abs(P_abs_indep - res.balance["P_abs"]) / max(abs(P_abs_indep), 1e-300),
            "alpha_spline_vs_direct_rel": float(np.max(np.abs(alpha / np.maximum(alpha_direct, 1e-300) - 1.0))),
            "I0_W_m2": cfg.I0_w_m2, "length_m": cfg.length_m}


def uniform_column_balance(phys, cfg, T=None) -> dict:
    """P_dep(T) and P_loss(T) for a column held uniformly at T [W]."""
    T = np.geomspace(1000.0, 40000.0, 400) if T is None else np.asarray(T, float)
    A, L = cfg.area_m2, cfg.length_m
    h_in = float(phys.h(np.array([cfg.T_inlet_K]))[0])
    alpha = phys.alpha(T)
    P_dep = cfg.laser.power_w * (-np.expm1(-alpha * L))
    P_enth = cfg.mdot_kg_s * (phys.h(T) - h_in)
    P_rad = phys.qrad(T) * A * L
    P_wall = phys.qwall(T) * A * L
    P_loss = P_enth + P_rad + P_wall
    k = int(np.argmax(P_dep - P_loss))
    return {"T": T, "P_dep": P_dep, "P_enthalpy": P_enth, "P_rad": P_rad, "P_wall": P_wall, "P_loss": P_loss,
            "max_surplus_W": float((P_dep - P_loss)[k]), "T_at_max_surplus_K": float(T[k]),
            "max_deposition_W": float(P_dep.max()), "T_at_max_deposition_K": float(T[int(np.argmax(P_dep))]),
            "sustainable_uniform_state_exists": bool(np.any(P_dep >= P_loss))}
