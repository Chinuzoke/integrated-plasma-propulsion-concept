"""
run.py -- run the LSP subsystem for one configuration and package an auditable result.

Strategy (bistability!): solve from a hot and a cold initial state and report every branch. The 'primary'
result is the converged ignited branch if one exists, else the converged cold branch, else the best
non-converged state clearly flagged as NOT CONVERGED. Nothing is hidden.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import absorption as ab
from . import argon_thermo as at
from . import constants as C
from . import diagnostics as dg
from . import energy as en
from . import solver as sv
from . import transport as tr
from .config import LSPConfig
from .physics import build_physics
from .tables import T_MAX


@dataclass
class LSPOutput:
    """Hand-off to a downstream (future) nozzle stage. NO nozzle/thrust quantities are computed here."""
    mdot_kg_s: float
    pressure_pa: float
    T_out_K: float
    h_out_J_kg: float
    h_in_J_kg: float
    gas_enthalpy_rise_W: float
    ion_fraction_out: float
    converged: bool
    ignited: bool
    note: str = "mixed-mean outlet state NOT computed: the 1-D plug model has no radial structure"


@dataclass
class LSPResult:
    cfg: LSPConfig
    branches: dict
    primary_kind: str
    primary: sv.SolveResult
    classification: str
    accounting: dict
    profiles: dict
    output: LSPOutput
    warnings: list = field(default_factory=list)
    phys: object = None
    diagnostics: dict = field(default_factory=dict)


def profiles_from(res: sv.SolveResult, phys, cfg) -> dict:
    T = res.T
    p = {"x": res.x, "T": T, "I": res.fields["I"], "alpha": res.fields["alpha"],
         "q_laser": res.fields["q_laser"], "q_rad": res.fields["q_rad"], "q_wall": res.fields["q_wall"],
         "k": res.fields["k"], "h": res.fields["h"]}
    comp = at.composition(T, cfg.pressure_pa)
    p.update({"n_e": comp["n_e"], "n_ion": comp["n_1"] + comp["n_2"], "n_neutral": comp["n_a"],
              "ion_fraction": comp["n_e"] / comp["n_heavy"], "rho": at.density(comp)})
    return p


def run_lsp(cfg: LSPConfig, seeds=("hot", "cold"), opts: sv.SolverOptions | None = None, phys=None) -> LSPResult:
    phys = phys or build_physics(cfg)
    geo = sv.make_geometry(cfg)
    branches = sv.solve_from_seeds(phys, geo, seeds, opts)
    ignited = [k for k, r in branches.items() if r.converged and r.T.max() > 6000.0]
    cold = [k for k, r in branches.items() if r.converged and r.T.max() <= 6000.0]
    if ignited:
        kind = ignited[0]
    elif cold:
        kind = cold[0]
    else:
        kind = min(branches, key=lambda k: branches[k].residual if np.isfinite(branches[k].residual) else 1e99)
    prim = branches[kind]
    warns = []
    if not prim.converged:
        warns.append(f"NOT CONVERGED: {prim.reason}. The reported state is NOT a steady solution.")
    if prim.T.max() > 6000 and bool(np.any(tr.tstar_extrapolated(prim.T))):
        warns.append("conductivity uses the Neufeld Omega22* fit beyond its stated range (T > ~9,300 K): extrapolation")
    if prim.T.max() > 6000:
        warns.append("LTE (Saha) assumed; one temperature for electrons and heavy particles is NOT validated here")
    comp_max = at.composition(np.array([prim.T.max()]), cfg.pressure_pa)
    if not ab.plasma_frequency_ok(comp_max["n_e"][0], cfg.laser.wavelength_m):
        warns.append("plasma frequency exceeds laser frequency: beam cut-off, IB model invalid")
    bc = cfg.beam_consistency()
    if bc["checked"] and not bc["collimated_ok"]:
        warns.append(f"beam Rayleigh range {bc['rayleigh_range_m']:.3g} m is < 10x chamber length: collimated top-hat assumption inconsistent with M^2")
    if prim.T.max() > 6000 and prim.T[0] > 6000:
        warns.append("the ignited plasma is ATTACHED TO THE INLET BOUNDARY (hottest cell at x = 0): in a uniform collimated beam the intensity is "
                     "highest at the entrance, nothing anchors the plasma inside the chamber, and the solution depends on the inlet (total-flux) "
                     "boundary condition. A focused beam (not modelled) is what normally locates the plasma.")
    warns.append("1-D collimated top-hat beam; no focusing, no radial structure (Jeng & Keefer 1987: 1-D LSP models are inherently limited)")
    warns.append("radiation = optically thin bremsstrahlung only (line/recombination omitted: loss underestimated)")
    warns.append("thermal conductivity omits reactive and electron contributions (conduction underestimated)")
    warns.append("operating inputs (p, mdot, R, L, T_in, T_wall) are ORIGINAL_MODEL_DEFAULT/illustrative, not experimental conditions")

    b = prim.balance
    acc = dict(b)
    acc["absorption_efficiency"] = b["P_abs"] / b["P_in"]
    acc["power_to_gas_enthalpy"] = b["dH"]
    acc["peak_T"] = float(prim.T.max())
    acc["outlet_T"] = float(prim.T[-1])
    acc["energy_balance_residual_W"] = b["imbalance"]
    acc["energy_balance_residual_rel"] = abs(b["imbalance"]) / b["P_in"]
    acc["beam_balance_residual_W"] = b["beam_balance"]
    prof = profiles_from(prim, phys, cfg)
    acc["peak_n_e"] = float(prof["n_e"].max())
    acc["peak_ion_fraction"] = float(prof["ion_fraction"].max())
    class_ = ("ignited (hot) branch" if (prim.converged and prim.T.max() > 6000) else
              "cold (non-ignited) branch" if prim.converged else "not converged")
    out = LSPOutput(cfg.mdot_kg_s, cfg.pressure_pa, float(prim.T[-1]), float(prim.fields["h"][-1]),
                    float(prim.fields["h_in"]), b["dH"], float(prof["ion_fraction"][-1]), prim.converged,
                    prim.converged and prim.T.max() > 6000)
    diag = {"absorption": dg.absorption_summary(prim, cfg), "uniform_column": dg.uniform_column_balance(phys, cfg),
            "seeds": {k: {"initial_T_max_K": r.initial_T_max, "initial_P_abs_W": r.initial_balance.get("P_abs", float("nan")),
                          "initial_P_abs_fraction": r.initial_balance.get("P_abs", float("nan")) / b["P_in"],
                          "final_T_max_K": float(r.T.max()), "final_P_abs_W": r.balance["P_abs"], "converged": r.converged}
                      for k, r in branches.items()}}
    ab_ = diag["absorption"]
    if prim.converged and ab_["absorbed_fraction_from_tau"] < 1e-6:
        warns.append(f"NEGLIGIBLE ABSORPTION is a model RESULT, not a fault: alpha_max = {ab_['alpha_max_1_m']:.3e} 1/m at T = {ab_['T_at_alpha_max_K']:.0f} K, "
                     f"optical depth over the chamber = {ab_['optical_depth_total']:.3e}")
    return LSPResult(cfg, branches, kind, prim, class_, acc, prof, out, warns, phys, diag)
