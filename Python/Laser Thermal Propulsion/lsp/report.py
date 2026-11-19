"""report.py -- LSP and integrated text reports; CSV exports."""
from __future__ import annotations

import csv

import numpy as np

from . import config as lc
from .validation import Check

LINE = "=" * 100


def _acc_lines(a: dict) -> list:
    return [
        f"   Laser power entering the chamber (P_in)        = {a['P_in']:12.3f} W",
        f"   Absorbed by the plasma (P_abs)                 = {a['P_abs']:12.3f} W   (absorption efficiency {100 * a['absorption_efficiency']:.3f} %)",
        f"   Residual laser beam leaving (P_res)            = {a['P_res']:12.3f} W",
        f"   beam balance P_in - P_abs - P_res              = {a['beam_balance_residual_W']:12.3e} W",
        f"   Optically thin bremsstrahlung loss (P_rad)     = {a['P_rad']:12.3f} W",
        f"   Radial/wall loss, closure-dependent (P_wall)   = {a['P_wall']:12.3f} W",
        f"   Gas enthalpy-flux increase mdot (h_out - h_in) = {a['dH']:12.3f} W",
        f"   Plasma energy-balance residual                 = {a['energy_balance_residual_W']:12.3e} W  ({a['energy_balance_residual_rel']:.2e} of P_in)",
    ]


def format_lsp_report(res, checks: list | None = None, title="LSP SUBSYSTEM REPORT") -> str:
    cfg, a, pr = res.cfg, res.accounting, res.profiles
    L = [LINE, title, LINE, ""]
    L.append("SCOPE: 1-D axial, steady-state, LTE argon plasma heated by an inverse-bremsstrahlung-absorbed laser beam. Ends at the LSP")
    L.append("outlet state. No nozzle, thrust, specific-impulse or trajectory calculation.")
    L.append("")
    L.append("1. INPUTS AND PROVENANCE")
    for name, val, unit, prov, note in lc.parameter_register(cfg):
        L.append(f"   {name:<42}{str(val):>16} {unit:<9}[{prov.name}]")
        if note:
            L.append(f"        - {note}")
    L.append("")
    L.append("2. SOLUTION STATUS")
    for k, r in res.branches.items():
        L.append(f"   seed '{k}': converged={r.converged}  iterations={r.n_iter}  rel. residual={r.residual:.2e}  "
                 f"energy-balance residual={r.balance_residual:.2e}  Tpeak={r.T.max():.0f} K  -> {('ignited' if r.converged and r.T.max() > 6000 else 'cold' if r.converged else 'NOT CONVERGED')}")
    L.append(f"   PRIMARY (reported) result: seed '{res.primary_kind}' -> {res.classification.upper()}")
    L.append(f"   convergence reason: {res.primary.reason}")
    L.append("   (convergence is judged on the governing-equation residual and global energy balance, not on peak temperature)")
    L.append("")
    L.append("3. LASER POWER AND ENERGY ACCOUNTING (laser optical balance is separate from the plasma thermal balance)")
    L += _acc_lines(a)
    L.append("")
    L.append("4. PLASMA STATE")
    L.append(f"   Peak plasma temperature                        = {a['peak_T']:12.1f} K   (outlet {a['outlet_T']:.1f} K)")
    L.append(f"   Peak electron density                          = {a['peak_n_e']:12.4e} m^-3")
    L.append(f"   Peak ionisation fraction (n_e / n_heavy)       = {a['peak_ion_fraction']:12.4f}")
    L.append(f"   Axial profile (every {max(1, len(pr['x']) // 10)}th cell):")
    L.append(f"      {'x [mm]':>8}{'T [K]':>10}{'n_e [m^-3]':>13}{'ion frac':>10}{'alpha [1/m]':>13}{'I/I0':>9}{'q_laser [W/m^3]':>17}")
    I0 = pr["I"][0]
    step = max(1, len(pr["x"]) // 10)
    for j in range(0, len(pr["x"]), step):
        L.append(f"      {1e3 * pr['x'][j]:8.2f}{pr['T'][j]:10.1f}{pr['n_e'][j]:13.3e}{pr['ion_fraction'][j]:10.4f}{pr['alpha'][j]:13.3e}{pr['I'][j] / I0:9.4f}{pr['q_laser'][j]:17.3e}")
    L.append("")
    L.append("5. HAND-OFF STATE (for a future nozzle stage; nothing downstream is computed)")
    o = res.output
    L.append(f"   mdot = {o.mdot_kg_s:.3e} kg/s, p = {o.pressure_pa:.3e} Pa, T_out = {o.T_out_K:.1f} K, h_out = {o.h_out_J_kg:.4e} J/kg "
             f"(inlet {o.h_in_J_kg:.4e}), ionisation fraction at outlet {o.ion_fraction_out:.4f}")
    L.append(f"   gas enthalpy-flux increase = {o.gas_enthalpy_rise_W:.3f} W; converged={o.converged}; ignited={o.ignited}")
    L.append(f"   NOTE: {o.note}")
    L.append("")
    L.append("6. WARNINGS AND KNOWN LIMITATIONS FOR THIS RESULT")
    for w in res.warnings:
        L.append(f"   - {w}")
    if checks is not None:
        L.append("")
        L.append("7. VERIFICATION STATUS")
        n_ok = sum(c.passed for c in checks)
        L.append(f"   {n_ok}/{len(checks)} verification checks passed")
        for cat in ("math", "consistency", "literature", "audit", "inputs"):
            sub = [c for c in checks if c.category == cat]
            if sub:
                L.append(f"   [{cat}] {sum(c.passed for c in sub)}/{len(sub)}")
        for c in checks:
            if not c.passed:
                L.append(f"   FAILED: {c.name} {c.detail}")
        L.append("   Verified (numerical/analytic): equation implementation, Jacobian, convergence, conservation.")
        L.append("   NOT validated experimentally: LTE closure, 1-D collimated geometry, transport/radiation closures, absorption (e-ion only).")
    L.append(LINE)
    return "\n".join(L)


def format_integrated_report(ir, lsp_checks=None) -> str:
    b = ir.budget
    L = [LINE, "INTEGRATED LASER-DELIVERY + LSP REPORT", LINE, ""]
    L.append("The delivery subsystem decides how much optical power reaches the LSP chamber; the LSP subsystem decides what")
    L.append("the plasma does with it. Delivered power is NOT assumed to be absorbed.")
    L.append("")
    L.append("A. POWER CHAIN")
    L.append(f"   {'stage':<58}{'power [W]':>14}  subsystem")
    for label, p, sub in ir.chain:
        L.append(f"   {label:<58}{p:14.3f}  {sub}")
    L.append("")
    L.append(f"   Laser-delivery efficiency (delivered / transmitter peak)  = {100 * b.efficiency:.4f} %")
    a = ir.lsp.accounting
    L.append(f"   Plasma absorption efficiency (absorbed / delivered)       = {100 * a['absorption_efficiency']:.4f} %")
    L.append(f"   End-to-end absorbed / transmitter peak                    = {100 * a['P_abs'] / b.inputs.power_w:.4f} %")
    L.append(f"   Laser power basis: {ir.laser.power_basis}")
    L.append("")
    L.append("B. INTERFACE CHECKS")
    for c in ir.checks:
        L.append(f"   [{'PASS' if c.passed else 'FAIL'}] {c.name}" + (f"  ({c.detail})" if c.detail else ""))
    L.append("")
    L.append("C. LSP SUBSYSTEM RESULT")
    L.append(format_lsp_report(ir.lsp, lsp_checks, title="LSP SUBSYSTEM REPORT (integrated run)"))
    return "\n".join(L)


def write_lsp_csv(res, path_budget: str, path_profiles: str) -> None:
    a = res.accounting
    with open(path_budget, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value", "unit"])
        for k, u in (("P_in", "W"), ("P_abs", "W"), ("P_res", "W"), ("P_rad", "W"), ("P_wall", "W"), ("dH", "W"),
                     ("absorption_efficiency", "-"), ("peak_T", "K"), ("outlet_T", "K"), ("peak_n_e", "m^-3"),
                     ("peak_ion_fraction", "-"), ("energy_balance_residual_W", "W"), ("beam_balance_residual_W", "W")):
            w.writerow([k, a[k], u])
        w.writerow(["classification", res.classification, ""])
        w.writerow(["converged", res.primary.converged, ""])
    p = res.profiles
    keys = ["x", "T", "n_e", "ion_fraction", "n_neutral", "n_ion", "rho", "alpha", "I", "q_laser", "q_rad", "q_wall", "k", "h"]
    n = len(p["x"])
    with open(path_profiles, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for j in range(n):
            w.writerow([p[k][j] for k in keys])


def write_register_csv(cfg, path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["parameter", "value", "unit", "provenance", "note"])
        for name, val, unit, prov, note in lc.parameter_register(cfg):
            w.writerow([name, val, unit, prov.name, note])
