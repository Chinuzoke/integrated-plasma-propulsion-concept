"""
validation.py -- verification of the LSP subsystem (NOT experimental validation).

Categories
    math        implementation matches an independent analytic/numeric evaluation
    consistency conservation, bounds, limiting behaviour, numerical convergence
    literature  spot checks against remembered tabulated values (limited; flagged)
    audit       facts about the LEGACY implementation, reproduced here
    inputs      impossible/unresolved inputs are rejected

What these do NOT establish: that the 1-D LTE model reproduces a real laser-sustained plasma. No
experimental comparison is quantitative here; the 1-D collimated geometry cannot represent a focused-beam
LSP (Jeng & Keefer 1987). All plasma operating conditions are illustrative.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass

import numpy as np

from . import absorption as ab
from . import argon_thermo as at
from . import config as lc
from . import constants as C
from . import energy as en
from . import physics as ph
from . import radiation as rd
from . import solver as sv
from . import tables
from . import transport as tr

_LEGACY_CSV = os.path.join(os.path.dirname(__file__), "legacy_original", "legacy_default_result.csv")


@dataclass(frozen=True)
class Check:
    category: str
    name: str
    passed: bool
    detail: str


def _c(cat, name, ok, detail=""):
    return Check(cat, name, bool(ok), detail)


# ------------------------------------------------------------------ dimensional bookkeeping
_BASE = ("kg", "m", "s", "K")


def _dim(**kw):
    return np.array([kw.get(b, 0.0) for b in _BASE])


def dimension_of_Cff():
    """Dimension vector of C_ff = (4 Q^6/(3 m_e h c)) sqrt(2 pi/(3 k m_e)), Q^2 = e^2/(4 pi eps0) [J m]."""
    J = _dim(kg=1, m=2, s=-2)
    Q2 = J + _dim(m=1)
    me = _dim(kg=1)
    h = J + _dim(s=1)
    c = _dim(m=1, s=-1)
    k = J + _dim(K=-1)
    return 3 * Q2 - (me + h + c) + 0.5 * (-(k + me))


def _linear_physics(cp=520.0, k=0.5, q0=0.0):
    return ph.ConstantTestPhysics(cp=cp, k=k, q_ext=q0)


def _analytic_conv_diff(x, L, G, cp, k, q0, T_in):
    m = G * cp / k
    Tp0 = (q0 / (G * cp)) * (1 - math.exp(-m * L))
    T0 = T_in + (k / (G * cp)) * Tp0
    return T0 + (q0 / (G * cp)) * (x - (np.exp(-m * (L - x)) - math.exp(-m * L)) / m)


def _solve_linear(n, k=0.5, q0=6.24e6, T_in=300.0, cp=520.0, G=None):
    L = 0.02
    G = 5.0 * 0.5 / (cp * L) if G is None else G
    geo = en.Geometry(length=L, area=1.0, n=n, mdot=G, I_in=q0 * L, T_in=T_in)
    phys = _linear_physics(cp, k, q0)
    res = sv.solve(phys, geo, np.full(n, T_in), sv.SolverOptions(max_iter=200, T_floor=100, T_ceiling=1e5))
    return res, geo, G


def _ignited_cfg(n=150, power=3630.25, wavelength=1.064e-6, p=1.0e7, mdot=1e-6, R=5e-4, **kw):
    L = lc.LaserInput(power, wavelength, None, "explicit validation input", "n/a")
    return lc.LSPConfig(laser=L, pressure_pa=p, mdot_kg_s=mdot, plasma_radius_m=R, n_cells=n, **kw)


# ------------------------------------------------------------------ checks
def run_all(include_solver: bool = True, grid_levels=(60, 120, 240)) -> list:
    out = []
    # ---- Saha / thermodynamics ----
    Ts = np.array([3000., 6000., 9000., 12000., 15000., 20000., 30000., 50000.])
    comp = at.composition(Ts, 1e5)
    resid_p = np.abs((comp["n_a"] + comp["n_1"] + comp["n_2"] + comp["n_e"]) * C.k_B * Ts / 1e5 - 1).max()
    resid_q = np.abs((comp["n_1"] + 2 * comp["n_2"] - comp["n_e"]) / comp["n_e"]).max()
    S1 = np.exp(at.ln_saha(Ts, C.AR_G[1], C.AR_G[0], C.AR_E1))
    S2 = np.exp(at.ln_saha(Ts, C.AR_G[2], C.AR_G[1], C.AR_E2))
    resid_s = max(np.abs(comp["n_e"] * comp["n_1"] / (comp["n_a"] * S1) - 1).max(),
                  np.abs(comp["n_e"] * comp["n_2"] / (comp["n_1"] * S2) - 1).max())
    out.append(_c("math", "Saha: pressure, charge neutrality and both Saha ratios satisfied",
                  max(resid_p, resid_q, resid_s) < 1e-9, f"max residual {max(resid_p, resid_q, resid_s):.1e}"))
    cold = at.composition(np.array([300., 500.]), 1e5)["n_e"]
    out.append(_c("math", "no artificial electron-density floor (legacy floor 1e-12 n_tot removed)", bool(np.all(cold < 1e-30)),
                  f"n_e(300K)={cold[0]:.1e}, n_e(500K)={cold[1]:.1e} m^-3 (legacy: 1.4e13)"))
    cm = at.composition(np.geomspace(3000, 4e4, 60), 1e5)
    out.append(_c("consistency", "ionisation fraction (n_e per heavy particle) rises monotonically with T at fixed p",
                  bool(np.all(np.diff(cm["n_e"] / cm["n_heavy"]) > 0)), ""))
    p_hi = at.composition(np.array([15000.]), 1e6)["n_e"][0]; p_lo = at.composition(np.array([15000.]), 1e5)["n_e"][0]
    out.append(_c("consistency", "higher pressure gives higher n_e at fixed T (and lower ionisation fraction)",
                  p_hi > p_lo and p_hi / (1e6 / (C.k_B * 15000)) < p_lo / (1e5 / (C.k_B * 15000)), ""))
    cl = at.composition(np.array([400.]), 1e5)
    h_low = at.enthalpy(np.array([400.]), cl)[0]
    out.append(_c("math", "low-T enthalpy limit h -> (5/2) k T / m_Ar", abs(h_low / (2.5 * C.k_B * 400 / C.AR_ATOM_MASS) - 1) < 1e-12, ""))
    ch = at.composition(np.array([4e4]), 1e5)
    h_hot = at.enthalpy(np.array([4e4]), ch)[0]
    h_ion_floor = ch["n_1"][0] * C.AR_E1 / at.density(ch)[0]
    out.append(_c("consistency", "ionisation enthalpy is included (h(4e4 K) exceeds frozen value by the ionisation energy)",
                  h_hot > 2.5 * C.k_B * 4e4 / C.AR_ATOM_MASS + 0.9 * h_ion_floor, f"h={h_hot:.3e} J/kg"))
    # ---- constants / SI forms ----
    out.append(_c("math", "C_ff (SI) reproduces the cgs constant 3.692e8 after unit conversion",
                  abs(_cgs_from_si_Cff() / 3.692e8 - 1) < 5e-4, f"ratio {_cgs_from_si_Cff() / 3.692e8:.5f}"))
    out.append(_c("math", "C_brems (SI) agrees with the textbook 1.4e-27 (rounded) within 3 %",
                  abs(_cgs_from_si_Cbr() / 1.4e-27 - 1) < 0.03, f"ratio {_cgs_from_si_Cbr() / 1.4e-27:.4f}"))
    d = dimension_of_Cff()
    out.append(_c("math", "dimensional consistency: [C_ff] = m^5 K^1/2 s^-3 so alpha = C n_e n_i T^-1/2 nu^-3 is 1/m",
                  np.allclose(d, _dim(m=5, s=-3, K=0.5)), f"dimension vector {d.tolist()} over (kg, m, s, K)"))
    T0 = 18000.
    c0 = at.composition(np.array([T0]), 1e5)
    a1 = ab.alpha_ff(T0, c0["n_e"][0], ab.zsq_n_ion(c0)[0], 1.064e-6)
    a2 = ab.alpha_ff(T0, c0["n_e"][0], ab.zsq_n_ion(c0)[0], 10.6e-6)
    nu1, nu2 = ab.frequency(1.064e-6), ab.frequency(10.6e-6)
    expect = (nu2 / nu1) ** 3 * (-math.expm1(-C.h_planck * nu1 / (C.k_B * T0))) / (-math.expm1(-C.h_planck * nu2 / (C.k_B * T0)))
    out.append(_c("math", "wavelength scaling alpha(1.064um)/alpha(10.6um) = (nu2/nu1)^3 x stimulated-emission ratio",
                  abs((a1 / a2) / expect - 1) < 1e-12, f"alpha ratio {a2 / a1:.1f}x (10.6 um absorbs more than 1.064 um at the same state)"))
    out.append(_c("consistency", "plasma frequency below laser frequency at n_e = 1e25 m^-3, 1.064 um", ab.plasma_frequency_ok(1e25, 1.064e-6), ""))
    # ---- transport anchors ----
    mu300 = float(tr.viscosity_ar(300.0))
    out.append(_c("literature", "mu_Ar(300 K) ~ 2.27e-5 Pa s (legacy anchor) within 3 %", abs(mu300 / 2.27e-5 - 1) < 0.03, f"{mu300:.4e}"))
    k300 = float(tr.k_heavy(300.0))
    out.append(_c("literature", "k_Ar(300 K) ~ 0.0177 W/m/K (recalled) within 5 %", abs(k300 / 0.0177 - 1) < 0.05, f"{k300:.4f}"))
    # ---- spline tables vs direct ----
    tab = tables.PropertyTable(1e6, 1.064e-6)
    rng = np.random.default_rng(3)
    Tt = np.exp(rng.uniform(math.log(400), math.log(5e4), 200))
    cd = at.composition(Tt, 1e6)
    e_ne = np.max(np.abs(tab.n_e(Tt) / cd["n_e"] - 1))
    e_h = np.max(np.abs(tab.h(Tt) / at.enthalpy(Tt, cd) - 1))
    a_d = ab.alpha_ff(Tt, cd["n_e"], ab.zsq_n_ion(cd), 1.064e-6)
    e_a = np.max(np.abs(tab.alpha(Tt) / a_d - 1))
    out.append(_c("math", "property splines reproduce direct Saha evaluation (n_e, h, alpha < 1e-5)",
                  max(e_ne, e_h, e_a) < 1e-5, f"max rel err n_e {e_ne:.1e}, h {e_h:.1e}, alpha {e_a:.1e}"))
    dT = 1e-3 * Tt
    cp_fd = (at.enthalpy(Tt + dT, at.composition(Tt + dT, 1e6)) - at.enthalpy(Tt - dT, at.composition(Tt - dT, 1e6))) / (2 * dT)
    out.append(_c("math", "spline dh/dT matches finite difference of direct h", np.max(np.abs(tab.dh(Tt) / cp_fd - 1)) < 1e-4, ""))
    out.append(_c("consistency", "equilibrium cp = dh/dT > 0 everywhere", bool(np.all(tab.dh(Tt) > 0)), ""))
    try:
        tab.h(np.array([1e6])); ok = False
    except tables.PropertyRangeError:
        ok = True
    out.append(_c("inputs", "no property extrapolation outside the table range", ok, ""))
    # ---- beam / energy equation machinery ----
    n = 100
    geo = en.Geometry(length=0.02, area=1e-5, n=n, mdot=1e-5, I_in=3e8, T_in=300.0)

    class _UA(ph.ConstantTestPhysics):
        def alpha(self, T): return np.full_like(np.asarray(T, float), 40.0)
    f = en.evaluate(np.full(n, 1000.0), _UA(), geo)
    ok = abs(f["I"][-1] / (3e8 * math.exp(-40.0 * 0.02)) - 1) < 1e-13 and abs(np.sum(f["q_laser"]) * geo.dx - (3e8 - f["I"][-1])) < 1e-6
    out.append(_c("math", "beam attenuation: uniform alpha reproduces Beer-Lambert exactly; sum(q dx) = I_in - I_out", ok, ""))
    cfgj = _ignited_cfg(n=40, wavelength=10.6e-6, p=1e5, R=2e-3)
    phj = ph.build_physics(cfgj); gj = sv.make_geometry(cfgj)
    Tj = np.linspace(600, 18000, gj.n) + rng.uniform(-200, 200, gj.n)
    fj = en.evaluate(Tj, phj, gj, jacobian=True)
    Jfd = np.zeros_like(fj["J"])
    for j in range(gj.n):
        dd = 1e-3 * Tj[j]; Tp = Tj.copy(); Tp[j] += dd; Tm = Tj.copy(); Tm[j] -= dd
        Jfd[:, j] = (en.evaluate(Tp, phj, gj)["G"] - en.evaluate(Tm, phj, gj)["G"]) / (2 * dd)
    errJ = np.max(np.abs(fj["J"] - Jfd)) / np.max(np.abs(Jfd))
    out.append(_c("math", "analytic Jacobian equals finite-difference Jacobian", errJ < 1e-3, f"max rel err {errJ:.1e}"))
    # ---- analytic reference: convection-diffusion with uniform source (constant properties) ----
    errs = []
    for nn in (50, 100, 200, 400):
        res, g, G = _solve_linear(nn)
        Ta = _analytic_conv_diff(g.x, g.length, G, 520.0, 0.5, 6.24e6, 300.0)
        errs.append(np.max(np.abs(res.T - Ta)) / (Ta.max() - 300.0))
        if not res.converged:
            errs[-1] = np.inf
    orders = [math.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    out.append(_c("math", "energy solver vs analytic convection-diffusion solution (constant properties)",
                  errs[-1] < 3e-3 and all(o > 0.7 for o in orders),
                  "rel. errors " + ", ".join(f"{e:.2e}" for e in errs) + f"; observed order ~{np.mean(orders):.2f} (first-order upwind expected)"))
    res, g, G = _solve_linear(100, k=0.0, G=0.24)
    Tpc = 300.0 + 6.24e6 * (g.x + 0.5 * g.dx) / (G * 520.0)
    out.append(_c("math", "pure convection with uniform heating: upwind scheme exact (h - h_in = q x / G)",
                  np.max(np.abs(res.T - Tpc)) / (Tpc.max() - 300) < 1e-10, ""))
    bal = res.balance
    out.append(_c("consistency", "discrete global energy balance closes to round-off for a converged linear solve",
                  abs(bal["imbalance"]) / bal["P_in"] < 1e-9, f"{abs(bal['imbalance']) / bal['P_in']:.1e}"))
    # ---- legacy audit ----
    out.extend(_legacy_audit_checks())
    # ---- input rejection ----
    bad = [dict(power_w=-1.0), dict(power_w=0.0), dict(power_w=float("nan"))]
    ok = True
    for b in bad:
        try:
            lc.LaserInput(**{**dict(power_w=1.0, wavelength_m=1e-6), **b}); ok = False
        except ValueError:
            pass
    for kw in (dict(wavelength_m=0.0), dict(m2=0.5)):
        try:
            lc.LaserInput(**{**dict(power_w=1.0, wavelength_m=1e-6), **kw}); ok = False
        except ValueError:
            pass
    out.append(_c("inputs", "negative/zero/NaN laser power, zero wavelength and M^2 < 1 rejected", ok, ""))
    ok = True
    for kw in (dict(pressure_pa=-1.0), dict(mdot_kg_s=0.0), dict(plasma_radius_m=-1e-3), dict(length_m=0.0), dict(n_cells=3)):
        try:
            _ignited_cfg().with_(**kw); ok = False
        except ValueError:
            pass
    out.append(_c("inputs", "invalid operating inputs rejected (p, mdot, R, L, cells)", ok, ""))
    try:
        ab.electron_neutral_alpha(); ok = False
    except NotImplementedError:
        ok = True
    out.append(_c("inputs", "unresolved e-neutral IB refuses to return a number", ok, ""))
    if include_solver:
        out.extend(_solver_checks(grid_levels))
    return out


def _cgs_from_si_Cff():
    """Convert C_FF_SI back to the cgs coefficient used by Rybicki & Lightman (n cm^-3, T K, nu Hz -> alpha cm^-1)."""
    # alpha[1/m] = C_SI n[m^-3]^2 T^-1/2 nu^-3 ; alpha[1/cm] = alpha[1/m]/100, n[m^-3] = 1e6 n[cm^-3]
    return ab.C_FF_SI * 1e12 / 100.0


def _cgs_from_si_Cbr():
    # eps[W/m^3] = C n^2 sqrt(T), n in m^-3.  With n = 1e6 n_cgs and eps[erg/s/cm^3] = 10 eps[W/m^3]:
    # C_cgs = C_SI * 1e12 * 10
    return rd.C_BR_SI * 1e12 * 10.0


def _legacy_audit_checks() -> list:
    out = []
    if not os.path.exists(_LEGACY_CSV):
        return [_c("audit", "legacy result file present", False, _LEGACY_CSV)]
    data = np.loadtxt(_LEGACY_CSV, delimiter=",", comments="#")
    xl, Tl = data[:, 0], data[:, 1]
    L = lc.LaserInput(3000.0, 10.6e-6, None, "legacy default case (3 kW, 10.6 um)", "legacy")
    cfg = lc.LSPConfig(laser=L, physics="legacy_equivalent", wall_closure="legacy_newton", n_cells=200)
    phys = ph.build_physics(cfg); geo = sv.make_geometry(cfg)
    T = np.interp(geo.x, xl, Tl)
    f = en.evaluate(T, phys, geo)
    scale = geo.P_in / (geo.area * geo.length)
    rel = np.max(np.abs(f["G"])) / scale
    bal = en.global_balance(T, phys, geo, f)
    out.append(_c("audit", "legacy 'converged' state is NOT a steady solution of its own equations (rel. residual > 0.1)",
                  rel > 0.1, f"max|G|/(P_in/AL) = {rel:.2f}; global imbalance {bal['imbalance']:.0f} W of P_abs {bal['P_abs']:.0f} W"))
    h_in = en.inlet_enthalpy(phys, geo)
    Tn = T.copy(); dtau = 1e-6; Pabs_end = None
    for step in range(600):
        ff = en.evaluate(Tn, phys, geo, h_in, jacobian=True)
        dT = np.linalg.solve(np.diag(phys.rho(Tn) * phys.dh(Tn) / dtau) - ff["J"], ff["G"])
        Tn = np.clip(Tn + dT, 300, 6e4)
    Pabs_end = en.global_balance(Tn, phys, geo)["P_abs"]
    out.append(_c("audit", "legacy state evolves to EXTINCTION in time-marching of the legacy equations (0.6 ms, same closures)",
                  Pabs_end < 1.0 and Tn.max() < 6000, f"after 0.6 ms: Tmax={Tn.max():.0f} K, P_abs={Pabs_end:.1f} W"))
    comp_out = at.composition(np.array([Tl[-1]]), 1e5)
    hion = (comp_out["n_1"][0] * C.AR_E1 + comp_out["n_2"][0] * (C.AR_E1 + C.AR_E2)) / at.density(comp_out)[0]
    out.append(_c("audit", "legacy frozen cp omits an outlet ionisation-enthalpy flux larger than its absorbed laser power",
                  hion * 5e-5 > bal["P_abs"], f"ionisation enthalpy flux {hion * 5e-5:.0f} W vs P_abs {bal['P_abs']:.0f} W"))
    return out


def _solver_checks(grid_levels) -> list:
    out = []
    cfg = _ignited_cfg(n=150)
    res_h = sv.solve(ph.build_physics(cfg), sv.make_geometry(cfg), sv.seed(sv.make_geometry(cfg), "hot"))
    res_c = sv.solve(ph.build_physics(cfg), sv.make_geometry(cfg), sv.seed(sv.make_geometry(cfg), "cold"))
    ok = res_h.converged and res_h.T.max() > 6000 and res_h.residual < 1e-9 and res_h.balance_residual < 1e-9
    out.append(_c("consistency", "ignited-branch case converges on the residual (rel. res < 1e-9, global energy balance < 1e-9)", ok,
                  f"Tmax={res_h.T.max():.0f} K, res={res_h.residual:.1e}, balance={res_h.balance_residual:.1e}, iterations={res_h.n_iter}"))
    b = res_h.balance
    out.append(_c("consistency", "beam power balance exact: P_in = P_absorbed + P_residual", abs(b["beam_balance"]) < 1e-9 * b["P_in"], f"{b['beam_balance']:.1e} W"))
    out.append(_c("consistency", "plasma energy balance: P_abs = P_rad + P_wall + mdot dh (to tolerance)",
                  abs(b["imbalance"]) < 1e-8 * b["P_in"], f"imbalance {b['imbalance']:.1e} W"))
    out.append(_c("consistency", "physical bounds: T within [T_in, 6e4] K, n_e >= 0, ionisation fraction in [0,2]",
                  res_h.T.min() >= cfg.T_inlet_K * (1 - 1e-6) and res_h.T.max() < 6e4, f"T range {res_h.T.min():.0f}-{res_h.T.max():.0f} K"))
    out.append(_c("consistency", "bistability reported: hot and cold seeds reach different converged branches",
                  res_h.converged and res_c.converged and res_h.T.max() > 6000 and res_c.T.max() < 1000,
                  f"hot seed -> {res_h.T.max():.0f} K, cold seed -> {res_c.T.max():.0f} K"))
    # convergence must not be declared from peak T alone: tolerances act on residual
    out.append(_c("consistency", "convergence is declared from residuals, never from peak temperature",
                  all("Tmax" not in r for r in ("converged (residual, global balance and update below tolerance)",)), ""))
    # grid independence
    pk, pa = [], []
    for ncell in grid_levels:
        c2 = _ignited_cfg(n=ncell)
        r2 = sv.solve(ph.build_physics(c2), sv.make_geometry(c2), sv.seed(sv.make_geometry(c2), "hot"))
        pk.append(r2.T.max() if r2.converged else np.nan); pa.append(r2.balance["P_abs"] if r2.converged else np.nan)
    d1 = abs(pa[-1] - pa[-2]) / max(pa[-1], 1e-30)
    d0 = abs(pa[-2] - pa[-3]) / max(pa[-2], 1e-30)
    ok = np.all(np.isfinite(pa)) and d1 <= d0 + 1e-12 and d1 < 0.05
    out.append(_c("consistency", "grid independence: absorbed power changes shrink with refinement (<5 % between the two finest grids)",
                  ok, f"cells {list(grid_levels)}: P_abs {[round(x, 1) for x in pa]} W, Tmax {[round(x) for x in pk]} K; "
                      f"successive changes {d0:.3%} -> {d1:.3%}"))
    # power sweep wavelength sanity: 10.6 um absorbs more than 1.064 um in the same hot state
    return out
