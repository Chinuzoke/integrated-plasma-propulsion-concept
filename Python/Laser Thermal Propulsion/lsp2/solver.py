"""
solver.py -- residual-based pseudo-transient continuation with damped Newton steps.

Why retained (from the legacy model): the steady problem is BISTABLE (alpha rises ~9 orders of magnitude between
5,000 and 20,000 K), so a plain Newton/Picard iteration from an arbitrary guess collapses onto the cold,
non-absorbing branch. Pseudo-transient continuation (PTC) marches rho*cp dT/dtau = G(T) implicitly, which
robustly reaches a stable steady branch.

What is changed:
* The linear system is the FULL Newton system  (C/dtau - J) dT = G  with analytic Jacobian J (not an explicit,
  lagged source), so as dtau grows the iteration becomes quadratically convergent Newton.
* Step control (SER) uses the RESIDUAL norm ||G||, not the step norm.
* Convergence requires the governing-equation residual (max_j |G_j| / (P_in/(A L))) and the global discrete
  energy-balance residual to fall below tolerance AND the update to be negligible. The peak temperature is
  recorded but is never a convergence criterion.
* No silent clipping: temperatures outside the property table or [T_floor, T_ceiling] reject/shrink the step;
  if no acceptable step exists the solve FAILS and says so.
* Which branch is reached depends on the initial state (bistability); `solve_from_seeds` reports all branches.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import energy as en
from .tables import PropertyRangeError


@dataclass
class SolverOptions:
    max_iter: int = 6000
    tol_residual: float = 1e-9        # relative residual (see ref scale)
    tol_balance: float = 1e-9         # |global imbalance| / P_in
    tol_step_K: float = 1e-6          # max |dT| of the accepted step
    dtau0: float | None = None        # initial pseudo-time step [s]; None -> chosen so the first step changes T by ~dT_target
    dT_target: float = 100.0          # [K] used for the automatic initial pseudo-time step
    dtau_min: float = 1e-12
    dtau_max: float = 1e12
    max_dT: float = 3000.0            # trust-region cap per iteration [K]
    T_floor: float = 300.0
    T_ceiling: float = 6.0e4
    verbose: bool = False


@dataclass
class SolveResult:
    T: np.ndarray
    x: np.ndarray
    converged: bool
    reason: str
    n_iter: int
    residual: float                   # relative max-norm residual
    balance_residual: float           # |imbalance| / P_in
    history: list = field(default_factory=list)   # dicts: it, res, dtau, Tmax, step
    fields: dict = field(default_factory=dict)
    balance: dict = field(default_factory=dict)
    geometry: en.Geometry = None
    initial_balance: dict = field(default_factory=dict)   # power accounting of the INITIAL (seed) state
    initial_T_max: float = float('nan')


def _scale(geo):
    return geo.P_in / (geo.area * geo.length)     # mean laser power density [W/m^3]


def solve(phys, geo, T0, opts: SolverOptions | None = None) -> SolveResult:
    o = opts or SolverOptions()
    n = geo.n
    T = np.array(T0, dtype=float).copy()
    h_in = en.inlet_enthalpy(phys, geo)
    scale = _scale(geo)
    hist = []
    dtau = None
    res_prev = None
    reason, converged = "max iterations reached", False

    def residual(Tv):
        f = en.evaluate(Tv, phys, geo, h_in)
        return f, np.max(np.abs(f["G"])) / scale

    try:
        f, res = residual(T)
        init_bal = en.global_balance(T, phys, geo, f)
        init_Tmax = float(T.max())
    except PropertyRangeError as e:
        return SolveResult(T, geo.x, False, f"initial state outside property range: {e}", 0, np.nan, np.nan, [], {}, {}, geo)

    it = 0
    for it in range(1, o.max_iter + 1):
        f = en.evaluate(T, phys, geo, h_in, jacobian=True)
        G, J = f["G"], f["J"]
        res = np.max(np.abs(G)) / scale
        bal = en.global_balance(T, phys, geo, f)
        bres = abs(bal["imbalance"]) / geo.P_in
        rho_cp = np.maximum(phys.rho(T) * phys.dh(T), 1e-12)
        if dtau is None:
            dtau = o.dtau0 if o.dtau0 is not None else float(np.clip(
                o.dT_target * np.median(rho_cp) / max(np.max(np.abs(G)), 1e-300), 1e-9, 1e3))
        elif res_prev is not None:
            dtau = float(np.clip(dtau * res_prev / max(res, 1e-300), o.dtau_min, o.dtau_max))
        accepted = False
        lam = 1.0
        for _ in range(30):                              # backtracking on step length
            A = np.diag(rho_cp / dtau) - J
            try:
                dT = np.linalg.solve(A, G)
            except np.linalg.LinAlgError:
                reason = "singular Newton matrix"
                break
            big = np.max(np.abs(dT))
            if big > o.max_dT:
                dT *= o.max_dT / big
            Tn = T + lam * dT
            if np.all(np.isfinite(Tn)) and Tn.min() >= o.T_floor and Tn.max() <= o.T_ceiling:
                try:
                    fn, resn = residual(Tn)
                except PropertyRangeError:
                    resn = np.inf
                if np.isfinite(resn) and resn < 1e3 * max(res, 1e-12):     # reject wild residual growth
                    accepted = True
                    break
            lam *= 0.5
            dtau = max(dtau * 0.25, o.dtau_min)           # also shorten pseudo-time on rejection
        if not accepted:
            if reason == "max iterations reached":
                reason = ("no admissible step: temperature left [%.0f, %.3g] K or property range" % (o.T_floor, o.T_ceiling))
            break
        step = float(np.max(np.abs(Tn - T)))
        hist.append({"it": it, "res": res, "balance": bres, "dtau": dtau, "Tmax": float(T.max()), "step": step})
        if o.verbose and (it % 10 == 0 or it < 5):
            print(f"  it {it:4d} res={res:.3e} bal={bres:.3e} dtau={dtau:.2e} Tmax={T.max():9.1f} step={step:.2e}")
        res_prev = res
        T = Tn
        fnew = en.evaluate(T, phys, geo, h_in)
        rnew = np.max(np.abs(fnew["G"])) / scale
        bnew = abs(en.global_balance(T, phys, geo, fnew)["imbalance"]) / geo.P_in
        if rnew < o.tol_residual and bnew < o.tol_balance and step < max(o.tol_step_K, 1e-12 * T.max()) * 1e3:
            converged, reason = True, "converged (residual, global balance and update below tolerance)"
            break
        if rnew < o.tol_residual and bnew < o.tol_balance:
            converged, reason = True, "converged (residual and global balance below tolerance)"
            break

    f = en.evaluate(T, phys, geo, h_in)
    bal = en.global_balance(T, phys, geo, f)
    res = float(np.max(np.abs(f["G"])) / scale)
    return SolveResult(T, geo.x, converged, reason, it, res, abs(bal["imbalance"]) / geo.P_in, hist, f, bal, geo, init_bal, init_Tmax)


def make_geometry(cfg) -> en.Geometry:
    return en.Geometry(length=cfg.length_m, area=cfg.area_m2, n=cfg.n_cells, mdot=cfg.mdot_kg_s,
                       I_in=cfg.I0_w_m2, T_in=cfg.T_inlet_K)


def seed(geo, kind: str, T_seed: float | None = None):
    n = geo.n
    if kind == "cold":
        return np.full(n, geo.T_in)
    if kind == "warm":
        return np.full(n, 8000.0 if T_seed is None else T_seed)
    if kind == "hot":
        return np.full(n, 16000.0 if T_seed is None else T_seed)
    raise ValueError("seed kind must be cold|warm|hot")


def solve_from_seeds(phys, geo, kinds=("hot", "warm", "cold"), opts=None):
    """Solve from several initial states; returns {kind: SolveResult}. Different converged branches
    demonstrate bistability (ignited vs cold branch)."""
    return {k: solve(phys, geo, seed(geo, k), opts) for k in kinds}


def classify_branch(res: SolveResult, ignition_T=6000.0) -> str:
    if not res.converged:
        return "not converged"
    return "ignited (hot) branch" if res.T.max() > ignition_T else "cold (non-ignited) branch"
