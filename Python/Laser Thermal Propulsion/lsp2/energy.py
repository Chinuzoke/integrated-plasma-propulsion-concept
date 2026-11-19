"""
energy.py -- finite-volume steady energy equation + Beer-Lambert beam, with analytic Jacobian.

Cell-centred finite volumes, N cells of width dx = L/N, uniform cross-section A, uniform pressure p (low Mach),
mass flux G = mdot/A constant. For cell j (T_j = cell temperature) the NET HEATING residual [W/m^3] is

    G_j = - G (h_j - h_{j-1})/dx                      axial convection of ENTHALPY, 1st-order upwind
          + (F_{j+1/2} - F_{j-1/2})/dx                axial conduction,  F = k_face (T_{j+1}-T_j)/dx
          + q_laser,j - q_rad,j - q_wall,j + q_ext    laser deposition, radiation, radial loss

with h_{-1} = h_in (inlet enthalpy; total-flux / Danckwerts inlet: the specified enthalpy flux mdot*h_in
equals convective + diffusive flux at the inlet face, so no extra diffusive flux crosses face -1/2) and zero
conductive flux through the outlet face (dT/dx = 0). Steady state  <=>  G_j = 0 for all j.

Laser (exact exponential per cell, face intensities I_j, I_0 = I_in):
    I_{j+1} = I_j exp(-alpha_j dx),  q_laser,j = I_j (1 - exp(-alpha_j dx))/dx = (I_j - I_{j+1})/dx   (-dI/dx averaged over the cell;
    evaluated with expm1 so that very small absorption is not lost to cancellation)
so that sum_j q_laser,j dx A = P_in - P_out EXACTLY (discrete energy conservation of the beam).
The global discrete balance, exact for the discretisation, is
    sum_j G_j A dx = P_abs - P_rad - P_wall + P_ext - mdot (h_{N-1} - h_in)

Jacobian dG/dT is assembled analytically (tridiagonal part from convection/conduction/local sources plus a
strictly lower-triangular part from the beam's upstream dependence); tests compare it with finite differences.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Geometry:
    length: float
    area: float
    n: int
    mdot: float
    I_in: float
    T_in: float

    @property
    def dx(self):
        return self.length / self.n

    @property
    def x(self):
        return (np.arange(self.n) + 0.5) * self.dx

    @property
    def P_in(self):
        return self.I_in * self.area


def inlet_enthalpy(phys, geo):
    return float(np.asarray(phys.h(np.array([geo.T_in])))[0])


def evaluate(T, phys, geo, h_in=None, jacobian=False):
    """Return dict with G (residual), fields, and optionally J = dG/dT (dense)."""
    n, dx = geo.n, geo.dx
    Gf = geo.mdot / geo.area
    if h_in is None:
        h_in = inlet_enthalpy(phys, geo)
    h = phys.h(T)
    k = phys.k(T)
    kf = 0.5 * (k[:-1] + k[1:])
    alpha = phys.alpha(T)
    I = np.empty(n + 1)
    I[0] = geo.I_in
    I[1:] = geo.I_in * np.exp(-np.cumsum(alpha * dx))
    # Per-cell absorbed intensity I_j (1 - exp(-alpha_j dx)) evaluated with expm1: the difference I_j - I_{j+1} would cancel to
    # EXACTLY 0.0 whenever alpha dx < ~1e-16 (e.g. cold gas, alpha ~ 1e-154 1/m) and erase the (tiny but real) deposition.
    q_l = I[:-1] * (-np.expm1(-alpha * dx)) / dx
    q_r = phys.qrad(T)
    q_w = phys.qwall(T)

    h_up = np.empty(n)
    h_up[0] = h_in
    h_up[1:] = h[:-1]
    conv = -Gf * (h - h_up) / dx
    E = kf * (T[1:] - T[:-1])                      # flux at faces 1/2 .. N-3/2 (N-1 faces)
    Fin = np.concatenate([[0.0], E])               # west face of each cell
    Fout = np.concatenate([E, [0.0]])              # east face of each cell
    diff = (Fout - Fin) / dx**2
    G = conv + diff + q_l - q_r - q_w + phys.q_ext
    out = {"G": G, "h": h, "k": k, "alpha": alpha, "I": I, "q_laser": q_l, "q_rad": q_r, "q_wall": q_w,
           "h_in": h_in, "conv": conv, "diff": diff}
    if jacobian:
        dh, dk = phys.dh(T), phys.dk(T)
        da = phys.dalpha(T)
        J = np.zeros((n, n))
        idx = np.arange(n)
        J[idx, idx] += -Gf * dh / dx
        J[idx[1:], idx[:-1]] += Gf * dh[:-1] / dx
        # conduction
        dE_dTj = 0.5 * dk[:-1] * (T[1:] - T[:-1]) - kf          # dE_j/dT_j      (face j between cells j, j+1)
        dE_dTjp1 = 0.5 * dk[1:] * (T[1:] - T[:-1]) + kf         # dE_j/dT_{j+1}
        j = idx[:-1]
        # east flux of cell j is E_j ; west flux of cell j+1 is E_j
        J[j, j] += dE_dTj / dx**2
        J[j, j + 1] += dE_dTjp1 / dx**2
        J[j + 1, j] -= dE_dTj / dx**2
        J[j + 1, j + 1] -= dE_dTjp1 / dx**2
        # laser: q_j = (I_j - I_{j+1})/dx
        J += np.tril(-dx * q_l[:, None] * da[None, :], -1)
        J[idx, idx] += I[1:] * da
        J[idx, idx] -= phys.dqrad(T) + phys.dqwall(T)
        out["J"] = J
    return out


def global_balance(T, phys, geo, fields=None):
    """Integrated power accounting [W]. Returns dict incl. the discrete-balance residual."""
    f = fields or evaluate(T, phys, geo)
    V = geo.area * geo.dx
    P_abs = float(np.sum(f["q_laser"]) * V)
    P_rad = float(np.sum(f["q_rad"]) * V)
    P_wall = float(np.sum(f["q_wall"]) * V)
    P_ext = float(phys.q_ext * geo.n * V)
    P_res = float(f["I"][-1] * geo.area)
    dH = geo.mdot * (float(f["h"][-1]) - f["h_in"])
    imbalance = P_abs - P_rad - P_wall + P_ext - dH
    return {"P_in": geo.P_in, "P_abs": P_abs, "P_res": P_res, "P_rad": P_rad, "P_wall": P_wall, "P_ext": P_ext,
            "dH": dH, "imbalance": imbalance, "beam_balance": geo.P_in - P_abs - P_res}
