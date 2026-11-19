"""
argon_thermo.py -- Saha-equilibrium composition and thermodynamics of argon at fixed pressure.

Model (LTE, ideal-gas mixture of Ar, Ar+, Ar2+, e-), single temperature T:
    n_e n_1 / n_0 = S1(T),   n_e n_2 / n_1 = S2(T)           (Saha, ground-term partition functions)
    S_z(T) = 2 (g_{z+1}/g_z) (2 pi m_e k T / h^2)^{3/2} exp(-E_z / kT)       [m^-3]
    n_e = n_1 + 2 n_2                                          (charge neutrality)
    p = (n_0 + n_1 + n_2 + n_e) k T                            (Dalton / ideal gas)
LTE is an ASSUMPTION, not an experimentally established property of these plasmas: it needs collisional
excitation/ionisation to dominate radiative processes, and T_e = T_heavy. At the low pressures, steep
gradients and laser heating of the electrons it may fail (see README, 'Known limitations').

Solution method: with r1 = S1/n_e, r12 = S1 S2/n_e^2 the mean charge is
    Zbar(n_e) = (r1 + 2 r12) / (1 + r1 + r12)   (monotonically decreasing in n_e)
and the root of  n_e - (p/kT - n_e) Zbar(n_e) = 0  is unique; it is found by bisection in ln n_e,
fully vectorised over T. Unlike the legacy brentq with a floor, no artificial n_e floor is imposed
(numerical lower bracket 1e-200 m^-3).

Thermodynamics (reference: neutral ground-state Ar at 0 K, h = 0):
    rho = n_heavy m_Ar
    h   = [ (5/2) k T (n_heavy + n_e) + n_1 E1 + n_2 (E1 + E2) ] / rho        [J/kg]
which INCLUDES the ionisation enthalpy and the electron translational enthalpy that the legacy frozen
cp omitted. Electronic excitation energy is neglected. cp_eq = dh/dT (equilibrium, includes reactive part).
"""
from __future__ import annotations

import numpy as np

from . import constants as C

_LN_NE_MIN = np.log(1e-200)


def ln_saha(T, g_up, g_low, E):
    """ln S(T) [ln m^-3] (log form avoids exp underflow at low T)."""
    T = np.asarray(T, dtype=float)
    return (np.log(2.0 * g_up / g_low)
            + 1.5 * np.log(2.0 * np.pi * C.m_e * C.k_B * T / C.h_planck**2) - E / (C.k_B * T))


def composition(T, p):
    """Vectorised Saha composition at temperature T [K] (scalar/array) and pressure p [Pa].
    Returns dict of arrays: n_a, n_1, n_2, n_e, n_heavy, n_tot (all m^-3)."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    if np.any(T <= 0) or not np.all(np.isfinite(T)):
        raise ValueError("temperature must be finite and > 0")
    if p <= 0:
        raise ValueError("pressure must be > 0")
    g0, g1, g2 = C.AR_G
    lnS1 = ln_saha(T, g1, g0, C.AR_E1)
    lnS2 = ln_saha(T, g2, g1, C.AR_E2)
    n_tot = p / (C.k_B * T)
    lo = np.full_like(T, _LN_NE_MIN)
    hi = np.log(n_tot)

    def f(y):
        ne = np.exp(y)
        r1 = np.exp(np.clip(lnS1 - y, -700, 700))
        r12 = np.exp(np.clip(lnS1 + lnS2 - 2.0 * y, -700, 700))
        zbar = (r1 + 2.0 * r12) / (1.0 + r1 + r12)
        return ne - (n_tot - ne) * zbar

    for _ in range(80):
        mid = 0.5 * (lo + hi)
        pos = f(mid) > 0.0
        hi = np.where(pos, mid, hi)
        lo = np.where(pos, lo, mid)
    y = 0.5 * (lo + hi)
    n_e = np.exp(y)
    r1 = np.exp(np.clip(lnS1 - y, -700, 700))
    r12 = np.exp(np.clip(lnS1 + lnS2 - 2.0 * y, -700, 700))
    n_h = n_tot - n_e
    n_a = n_h / (1.0 + r1 + r12)
    return {"n_a": n_a, "n_1": n_a * r1, "n_2": n_a * r12, "n_e": n_e, "n_heavy": n_h, "n_tot": n_tot}


def ionisation_fraction(comp):
    """n_e / n_heavy (electrons per heavy particle; = mean charge Zbar)."""
    return comp["n_e"] / comp["n_heavy"]


def density(comp):
    """Mass density [kg/m^3] = n_heavy m_Ar."""
    return comp["n_heavy"] * C.AR_ATOM_MASS


def enthalpy(T, comp):
    """Specific enthalpy [J/kg] including translational (heavy+electron) and ionisation contributions."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    rho = density(comp)
    h_trans = 2.5 * C.k_B * T * (comp["n_heavy"] + comp["n_e"]) / rho
    h_ion = (comp["n_1"] * C.AR_E1 + comp["n_2"] * (C.AR_E1 + C.AR_E2)) / rho
    return h_trans + h_ion


def frozen_cp():
    """Frozen monatomic cp = 5/2 k/m_Ar [J/(kg K)] (legacy value)."""
    return 2.5 * C.k_B / C.AR_ATOM_MASS
