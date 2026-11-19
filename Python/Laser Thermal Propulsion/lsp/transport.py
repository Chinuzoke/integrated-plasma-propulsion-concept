"""
transport.py -- thermal conductivity closures.  THIS IS THE WEAKEST CLOSURE: read the limits.

Baseline (reconstructed) conductivity:
    k = k_multiplier * k_heavy(T)
    k_heavy = (15/4) (k_B / m_Ar) mu(T)                  monatomic Eucken / Chapman-Enskog 1st approximation
    mu(T)   = 2.6693e-6 sqrt(M[g/mol] T) / (sigma[A]^2 Omega22*(T*)),  T* = T/(eps/k)   [Pa s]
    Omega22* = 1.16145 T*^-0.14874 + 0.52487 exp(-0.77320 T*) + 2.16178 exp(-2.43787 T*)   (Neufeld et al. 1972 fit)
Ar-Ar Lennard-Jones (sigma = 3.542 A, eps/k = 93.3 K, recalled from Hirschfelder-Curtiss-Bird). The check
mu(300 K) ~ 2.27e-5 Pa s (the value the legacy code anchored to) is in validation.py.

LIMITS (all documented, none hidden):
* Neufeld's fit is stated for T* in [0.3, 100] i.e. T <~ 9,300 K for argon; above that it is an EXTRAPOLATION.
* Heavy species (Ar, Ar+, Ar2+) are all treated as LJ argon: ion-atom charge exchange and Coulomb
  interactions are ignored.
* OMITTED: (i) REACTIVE (ionisation-recombination) thermal conductivity and (ii) ELECTRON translational
  conductivity. In real argon plasma these produce a pronounced conductivity peak in the ionisation range
  (Murphy & Arundell 1994; Boulos-Fauchais-Pfender 1994). Their magnitude requires sourced e-Ar and
  ion-atom collision data which are UNRESOLVED here. The baseline k therefore UNDERESTIMATES conduction
  (axial and, through the wall closure, radial) in roughly the 8-20 kK band. `k_multiplier` brackets this
  and `TabulatedConductivity` lets you replace k by a sourced table (e.g. from Murphy & Arundell).
* This is not a validated property set.
"""
from __future__ import annotations

import csv

import numpy as np

from . import constants as C


def omega22_star(Tstar):
    Tstar = np.asarray(Tstar, dtype=float)
    return 1.16145 * Tstar**-0.14874 + 0.52487 * np.exp(-0.77320 * Tstar) + 2.16178 * np.exp(-2.43787 * Tstar)


def viscosity_ar(T):
    """Ar-Ar LJ Chapman-Enskog viscosity [Pa s]."""
    T = np.asarray(T, dtype=float)
    Tstar = T / C.AR_LJ_EPS_OVER_K
    M = C.AR_MOLAR_MASS * 1e3
    return 2.6693e-6 * np.sqrt(M * T) / (C.AR_LJ_SIGMA_ANGSTROM**2 * omega22_star(Tstar))


def k_heavy(T):
    """Heavy-particle translational conductivity [W/(m K)] = (15/4)(k/m) mu."""
    return 3.75 * (C.k_B / C.AR_ATOM_MASS) * viscosity_ar(T)


def tstar_extrapolated(T) -> np.ndarray:
    """Boolean mask: temperatures beyond the stated range of the Neufeld fit (T* > 100)."""
    return np.asarray(T, dtype=float) / C.AR_LJ_EPS_OVER_K > 100.0


# ---- legacy closure (for the legacy-equivalent physics only) ----
def legacy_viscosity(T):
    return 2.27e-5 * (np.asarray(T, dtype=float) / 300.0) ** 0.72


def legacy_k(T):
    gamma = 5.0 / 3.0
    cv = (C.R_gas / C.AR_MOLAR_MASS) / (gamma - 1.0)
    return legacy_viscosity(T) * cv * (9.0 * gamma - 5.0) / 4.0


class TabulatedConductivity:
    """User-supplied total thermal conductivity table (CSV: T_K, k_W_mK) at the operating pressure.
    Nothing is bundled: obtain values from a sourced reference (e.g. Murphy & Arundell 1994, Plasma Chem.
    Plasma Process. 14, 451) and record the source here. No extrapolation."""

    def __init__(self, T_K, k_W_mK, source):
        self.T = np.asarray(T_K, float)
        self.k = np.asarray(k_W_mK, float)
        self.source = source
        if np.any(np.diff(self.T) <= 0) or np.any(self.k <= 0):
            raise ValueError("table needs increasing T and positive k")

    @classmethod
    def from_csv(cls, path, source):
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f))
        return cls([float(r["T_K"]) for r in rows], [float(r["k_W_mK"]) for r in rows], source)

    def __call__(self, T):
        T = np.asarray(T, float)
        if np.any(T < self.T[0]) or np.any(T > self.T[-1]):
            raise ValueError("temperature outside conductivity-table range (no extrapolation)")
        return np.interp(T, self.T, self.k)
