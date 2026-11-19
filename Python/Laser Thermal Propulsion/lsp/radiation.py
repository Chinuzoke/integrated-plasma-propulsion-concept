"""
radiation.py -- optically thin bremsstrahlung emission (retained from the legacy model, rewritten in SI).

    eps_ff = C_br * g_B * (sum_s Z_s^2 n_s) * n_e * T^(1/2)                    [W/m^3]
    C_br   = (32 pi q^6 / (3 h m_e c^3)) * sqrt(2 pi k / (3 m_e)),   q^2 = e^2/(4 pi eps0)
(Rybicki & Lightman eq. 5.15b: 1.4e-27 sqrt(T) n_e n_i Z^2 g_B erg s^-1 cm^-3, g_B ~ 1.2, which is the
legacy constant. validation.py verifies C_br against 1.4e-27.)

OMITTED (documented limitation, NOT modelled): recombination continuum and line radiation, and any
reabsorption. For argon near 1e4-2e4 K these dominate the emission by a large factor (use a net-emission-
coefficient table; the interface `NECTable` accepts a user-supplied sourced table). Treating the plasma as
optically thin brems-only therefore UNDERESTIMATES radiative loss. The legacy grey-body T^4 'calibrated'
term is NOT used in the reconstructed model: it had no scientific basis (it was added to bound the solution).
"""
from __future__ import annotations

import csv
import math

import numpy as np

from . import constants as C

_Q2 = C.e_charge**2 / (4.0 * math.pi * C.eps0)
C_BR_SI = (32.0 * math.pi * _Q2**3 / (3.0 * C.h_planck * C.m_e * C.c_light**3)) * math.sqrt(2.0 * math.pi * C.k_B / (3.0 * C.m_e))


def bremsstrahlung_loss(T, n_e, zsq_n_ion, gaunt_b=1.2):
    """Volumetric optically thin free-free emission [W/m^3]."""
    return C_BR_SI * gaunt_b * zsq_n_ion * n_e * np.sqrt(np.asarray(T, dtype=float))


class NECTable:
    """User-supplied net emission coefficient table: CSV with columns T_K, nec_W_m3_sr at a stated pressure.
    Volumetric loss = 4 pi * NEC. NO table is bundled (no sourced data available); provenance is the user's."""

    def __init__(self, T_K, nec_W_m3_sr, pressure_pa, source):
        self.T = np.asarray(T_K, float)
        self.nec = np.asarray(nec_W_m3_sr, float)
        self.p = pressure_pa
        self.source = source
        if np.any(np.diff(self.T) <= 0):
            raise ValueError("NEC table temperatures must be increasing")

    @classmethod
    def from_csv(cls, path, pressure_pa, source):
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f))
        return cls([float(r["T_K"]) for r in rows], [float(r["nec_W_m3_sr"]) for r in rows], pressure_pa, source)

    def loss(self, T):
        T = np.asarray(T, float)
        if np.any(T < self.T[0]) or np.any(T > self.T[-1]):
            raise ValueError("temperature outside NEC table range (no extrapolation)")
        return 4.0 * math.pi * np.interp(T, self.T, self.nec)
