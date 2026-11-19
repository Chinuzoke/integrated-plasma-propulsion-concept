"""
absorption.py -- inverse-bremsstrahlung (free-free) absorption by electron-ion collisions.

Classical/Kramers hydrogenic result with stimulated-emission correction (Rybicki & Lightman, Radiative
Processes in Astrophysics, 1979, eq. 5.19b), written here in SI from fundamentals:

    alpha = C_ff * g_ff * (sum_s Z_s^2 n_s) * n_e * T^(-1/2) * nu^(-3) * (1 - exp(-h nu / kT))      [1/m]
    C_ff  = (4 q^6 / (3 m_e h c)) * sqrt(2 pi / (3 k m_e)),   q^2 -> e^2/(4 pi eps0)  (cgs -> SI)

Units: n in m^-3, T in K, nu in Hz  ->  alpha in m^-1.  The local deposition is q_laser = alpha I [W/m^3].
`validation.py` checks C_ff against the legacy cgs constant 3.692e8 (the two must agree after unit conversion).

Validity / limits (important):
* Electron-ION IB only. Electron-NEUTRAL IB, which dominates at low ionisation and therefore controls the
  ignition threshold, is NOT included: it needs the e-Ar momentum-transfer cross-section, which is an
  UNRESOLVED input here (no sourced data). `electron_neutral_alpha()` raises until one is supplied.
* Gaunt factor g_ff is an adjustable ASSUMPTION (default 1.0, as in the legacy model); it differs from 1 at
  1 um / 1e4 K. Treat absorption as accurate to a factor of order unity at best.
* Requires omega >> collision frequency and omega > plasma frequency (checked by `plasma_frequency_ok`).
* Strong-coupling / Debye-screening corrections ignored.
"""
from __future__ import annotations

import math

import numpy as np

from . import constants as C

_Q2 = C.e_charge**2 / (4.0 * math.pi * C.eps0)          # J m
C_FF_SI = (4.0 * _Q2**3 / (3.0 * C.m_e * C.h_planck * C.c_light)) * math.sqrt(2.0 * math.pi / (3.0 * C.k_B * C.m_e))


def frequency(wavelength_m: float) -> float:
    if not (wavelength_m > 0 and math.isfinite(wavelength_m)):
        raise ValueError("wavelength must be finite and > 0")
    return C.c_light / wavelength_m


def alpha_ff(T, n_e, zsq_n_ion, wavelength_m, gaunt=1.0):
    """Electron-ion inverse-bremsstrahlung absorption coefficient [1/m]."""
    T = np.asarray(T, dtype=float)
    nu = frequency(wavelength_m)
    stim = -np.expm1(-C.h_planck * nu / (C.k_B * T))
    return C_FF_SI * gaunt * zsq_n_ion * n_e * T**-0.5 * nu**-3 * stim


def zsq_n_ion(comp):
    """sum_s Z_s^2 n_s = n_1 + 4 n_2  [m^-3]."""
    return comp["n_1"] + 4.0 * comp["n_2"]


def plasma_frequency(n_e):
    return np.sqrt(np.asarray(n_e, dtype=float) * C.e_charge**2 / (C.eps0 * C.m_e)) / (2.0 * math.pi)


def plasma_frequency_ok(n_e_max, wavelength_m) -> bool:
    """True if the laser frequency exceeds the electron plasma frequency (beam can propagate, not cut off)."""
    return bool(frequency(wavelength_m) > plasma_frequency(n_e_max))


def electron_neutral_alpha(*_args, **_kw):
    raise NotImplementedError("electron-neutral inverse bremsstrahlung needs a sourced e-Ar momentum-transfer "
                              "cross-section: UNRESOLVED input (see lsp/references.py).")
