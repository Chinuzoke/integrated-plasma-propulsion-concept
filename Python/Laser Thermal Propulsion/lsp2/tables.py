"""
tables.py -- pre-tabulated property splines at fixed pressure (fast, smooth, differentiable).

At the (uniform) operating pressure every thermodynamic/radiative property is a function of T only. They are
computed once on a fine log-spaced temperature grid by the vectorised Saha solver and interpolated with cubic
splines in ln T (log-values for positive, widely varying quantities). Derivatives come from the spline, which
the Newton Jacobian uses. validation.py compares the splines with direct Saha evaluation.
No extrapolation: an out-of-range temperature raises PropertyRangeError (the solver treats it as a failed step).
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from . import absorption as ab
from . import argon_thermo as at
from . import radiation as rd
from . import transport as tr

T_MIN, T_MAX, N_GRID = 250.0, 8.0e4, 9000
_LN_FLOOR = np.log(1e-300)


class PropertyRangeError(ValueError):
    pass


class _LogSpline:
    """Spline of ln(y) vs ln(T); returns y and dy/dT."""

    def __init__(self, lnT, y):
        self.cs = CubicSpline(lnT, np.log(np.maximum(y, 1e-300)))
        self.dcs = self.cs.derivative()

    def __call__(self, lnT):
        return np.exp(self.cs(lnT))

    def deriv(self, lnT, T):
        return np.exp(self.cs(lnT)) * self.dcs(lnT) / T


class PropertyTable:
    def __init__(self, pressure_pa, wavelength_m, gaunt_ff=1.0, gaunt_brems=1.2, k_multiplier=1.0,
                 k_override=None, nec_table=None):
        self.p = pressure_pa
        self.T = np.geomspace(T_MIN, T_MAX, N_GRID)
        lnT = np.log(self.T)
        comp = at.composition(self.T, pressure_pa)
        zn = ab.zsq_n_ion(comp)
        self.comp_grid = comp
        h = at.enthalpy(self.T, comp)
        self._h = CubicSpline(lnT, h)
        self._dh = self._h.derivative()
        self._rho = _LogSpline(lnT, at.density(comp))
        self._alpha = _LogSpline(lnT, ab.alpha_ff(self.T, comp["n_e"], zn, wavelength_m, gaunt_ff))
        if nec_table is not None:
            qr = nec_table.loss(np.clip(self.T, nec_table.T[0], nec_table.T[-1]))
        else:
            qr = rd.bremsstrahlung_loss(self.T, comp["n_e"], zn, gaunt_brems)
        self._qrad = _LogSpline(lnT, qr)
        k = (k_override(self.T) if k_override is not None else tr.k_heavy(self.T)) * k_multiplier
        self._k = _LogSpline(lnT, k)
        # Kirchhoff conduction potential S0(T) = int_{T_MIN}^{T} k dT (trapezoid on the dense grid)
        S0 = np.concatenate([[0.0], np.cumsum(0.5 * (k[1:] + k[:-1]) * np.diff(self.T))])
        self._S0 = CubicSpline(lnT, S0)
        self._ne = _LogSpline(lnT, comp["n_e"])
        self._ion = CubicSpline(lnT, comp["n_e"] / comp["n_heavy"])

    # ---- guards ----
    def _ln(self, T):
        T = np.asarray(T, dtype=float)
        if not np.all(np.isfinite(T)) or np.any(T < T_MIN) or np.any(T > T_MAX):
            raise PropertyRangeError(f"temperature outside table range [{T_MIN}, {T_MAX}] K")
        return np.log(T), T

    def h(self, T):
        l, _ = self._ln(T); return self._h(l)

    def dh(self, T):
        l, T = self._ln(T); return self._dh(l) / T

    def rho(self, T):
        l, _ = self._ln(T); return self._rho(l)

    def alpha(self, T):
        l, _ = self._ln(T); return self._alpha(l)

    def dalpha(self, T):
        l, T = self._ln(T); return self._alpha.deriv(l, T)

    def qrad(self, T):
        l, _ = self._ln(T); return self._qrad(l)

    def dqrad(self, T):
        l, T = self._ln(T); return self._qrad.deriv(l, T)

    def k(self, T):
        l, _ = self._ln(T); return self._k(l)

    def dk(self, T):
        l, T = self._ln(T); return self._k.deriv(l, T)

    def S(self, T, T_ref):
        """Kirchhoff potential int_{T_ref}^{T} k dT'."""
        l, _ = self._ln(T)
        lr, _ = self._ln(np.asarray(T_ref, float))
        return self._S0(l) - self._S0(lr)

    def n_e(self, T):
        l, _ = self._ln(T); return self._ne(l)

    def ion_fraction(self, T):
        l, _ = self._ln(T); return self._ion(l)
