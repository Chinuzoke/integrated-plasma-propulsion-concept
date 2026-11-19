"""
physics.py -- constitutive closures bundled for the energy equation.

A Physics object exposes, as functions of the cell temperatures T (arrays), the quantities and their
T-derivatives the finite-volume residual/Jacobian need:

    h, dh        specific enthalpy [J/kg] and dh/dT
    k, dk        thermal conductivity [W/m/K] and dk/dT
    alpha, dalpha  IB absorption [1/m] and derivative
    qrad, dqrad    volumetric radiative loss [W/m^3]
    qwall, dqwall  radial/wall loss [W/m^3]
    q_ext          uniform external source [W/m^3] (0 except in analytic verification tests)
    rho, cp        for the pseudo-time scaling and diagnostics

Three implementations:
    ReconstructedPhysics   enthalpy incl. ionisation; LJ heavy-particle k; IB alpha; brems; wall closure
    LegacyEquivalentPhysics  the earlier model's constitutive set (frozen cp, power-law mu Eucken k, tuned wall
                             terms) so the NEW residual-based solver can be applied to the OLD physics
    ConstantTestPhysics    constant cp and k, no radiation/absorption (analytic verification)

Wall closures (reconstructed physics):
    kirchhoff   q_w = 4 S(T)/R^2,  S(T) = int_{T_w}^{T} k dT'.  Derived from steady radial conduction in a
                column of radius R with uniform volumetric source and T(R) = T_w (Elenbaas-Heller type
                conduction potential). Parameter-free GIVEN k(T), but k is provisional (see transport.py).
    adiabatic   q_w = 0 (upper bound on gas heating; no wall loss).
    legacy_newton  legacy tuned lumped terms (h, emissivity): NOT physically justified; legacy only.
"""
from __future__ import annotations

import numpy as np

from . import constants as C
from . import argon_thermo as at
from . import transport as tr
from .config import LSPConfig
from .tables import PropertyTable


class Physics:
    name = "abstract"
    q_ext = 0.0
    nec_used = False

    def cp(self, T):
        return self.dh(T)


class ReconstructedPhysics(Physics):
    name = "reconstructed"

    def __init__(self, cfg: LSPConfig, k_override=None, nec_table=None):
        self.cfg = cfg
        self.tab = PropertyTable(cfg.pressure_pa, cfg.laser.wavelength_m, cfg.gaunt_ff, cfg.gaunt_brems,
                                 cfg.k_multiplier, k_override, nec_table)
        self.nec_used = nec_table is not None
        self.R = cfg.plasma_radius_m
        self.Tw = cfg.T_wall_K
        self.closure = cfg.wall_closure

    h = property(lambda self: self.tab.h)
    dh = property(lambda self: self.tab.dh)
    rho = property(lambda self: self.tab.rho)
    k = property(lambda self: self.tab.k)
    dk = property(lambda self: self.tab.dk)
    alpha = property(lambda self: self.tab.alpha)
    dalpha = property(lambda self: self.tab.dalpha)
    qrad = property(lambda self: self.tab.qrad)
    dqrad = property(lambda self: self.tab.dqrad)

    def qwall(self, T):
        T = np.asarray(T, float)
        if self.closure == "adiabatic":
            return np.zeros_like(T)
        if self.closure == "kirchhoff":
            return 4.0 * self.tab.S(T, self.Tw) / self.R**2
        h, eps = self.cfg.legacy_wall_h, self.cfg.legacy_wall_emissivity
        return (2 * h / self.R) * (T - self.Tw) + (2 * eps * C.sigma_SB / self.R) * (T**4 - self.Tw**4)

    def dqwall(self, T):
        T = np.asarray(T, float)
        if self.closure == "adiabatic":
            return np.zeros_like(T)
        if self.closure == "kirchhoff":
            return 4.0 * self.tab.k(T) / self.R**2
        h, eps = self.cfg.legacy_wall_h, self.cfg.legacy_wall_emissivity
        return 2 * h / self.R + (8 * eps * C.sigma_SB / self.R) * T**3


class LegacyEquivalentPhysics(Physics):
    """Earlier model's constitutive set; absorption and brems come from the same (verified) SI formulas."""
    name = "legacy_equivalent"

    def __init__(self, cfg: LSPConfig):
        self.cfg = cfg
        base = PropertyTable(cfg.pressure_pa, cfg.laser.wavelength_m, cfg.gaunt_ff, 1.2, 1.0)
        self.tab = base
        self.cp_f = at.frozen_cp()
        self.R, self.Tw = cfg.plasma_radius_m, cfg.T_wall_K

    def h(self, T):
        return self.cp_f * np.asarray(T, float)

    def dh(self, T):
        return np.full_like(np.asarray(T, float), self.cp_f)

    def rho(self, T):
        return self.cfg.pressure_pa * C.AR_MOLAR_MASS / (C.R_gas * np.asarray(T, float))

    def k(self, T):
        return tr.legacy_k(T)

    def dk(self, T):
        T = np.asarray(T, float)
        return 0.72 * tr.legacy_k(T) / T      # power-law in T

    alpha = property(lambda self: self.tab.alpha)
    dalpha = property(lambda self: self.tab.dalpha)
    qrad = property(lambda self: self.tab.qrad)
    dqrad = property(lambda self: self.tab.dqrad)

    def qwall(self, T):
        T = np.asarray(T, float)
        h, eps = self.cfg.legacy_wall_h, self.cfg.legacy_wall_emissivity
        return (2 * h / self.R) * (T - self.Tw) + (2 * eps * C.sigma_SB / self.R) * (T**4 - self.Tw**4)

    def dqwall(self, T):
        T = np.asarray(T, float)
        h, eps = self.cfg.legacy_wall_h, self.cfg.legacy_wall_emissivity
        return 2 * h / self.R + (8 * eps * C.sigma_SB / self.R) * T**3


class ConstantTestPhysics(Physics):
    """Constant cp and k, no absorption/radiation/wall; optional uniform source q_ext. For analytic tests."""
    name = "constant_test"

    def __init__(self, cp=520.0, k=0.1, q_ext=0.0, rho=0.02):
        self.cp_c, self.k_c, self.q_ext, self.rho_c = cp, k, q_ext, rho

    def h(self, T): return self.cp_c * np.asarray(T, float)
    def dh(self, T): return np.full_like(np.asarray(T, float), self.cp_c)
    def rho(self, T): return np.full_like(np.asarray(T, float), self.rho_c)
    def k(self, T): return np.full_like(np.asarray(T, float), self.k_c)
    def dk(self, T): return np.zeros_like(np.asarray(T, float))
    def alpha(self, T): return np.zeros_like(np.asarray(T, float))
    def dalpha(self, T): return np.zeros_like(np.asarray(T, float))
    def qrad(self, T): return np.zeros_like(np.asarray(T, float))
    def dqrad(self, T): return np.zeros_like(np.asarray(T, float))
    def qwall(self, T): return np.zeros_like(np.asarray(T, float))
    def dqwall(self, T): return np.zeros_like(np.asarray(T, float))


def build_physics(cfg: LSPConfig, **kw) -> Physics:
    return LegacyEquivalentPhysics(cfg) if cfg.physics == "legacy_equivalent" else ReconstructedPhysics(cfg, **kw)
