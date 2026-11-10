"""
optical_components.py -- component definitions and individual transmission models.

Design rule (no double counting, no invented splits)
----------------------------------------------------
A component is modelled at the finest level for which INDEPENDENT data exist.

* LumpedComponent: only an assembly/component-level total (reflectance or
  transmission) is available (receiving mirror, beam-conditioning assembly,
  catalogue window). The loss 1 - throughput is reported as ONE combined
  channel whose split into reflection/absorption/scatter/transmission is
  explicitly "unavailable". Nothing is divided among surfaces.
* FresnelWindow: physical model for an UNCOATED fused-silica plate: two
  Fresnel interfaces with incoherent multiple reflection plus optional bulk
  absorption exp(-alpha t). Scattering is not modelled (no data).

Catalogue values are specifications for the project model, NOT measured
performance of the assembled system.

Fresnel / window equations (normal incidence, incoherent, no scatter):
    R   = ((n - 1)/(n + 1))^2                  per interface (air-glass)
    tau = exp(-alpha t)                         single-pass internal transmittance
    T   = (1-R)^2 tau / (1 - R^2 tau^2)
    Rtot= R + (1-R)^2 R tau^2 / (1 - R^2 tau^2)
    A   = 1 - T - Rtot                          (>= 0)
Fused-silica index from the Sellmeier dispersion formula of Malitson (1965),
20 degC, valid 0.21-3.71 um. Coefficients as implemented in the PyXOpto
library (secondary source); verify against the original paper (see references).
Neglected: coherent (etalon) effects, polarisation (irrelevant at normal
incidence), angle of incidence, temperature dependence, thermal lensing.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import config as cfg


@dataclass(frozen=True)
class StageResult:
    """Power accounting for one stage. Conservation: p_in = p_out + known + unresolved."""
    name: str
    p_in: float
    p_out: float
    # channel -> W (independently modelled)
    known_losses: dict = field(default_factory=dict)
    combined_unresolved_loss: float = 0.0              # W, split unavailable
    unresolved_note: str = ""
    basis: str = ""                                    # provenance of the numbers
    notes: tuple = ()

    @property
    def total_loss(self) -> float:
        return sum(self.known_losses.values()) + self.combined_unresolved_loss

    @property
    def residual(self) -> float:
        """p_in - (p_out + all losses); must be ~0."""
        return self.p_in - (self.p_out + self.total_loss)


@dataclass(frozen=True)
class LumpedComponent:
    """
    Component described only by a total throughput. Optionally, independently
    supported loss channels (fractions of incident power) may be given in
    `known_loss_fractions`; the rest of 1 - throughput stays one COMBINED,
    unattributed channel. Nothing is split without a source.
    """
    name: str
    throughput: float          # fraction of incident power delivered onward
    basis: str                 # provenance of the number
    unresolved_note: str
    notes: tuple = ()
    # channel -> fraction of incident power
    known_loss_fractions: dict = field(default_factory=dict)

    def __post_init__(self):
        cfg.require_fraction(f"{self.name} throughput", self.throughput)
        for k, f in self.known_loss_fractions.items():
            cfg.require_nonnegative(f"{self.name} known channel '{k}'", f)
        if sum(self.known_loss_fractions.values()) > (1.0 - self.throughput) + 1e-15:
            raise ValueError(
                f"{self.name}: known loss channels exceed the total loss 1 - throughput")

    def propagate(self, p_in: float) -> StageResult:
        cfg.require_positive(f"{self.name} incident power", p_in)
        p_out = p_in * self.throughput
        known = {k: p_in * f for k, f in self.known_loss_fractions.items()}
        combined = max(p_in - p_out - sum(known.values()), 0.0)
        return StageResult(self.name, p_in, p_out, known, combined,
                           self.unresolved_note, self.basis, self.notes)


# ----------------------------------------------------------------- fused silica
_SELLMEIER = ((0.6961663, 0.0684043), (0.4079426,
              0.1162414), (0.8974794, 9.896161))


def fused_silica_index(wavelength_m: float) -> float:
    """Refractive index of fused silica, Malitson (1965) Sellmeier form, 20 degC."""
    cfg.require_positive("wavelength", wavelength_m)
    lam_um = wavelength_m * 1e6
    if not (0.21 <= lam_um <= 3.71):
        raise ValueError(
            f"wavelength {lam_um} um outside Sellmeier validity (0.21-3.71 um)")
    l2 = lam_um**2
    n2 = 1.0 + sum(b * l2 / (l2 - c**2) for b, c in _SELLMEIER)
    return math.sqrt(n2)


def fresnel_reflectance_normal(n: float) -> float:
    """Single air-glass interface reflectance at normal incidence."""
    if n < 1.0 or not math.isfinite(n):
        raise ValueError("refractive index must be >= 1")
    return ((n - 1.0) / (n + 1.0)) ** 2


def plate_fractions(n: float, thickness_m: float | None, alpha_per_m: float | None):
    """(T, R_total, A) for an uncoated plate. If alpha or thickness is None: tau = 1."""
    R = fresnel_reflectance_normal(n)
    tau = 1.0
    if thickness_m is not None and alpha_per_m is not None:
        tau = math.exp(-alpha_per_m * thickness_m)
    denom = 1.0 - (R * tau) ** 2
    T = (1.0 - R) ** 2 * tau / denom
    Rt = R + (1.0 - R) ** 2 * R * tau**2 / denom
    A = 1.0 - T - Rt
    if A < -1e-12:
        raise ArithmeticError("negative absorption: inconsistent plate model")
    return T, Rt, max(A, 0.0)


@dataclass(frozen=True)
class FresnelWindow:
    """Uncoated fused-silica window, physical model (no catalogue number used)."""
    name: str
    wavelength_m: float
    thickness_m: float | None = None
    alpha_per_m: float | None = None

    def propagate(self, p_in: float) -> StageResult:
        cfg.require_positive(f"{self.name} incident power", p_in)
        n = fused_silica_index(self.wavelength_m)
        T, Rt, A = plate_fractions(n, self.thickness_m, self.alpha_per_m)
        notes = [f"n = {n:.5f} (Sellmeier), R per surface = {fresnel_reflectance_normal(n)*100:.3f} %",
                 "UNCOATED fused silica; scattering not modelled (no data)"]
        if self.thickness_m is None or self.alpha_per_m is None:
            notes.append("bulk absorption NOT included (thickness and/or alpha unresolved): "
                         "this is a lower bound on window loss")
        return StageResult(self.name, p_in, p_in * T,
                           {"reflected (both surfaces, multiple pass)": p_in * Rt,
                            "bulk absorbed": p_in * A},
                           0.0, "", "physical model (Fresnel + optional bulk absorption)",
                           tuple(notes))
