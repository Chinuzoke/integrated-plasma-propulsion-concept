"""
beam_propagation.py -- Gaussian beam propagation with M^2 correction (SI units).

Model: the "embedded Gaussian" description of a real beam (Siegman). A beam
with beam-propagation ratio M^2 propagates like an ideal Gaussian whose
wavelength is replaced by M^2 * lambda. The beam radius is the second-moment
(1/e^2-equivalent) radius.

Equations (all in SI)
---------------------
Rayleigh range     z_R   = pi w0^2 / (M^2 lambda)                      [m]
    Distance over which the beam area doubles. Marks the near/far-field
    transition.
Beam radius        w(z)  = w0 sqrt(1 + (z/z_R)^2)                      [m]
    Exact hyperbola for the embedded Gaussian. Valid in vacuum, paraxial.
Far-field half-angle  theta = M^2 lambda / (pi w0)                     [rad]
    Asymptote of w(z)/z for z >> z_R (small-angle). BPP = w0*theta.
Intensity          I(r,z) = (2P / pi w(z)^2) exp(-2 r^2 / w(z)^2)      [W/m^2]
    Integrates to P over the plane. On-axis I0 = 2P / (pi w^2).

Validity limits / caveats
-------------------------
* Paraxial: theta << 1 rad (here ~1.6e-5 rad: excellent).
* M^2 = 11.8 describes a strongly multimode beam. The Gaussian *transverse
  profile* is then an assumption, not a measurement: a real multimode beam
  can be flat-topped or structured, which changes captured power for an
  aperture smaller than the beam. Second-moment radius propagation (w(z)) is
  robust; the profile shape is the model-form uncertainty.
* Scalar M^2 assumes identical x and y behaviour (stigmatic beam).
* Vacuum, no turbulence, no thermal blooming, no optics-induced aberration.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import config as cfg


def beam_parameter_product(lam: float, m2: float) -> float:
    """BPP = M^2 lambda / pi = w0 * theta   [m rad]. Conserved in free propagation."""
    cfg.require_positive("wavelength", lam)
    cfg.require_positive("M^2", m2)
    return m2 * lam / math.pi


def rayleigh_range(w0: float, lam: float, m2: float) -> float:
    """z_R = pi w0^2 / (M^2 lambda)  [m]."""
    cfg.require_positive("waist", w0)
    cfg.require_positive("wavelength", lam)
    cfg.require_positive("M^2", m2)
    return math.pi * w0**2 / (m2 * lam)


def divergence_half_angle(w0: float, lam: float, m2: float) -> float:
    """Far-field half-angle theta = M^2 lambda / (pi w0)  [rad]."""
    cfg.require_positive("waist", w0)
    return beam_parameter_product(lam, m2) / w0


def beam_radius(z, w0: float, lam: float, m2: float):
    """w(z) = w0 sqrt(1 + (z/z_R)^2)  [m]. Accepts scalars or numpy arrays (z >= 0)."""
    zr = rayleigh_range(w0, lam, m2)
    z = np.asarray(z, dtype=float)
    if np.any(z < 0) or not np.all(np.isfinite(z)):
        raise ValueError("propagation distance must be finite and >= 0")
    w = w0 * np.sqrt(1.0 + (z / zr) ** 2)
    return float(w) if w.ndim == 0 else w


def beam_diameter(z, w0: float, lam: float, m2: float):
    """2 w(z)  [m]."""
    return 2.0 * beam_radius(z, w0, lam, m2)


def on_axis_intensity(power_w: float, w: float) -> float:
    """I0 = 2P / (pi w^2)  [W/m^2]."""
    cfg.require_positive("power", power_w)
    cfg.require_positive("beam radius", w)
    return 2.0 * power_w / (math.pi * w**2)


def intensity(r, power_w: float, w: float):
    """Gaussian intensity at radius r (scalar/array) in a plane where the radius is w."""
    r = np.asarray(r, dtype=float)
    return on_axis_intensity(power_w, w) * np.exp(-2.0 * r**2 / w**2)


@dataclass(frozen=True)
class BeamState:
    """Beam parameters at the receiver plane."""
    distance_m: float
    rayleigh_range_m: float
    divergence_half_angle_rad: float
    beam_radius_m: float
    beam_diameter_m: float
    on_axis_intensity_w_m2: float
    z_over_zR: float
    far_field: bool  # True if z >= 10 z_R (asymptotic regime)


def beam_state(power_w: float, lam: float, w0: float, m2: float, z: float) -> BeamState:
    zr = rayleigh_range(w0, lam, m2)
    w = beam_radius(z, w0, lam, m2)
    return BeamState(
        distance_m=z,
        rayleigh_range_m=zr,
        divergence_half_angle_rad=divergence_half_angle(w0, lam, m2),
        beam_radius_m=w,
        beam_diameter_m=2.0 * w,
        on_axis_intensity_w_m2=on_axis_intensity(power_w, w),
        z_over_zR=z / zr,
        far_field=(z / zr) >= 10.0,
    )
