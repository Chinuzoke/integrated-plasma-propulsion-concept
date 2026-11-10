"""
receiver.py -- power captured by a circular aperture from a (displaced) Gaussian beam.

Physics
-------
Plane intensity, beam centred at (d, 0) relative to the aperture centre:

    I(x,y) = (2P / pi w^2) exp( -2[(x-d)^2 + y^2] / w^2 )

Captured power  P_cap = integral of I over the disc x^2+y^2 <= a^2.

Three independent evaluations are implemented:

1. PRODUCTION (capture_fraction_radial): polar coordinates (r,phi) centred on
   the aperture. The angular integral has the closed form
       int_0^{2pi} exp(4 r d cos(phi)/w^2) dphi = 2 pi I0(4 r d / w^2),
   leaving a 1-D radial integral, evaluated by composite Gauss-Legendre
   quadrature:
       eta = (4/w^2) int_0^a r exp(-2 (r-d)^2/w^2) i0e(4 r d/w^2) dr
   (i0e(x) = exp(-x) I0(x), exponentially scaled Bessel: no overflow).
2. CROSS-CHECK (capture_fraction_polar): full 2-D tensor-product quadrature
   (Gauss-Legendre in r, periodic trapezoid in phi) with NO Bessel function.
3. ANALYTIC REFERENCE (capture_fraction_reference): per-axis sigma = w/2, so
   |X|^2/sigma^2 is non-central chi-square (2 dof, non-centrality (d/sigma)^2);
   eta = ncx2.cdf((a/sigma)^2, 2, (d/sigma)^2) (a Marcum-Q-function form).
   For d = 0 this reduces exactly to eta = 1 - exp(-2 a^2 / w^2).

The centred closed form is never used for displaced beams.

Assumptions: hard-edged, unobstructed aperture; perfect Gaussian profile; the
aperture is perpendicular to the beam (no cos(theta) projection: theta ~ 1e-5).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import i0e
from scipy.stats import ncx2

import config as cfg


# ------------------------------------------------------------------ quadrature
def _gl_panels(a: float, n_panels: int, n_nodes: int):
    """Composite Gauss-Legendre nodes/weights on [0, a]."""
    if n_panels < 1 or n_nodes < 2:
        raise ValueError("need n_panels >= 1 and n_nodes >= 2")
    x, wt = np.polynomial.legendre.leggauss(n_nodes)
    edges = np.linspace(0.0, a, n_panels + 1)
    half = 0.5 * (edges[1:] - edges[:-1])
    mid = 0.5 * (edges[1:] + edges[:-1])
    r = (mid[:, None] + half[:, None] * x[None, :]).ravel()
    wr = (half[:, None] * wt[None, :]).ravel()
    return r, wr


def _check(a: float, w: float, d: float) -> None:
    cfg.require_positive("aperture radius", a)
    cfg.require_positive("beam radius", w)
    cfg.require_nonnegative("displacement", d)


def capture_fraction_radial(a: float, w: float, d: float = 0.0,
                            n_panels: int = 4, n_nodes: int = 64) -> float:
    """Production method: 1-D radial Gauss-Legendre with analytic angular integral."""
    _check(a, w, d)
    r, wr = _gl_panels(a, n_panels, n_nodes)
    integrand = (4.0 / w**2) * r * np.exp(-2.0 * (r - d) ** 2 / w**2) * i0e(4.0 * r * d / w**2)
    return float(np.sum(wr * integrand))


def capture_fraction_polar(a: float, w: float, d: float = 0.0,
                           n_panels: int = 4, n_nodes: int = 64, n_phi: int = 256) -> float:
    """Independent full 2-D polar quadrature (no Bessel function)."""
    _check(a, w, d)
    r, wr = _gl_panels(a, n_panels, n_nodes)
    phi = np.arange(n_phi) * (2.0 * math.pi / n_phi)
    R, PHI = np.meshgrid(r, phi, indexing="ij")
    inten = (2.0 / (math.pi * w**2)) * np.exp(-2.0 * (R**2 + d**2 - 2.0 * R * d * np.cos(PHI)) / w**2)
    return float(np.sum(wr[:, None] * R * inten) * (2.0 * math.pi / n_phi))


def capture_fraction_centered_analytic(a: float, w: float) -> float:
    """eta = 1 - exp(-2 a^2 / w^2). Valid ONLY for a centred beam."""
    cfg.require_positive("aperture radius", a)
    cfg.require_positive("beam radius", w)
    return -math.expm1(-2.0 * a**2 / w**2)


def capture_fraction_reference(a: float, w: float, d: float = 0.0) -> float:
    """Analytic reference (non-central chi-square / Marcum-Q form), valid for any d."""
    _check(a, w, d)
    if d == 0.0:
        return capture_fraction_centered_analytic(a, w)
    sigma = 0.5 * w
    return float(ncx2.cdf((a / sigma) ** 2, 2, (d / sigma) ** 2))


# ------------------------------------------------------------------ result
@dataclass(frozen=True)
class CaptureResult:
    p_incident_w: float
    p_captured_w: float
    p_missed_w: float
    efficiency: float
    displacement_m: float
    beam_radius_m: float
    aperture_radius_m: float
    quadrature_self_error: float   # |eta(N) - eta(N/2 nodes)|: numerical convergence evidence
    method_discrepancy: float      # |radial - 2-D polar|: independent-method agreement


def capture_power(p_in_w: float, w: float, a: float, d: float = 0.0,
                  n_panels: int | None = None, n_nodes: int | None = None,
                  n_phi: int | None = None) -> CaptureResult:
    cfg.require_positive("incident power", p_in_w)
    n_panels = n_panels or cfg.NUMERICS["n_panels"]
    n_nodes = n_nodes or cfg.NUMERICS["n_nodes"]
    n_phi = n_phi or cfg.NUMERICS["n_phi"]
    eta = capture_fraction_radial(a, w, d, n_panels, n_nodes)
    eta_coarse = capture_fraction_radial(a, w, d, n_panels, max(2, n_nodes // 2))
    eta_polar = capture_fraction_polar(a, w, d, n_panels, n_nodes, n_phi)
    if not (-1e-12 <= eta <= 1.0 + 1e-9):
        raise ArithmeticError(f"capture fraction out of physical bounds: {eta}")
    eta = min(max(eta, 0.0), 1.0)  # remove round-off only (checked above)
    return CaptureResult(
        p_incident_w=p_in_w,
        p_captured_w=p_in_w * eta,
        p_missed_w=p_in_w * (1.0 - eta),
        efficiency=eta,
        displacement_m=d,
        beam_radius_m=w,
        aperture_radius_m=a,
        quadrature_self_error=abs(eta - eta_coarse),
        method_discrepancy=abs(eta - eta_polar),
    )


# ------------------------------------------------------------------ studies
def convergence_study(a: float, w: float, d: float, node_list=(4, 8, 16, 32, 64, 128),
                      n_panels: int = 1):
    """Return rows (n_nodes, eta_radial, |eta - reference|, eta_polar, |polar - reference|)."""
    ref = capture_fraction_reference(a, w, d)
    rows = []
    for n in node_list:
        e1 = capture_fraction_radial(a, w, d, n_panels, n)
        e2 = capture_fraction_polar(a, w, d, n_panels, n, n_phi=max(16, 4 * n))
        rows.append((n, e1, abs(e1 - ref), e2, abs(e2 - ref)))
    return ref, rows


def mean_capture_with_jitter(a: float, w: float, sigma_axis_m: float,
                             n_nodes: int = 96, d_max_sigmas: float = 9.0) -> float:
    """
    Expected capture if the pointing error is zero-mean Gaussian jitter with
    RMS displacement sigma_axis_m PER AXIS (so the radial displacement d has a
    Rayleigh pdf  p(d) = d/s^2 exp(-d^2/2s^2)). Informational alternative to the
    fixed-offset assumption.
    """
    cfg.require_positive("aperture radius", a)
    cfg.require_positive("beam radius", w)
    cfg.require_nonnegative("sigma", sigma_axis_m)
    if sigma_axis_m == 0.0:
        return capture_fraction_radial(a, w, 0.0)
    s = sigma_axis_m
    x, wt = np.polynomial.legendre.leggauss(n_nodes)
    dmax = d_max_sigmas * s
    d = 0.5 * dmax * (x + 1.0)
    wd = 0.5 * dmax * wt
    pdf = d / s**2 * np.exp(-0.5 * d**2 / s**2)
    eta = np.array([capture_fraction_radial(a, w, float(di)) for di in d])
    norm = float(np.sum(wd * pdf))  # ~1 - exp(-d_max^2/2s^2)
    return float(np.sum(wd * pdf * eta) / norm)
