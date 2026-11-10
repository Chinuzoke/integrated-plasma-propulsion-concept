"""
validation.py -- analytical comparisons, numerical convergence, physical consistency.

WHAT THESE CHECKS ESTABLISH (and what they do not)
--------------------------------------------------
category 'math'        : the code correctly implements the stated equations
                         (agreement with independent analytic/numeric evaluations).
category 'consistency' : energy conservation, bounds, monotonic/limiting behaviour.
category 'literature'  : spot-checks of constants against tabulated values
                         (LIMITED; see notes in each check).
category 'inputs'      : impossible inputs are rejected.

Passing them does NOT validate the physical model against experiment. No
experimental data exist in this project: the Gaussian-profile assumption for
M^2 = 11.8, the pointing model, the catalogue optical values and the neglect of
jitter/turbulence/aberrations all remain UNVALIDATED.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import beam_propagation as bp
import config as cfg
import optical_components as oc
import pointing_error as pe
import power_budget as pb
import receiver as rx


@dataclass(frozen=True)
class Check:
    category: str
    name: str
    passed: bool
    detail: str


def _chk(cat, name, ok, detail):
    return Check(cat, name, bool(ok), detail)


def run_all(inp: cfg.SystemInputs | None = None) -> list:
    inp = inp or cfg.SystemInputs.from_config()
    out = []
    bud = pb.compute_budget(inp)
    w, a, d = bud.beam.beam_radius_m, inp.aperture_radius_m, bud.capture.displacement_m

    # ---------------- beam equations
    w0, lam, m2 = inp.waist_m, inp.wavelength_m, inp.m2
    zr = bp.rayleigh_range(w0, lam, m2)
    out.append(_chk("math", "w(0) = w0", math.isclose(
        bp.beam_radius(0.0, w0, lam, m2), w0, rel_tol=1e-15), ""))
    r = bp.beam_radius(zr, w0, lam, m2) / w0
    out.append(_chk("math", "w(z_R) = sqrt(2) w0", math.isclose(
        r, math.sqrt(2), rel_tol=1e-14), f"ratio={r:.15f}"))
    z_far = 1e4 * zr
    ratio = bp.beam_radius(z_far, w0, lam, m2) / \
        (bp.divergence_half_angle(w0, lam, m2) * z_far)
    out.append(_chk("math", "far field: w(z) -> theta z",
               abs(ratio - 1) < 1e-7, f"ratio={ratio:.10f}"))
    bpp = w0 * bp.divergence_half_angle(w0, lam, m2)
    out.append(_chk("math", "BPP = w0*theta = M^2 lambda/pi",
               math.isclose(bpp, m2 * lam / math.pi, rel_tol=1e-14), f"{bpp:.6e} m rad"))
    # M^2 scaling: M^2 = 1 gives ideal Gaussian z_R = pi w0^2 / lambda
    out.append(_chk("math", "M^2=1 reduces to ideal Gaussian",
                    math.isclose(bp.rayleigh_range(w0, lam, 1.0), math.pi * w0**2 / lam, rel_tol=1e-14), ""))
    out.append(_chk("math", "M^2 widens the beam (w(M2) > w(1))",
                    bp.beam_radius(inp.distance_m, w0, lam, m2) > bp.beam_radius(inp.distance_m, w0, lam, 1.0), ""))
    # intensity integrates to P (numerical 2-D check over a large disc)
    wt = bud.beam.beam_radius_m
    tot = rx.capture_fraction_radial(8.0 * wt, wt, 0.0, 8, 64)
    out.append(_chk("math", "intensity integrates to P (8w disc)",
               abs(tot - 1) < 1e-12, f"1-integral={1-tot:.2e}"))
    out.append(_chk("math", "I0 = 2P/(pi w^2)", math.isclose(
        bp.on_axis_intensity(1.0, wt) * math.pi * wt**2 / 2, 1.0, rel_tol=1e-14), ""))

    # ---------------- receiver: analytical vs numerical (centred)
    ana = rx.capture_fraction_centered_analytic(a, w)
    num = rx.capture_fraction_radial(
        a, w, 0.0, **{"n_panels": cfg.NUMERICS["n_panels"], "n_nodes": cfg.NUMERICS["n_nodes"]})
    out.append(_chk("math", "centred capture: numeric vs 1-exp(-2a^2/w^2)", abs(ana - num) < 1e-12,
                    f"analytic={ana:.12f} numeric={num:.12f} diff={abs(ana-num):.2e}"))
    # displaced: three independent methods
    for (aa, ww, dd) in [(a, w, d), (5.0, 3.0, 6.0), (5.0, 8.0, 4.0)]:
        e1 = rx.capture_fraction_radial(aa, ww, dd)
        e2 = rx.capture_fraction_polar(aa, ww, dd)
        e3 = rx.capture_fraction_reference(aa, ww, dd)
        out.append(_chk("math", f"displaced capture a={aa} w={ww} d={dd:.3g}: radial/polar/Marcum agree",
                        max(abs(e1 - e3), abs(e2 - e3)) < 1e-10,
                        f"radial={e1:.12f} polar={e2:.12f} ref={e3:.12f}"))
    # convergence
    ref, rows = rx.convergence_study(a, w, d, node_list=(2, 4, 8, 16, 32, 64))
    errs = [r_[2] for r_ in rows]
    out.append(_chk("math", "quadrature convergence (error -> <1e-12 as nodes grow)",
                    errs[-1] < 1e-12 and errs[-1] <= errs[0], "errors: " + ", ".join(f"{e:.1e}" for e in errs)))
    # grid sensitivity of the final answer
    base = rx.capture_fraction_radial(a, w, d, 4, 64)
    variants = [rx.capture_fraction_radial(a, w, d, p, n) for p, n in [
        (1, 32), (2, 48), (8, 96), (16, 128)]]
    out.append(_chk("consistency", "result insensitive to grid resolution (<1e-10)",
                    max(abs(v - base) for v in variants) < 1e-10, f"max change={max(abs(v-base) for v in variants):.2e}"))
    pv = rx.capture_fraction_polar(
        a, w, d, 4, 64, 128) - rx.capture_fraction_polar(a, w, d, 4, 64, 512)
    out.append(_chk("consistency", "2-D polar result insensitive to phi resolution",
               abs(pv) < 1e-10, f"diff={abs(pv):.2e}"))

    # ---------------- receiver: behaviour
    z0 = rx.capture_fraction_radial(
        a, w, pe.lateral_displacement(inp.distance_m, 0.0))
    out.append(_chk("consistency", "zero pointing error = centred result", abs(
        z0 - ana) < 1e-12, ""))
    ds = np.linspace(0, 3 * w, 40)
    e = [rx.capture_fraction_radial(a, w, float(x)) for x in ds]
    out.append(_chk("consistency", "capture decreases monotonically with pointing displacement",
                    all(e[i] >= e[i + 1] - 1e-15 for i in range(len(e) - 1)), ""))
    As = np.linspace(0.1, 6 * w, 40)
    e = [rx.capture_fraction_radial(float(x), w, d, 8, 64) for x in As]
    out.append(_chk("consistency", "capture increases monotonically with aperture and -> 1",
                    all(e[i] <= e[i + 1] + 1e-15 for i in range(len(e) - 1)) and e[-1] > 1 - 1e-9, f"final={e[-1]:.12f}"))
    out.append(_chk("consistency", "capture efficiency within [0,1]", 0.0 <=
               bud.capture.efficiency <= 1.0, f"{bud.capture.efficiency:.6f}"))
    out.append(_chk("math", "small aperture limit: eta -> 2a^2/w^2",
                    abs(rx.capture_fraction_radial(1e-3 * w, w, 0.0) / (2e-6) - 1) < 1e-5, ""))

    # ---------------- budget: conservation, chaining, bounds
    res = [abs(s.residual) / s.p_in for s in bud.stages]
    out.append(_chk("consistency", "energy conservation at every stage (rel. residual < 1e-12)",
                    max(res) < 1e-12, f"max rel residual={max(res):.1e}"))
    chained = all(bud.stages[i].p_out == bud.stages[i +
                  1].p_in for i in range(len(bud.stages) - 1))
    out.append(_chk("consistency",
               "stages chained in series (each input = previous output)", chained, ""))
    out.append(_chk("consistency", "global balance P_tx = P_delivered + sum of losses",
                    math.isclose(inp.power_w, bud.delivered_w + sum(s.total_loss for s in bud.stages), rel_tol=1e-12), ""))
    out.append(_chk("consistency", "overall efficiency in (0,1]",
               0.0 < bud.efficiency <= 1.0, f"{bud.efficiency:.6f}"))
    opt = [s for s in bud.stages if s.name in (
        "Receiving mirror", "Beam-conditioning assembly", "Fused-silica window")]
    prod = math.prod(s.p_out / s.p_in for s in opt)
    out.append(_chk("math", "optical train throughput = product of component throughputs",
                    math.isclose(prod, opt[-1].p_out / opt[0].p_in, rel_tol=1e-12), f"{prod:.6f}"))

    # ---------------- window physics
    n = oc.fused_silica_index(1064e-9)
    n_d = oc.fused_silica_index(587.56e-9)
    out.append(_chk("literature", "fused-silica n(587.56 nm) close to commonly tabulated 1.4585",
                    abs(n_d -
                        1.4585) < 5e-4, f"computed {n_d:.5f}; tabulated value, "
                    "coefficients from a secondary source: confirm against the paper"))
    T, Rt, A = oc.plate_fractions(n, None, None)
    out.append(_chk("math", "uncoated plate: T + R + A = 1 and T = (1-R)/(1+R) when tau=1",
                    math.isclose(T + Rt + A, 1.0, abs_tol=1e-14) and math.isclose(
                        T, (1 - oc.fresnel_reflectance_normal(n)) / (1 + oc.fresnel_reflectance_normal(n)), rel_tol=1e-13),
                    f"T={T:.5f} R={Rt:.5f}"))
    T2, R2, A2 = oc.plate_fractions(n, 0.01, 1.0)
    out.append(_chk("consistency", "bulk absorption lowers T and raises A",
               T2 < T and A2 > A, f"A={A2:.5f}"))

    # ---------------- register v1 consistency
    R_s = cfg.value("window_single_sided_reflectance")
    T_from_R = (1.0 - R_s) / (1.0 + R_s)
    out.append(_chk("literature", "window spec self-consistent: (1-R)/(1+R) with R=0.2 % >= stated T",
                    T_from_R >= cfg.value("window_transmission") - 1e-12,
                    f"derived T={T_from_R:.6f} vs stated {cfg.value('window_transmission')}; consistent with the "
                    "manufacturer deriving T from single-sided reflectance (manufacturer text as supplied in the register)"))
    typ = pb.compute_budget(inp.with_(mirror_R=cfg.value("mirror_typical_reflectance"),
                                      mirror_abs_scatter_ppm=cfg.value("mirror_typical_abs_scatter_ppm")))
    ms = typ.stages[2]
    out.append(_chk("consistency", "typical-mirror scenario conserves energy with a combined abs+scatter channel",
                    abs(ms.residual) < 1e-12 * ms.p_in and len(
                        ms.known_losses) == 1 and ms.combined_unresolved_loss >= 0,
                    f"combined A+S = {list(ms.known_losses.values())[0]:.4f} W; unattributed remainder = {ms.combined_unresolved_loss:.4f} W"))
    out.append(_chk("consistency", "catalogue window applies stated T once (no extra single-sided reflectance)",
                    math.isclose(bud.stages[-1].p_out / bud.stages[-1].p_in, inp.window_T, rel_tol=1e-14), ""))
    out.append(_chk("consistency", "conditioning assembly loss not split among surfaces",
                    bud.stages[3].known_losses == {}, ""))

    # ---------------- input rejection
    bad = [dict(power_w=-1.0), dict(aperture_radius_m=-5.0), dict(cond_T=1.01), dict(mirror_R=0.0),
           dict(m2=0.5), dict(window_T=1.5), dict(waist_m=0.0), dict(
               distance_m=-1.0), dict(pointing_arcsec=-0.1),
           dict(wavelength_m=float("nan"))]
    ok = True
    for b in bad:
        try:
            inp.with_(**b)
            ok = False
        except ValueError:
            pass
    out.append(_chk(
        "inputs", f"{len(bad)} impossible inputs rejected (neg power/aperture, T>1, M^2<1, ...)", ok, ""))
    try:
        cfg.value("window_thickness")
        ok = False
    except ValueError:
        ok = True
    out.append(
        _chk("inputs", "UNRESOLVED parameters refuse to return a value", ok, ""))
    return out


def summary(checks: list) -> str:
    n_ok = sum(c.passed for c in checks)
    return f"{n_ok}/{len(checks)} checks passed"
