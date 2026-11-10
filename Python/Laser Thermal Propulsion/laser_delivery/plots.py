"""plots.py -- figures (matplotlib, Agg backend). All axes in SI units."""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import beam_propagation as bp
import config as cfg
import receiver as rx
import uncertainty as un


def _save(fig, outdir, name):
    path = os.path.join(outdir, name)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_intensity(bud, outdir):
    inp, w, a = bud.inputs, bud.beam.beam_radius_m, bud.inputs.aperture_radius_m
    d = bud.capture.displacement_m
    P = bud.stages[0].p_out
    lim = 2.0 * w
    x = np.linspace(-lim, lim, 401)
    X, Y = np.meshgrid(x, x)
    I = bp.on_axis_intensity(P, w) * np.exp(-2 * ((X - d) ** 2 + Y**2) / w**2)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    im = ax[0].pcolormesh(x, x, I, shading="auto", cmap="inferno")
    ax[0].add_patch(plt.Circle((0, 0), a, fill=False, color="cyan", lw=1.5, label=f"receiver aperture (a = {a:g} m)"))
    ax[0].plot([d], [0], "w+", ms=10, label=f"beam centre (d = {d:.3f} m)")
    ax[0].set_aspect("equal"); ax[0].set_xlabel("x at receiver plane [m]"); ax[0].set_ylabel("y [m]")
    ax[0].set_title("Beam intensity at receiver plane"); ax[0].legend(loc="upper right", fontsize=7)
    fig.colorbar(im, ax=ax[0], label="intensity [W/m$^2$]")
    r = np.linspace(0, lim, 400)
    ax[1].plot(r, bp.intensity(r, P, w), label="Gaussian, M$^2$-corrected radius")
    ax[1].axvline(a, color="c", ls="--", label="aperture radius")
    ax[1].axvline(w, color="gray", ls=":", label="w (1/e$^2$ radius)")
    ax[1].set_xlabel("radius r [m]"); ax[1].set_ylabel("intensity [W/m$^2$]")
    ax[1].set_title("Radial profile (centred)"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    return _save(fig, outdir, "01_intensity_at_receiver.png")


def plot_capture_vs_pointing(bud, outdir):
    inp, w, a = bud.inputs, bud.beam.beam_radius_m, bud.inputs.aperture_radius_m
    ds = np.linspace(0, 3.0 * w, 300)
    eta = [100 * rx.capture_fraction_radial(a, w, float(x)) for x in ds]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.plot(ds, eta)
    ax.plot([bud.capture.displacement_m], [100 * bud.capture.efficiency], "ro", label="operating point")
    ax.set_xlabel("lateral pointing displacement at receiver, d [m]"); ax.set_ylabel("capture efficiency [%]")
    ax.set_title("Receiver capture vs pointing displacement"); ax.grid(alpha=.3); ax.legend()
    k = inp.distance_m * cfg.ARCSEC_TO_RAD
    sec = ax.secondary_xaxis("top", functions=(lambda x: x / k, lambda y: y * k))
    sec.set_xlabel("equivalent pointing error [arcsec]")
    return _save(fig, outdir, "02_capture_vs_pointing.png")


def plot_capture_vs_aperture(bud, outdir):
    w, d = bud.beam.beam_radius_m, bud.capture.displacement_m
    a_s = np.linspace(0.1, 3.0 * w, 300)
    cen = [100 * rx.capture_fraction_centered_analytic(float(x), w) for x in a_s]
    dis = [100 * rx.capture_fraction_radial(float(x), w, d, 8, 64) for x in a_s]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.plot(a_s, cen, label="centred (analytic)")
    ax.plot(a_s, dis, "--", label=f"displaced d = {d:.3f} m (numerical)")
    ax.axvline(bud.inputs.aperture_radius_m, color="gray", ls=":", label="baseline aperture")
    ax.set_xlabel("receiver aperture radius a [m]"); ax.set_ylabel("capture efficiency [%]")
    ax.set_title("Capture efficiency vs receiver aperture radius"); ax.grid(alpha=.3); ax.legend()
    return _save(fig, outdir, "03_capture_vs_aperture.png")


def plot_loss_chart(bud, outdir):
    ptx = bud.inputs.power_w
    names = ["Transmitter"] + [s.name for s in bud.stages]
    remaining = [ptx] + [s.p_out for s in bud.stages]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    ax[0].bar(range(len(names)), [x / 1e3 for x in remaining], color="tab:blue")
    ax[0].set_xticks(range(len(names))); ax[0].set_xticklabels(names, rotation=30, ha="right", fontsize=8)
    ax[0].set_ylabel("power after stage [kW]"); ax[0].set_title("Power remaining along the chain"); ax[0].grid(alpha=.3, axis="y")
    opt = [s for s in bud.stages if s.name not in ("Free-space propagation", "Receiver aperture capture")]
    ax[1].bar([s.name for s in opt], [s.total_loss for s in opt], color="tab:red")
    ax[1].set_ylabel("power lost in stage [W]"); ax[1].set_title("Optical-component losses (given captured power)")
    ax[1].tick_params(axis="x", rotation=15, labelsize=8); ax[1].grid(alpha=.3, axis="y")
    return _save(fig, outdir, "04_power_loss_chart.png")


def plot_sensitivity(inp, outdir):
    el = un.local_elasticities(inp)
    tor = un.tornado(inp)
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    names = sorted(el, key=lambda k: abs(el[k]))
    ax[0].barh([un.PARAM_LABELS[k] for k in names], [el[k] for k in names], color="tab:green")
    ax[0].set_xlabel("elasticity  d ln P$_{del}$ / d ln x  [-]"); ax[0].set_title("Local sensitivity at baseline"); ax[0].grid(alpha=.3, axis="x")
    base = un.delivered_power(inp)
    rows = tor[::-1]
    for i, (n, lo, hi, plo, phi, sw) in enumerate(rows):
        ax[1].barh(i, plo - base, left=base, color="tab:orange", alpha=.8)
        ax[1].barh(i, phi - base, left=base, color="tab:blue", alpha=.8)
    ax[1].set_yticks(range(len(rows))); ax[1].set_yticklabels([un.PARAM_LABELS[r[0]] for r in rows])
    ax[1].axvline(base, color="k", lw=.8)
    ax[1].set_xlabel("delivered power [W]"); ax[1].set_title("One-at-a-time sweep, ILLUSTRATIVE ranges (not measured uncertainties)", fontsize=9)
    ax[1].grid(alpha=.3, axis="x")
    return _save(fig, outdir, "05_sensitivity.png")


def make_all(bud, outdir):
    os.makedirs(outdir, exist_ok=True)
    return [plot_intensity(bud, outdir), plot_capture_vs_pointing(bud, outdir),
            plot_capture_vs_aperture(bud, outdir), plot_loss_chart(bud, outdir),
            plot_sensitivity(bud.inputs, outdir)]
