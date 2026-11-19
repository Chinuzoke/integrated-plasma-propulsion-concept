"""plots.py -- LSP figures (SI units)."""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


from .outputs import save_figure


def _save(fig, outdir, name):
    """Save + verify (exists, non-empty); returns the absolute path or raises OutputError."""
    return save_figure(fig, os.path.join(outdir, name), dpi=140)


def plot_profiles(res, outdir, prefix="lsp"):
    p = res.profiles
    x = 1e3 * p["x"]
    fig, ax = plt.subplots(2, 2, figsize=(11.5, 7.5))
    ax[0, 0].plot(x, p["T"]); ax[0, 0].set_ylabel("temperature [K]"); ax[0, 0].set_title("Axial temperature")
    ax[0, 1].semilogy(x, np.maximum(p["n_e"], 1e-10)); ax[0, 1].set_ylabel("electron density [m$^{-3}$]")
    ax2 = ax[0, 1].twinx(); ax2.plot(x, p["ion_fraction"], "r--"); ax2.set_ylabel("ionisation fraction n$_e$/n$_{heavy}$", color="r")
    ax[0, 1].set_title("Electron density and ionisation")
    ax[1, 0].plot(x, p["I"][:-1] / p["I"][0], label="I(x)/I$_0$ (cell entry)"); ax[1, 0].set_ylabel("laser intensity / entrance intensity")
    ax3 = ax[1, 0].twinx(); ax3.semilogy(x, np.maximum(p["alpha"], 1e-12), "g:"); ax3.set_ylabel("absorption coefficient [1/m]", color="g")
    ax[1, 0].set_title("Laser attenuation (Beer-Lambert)")
    ax[1, 1].plot(x, p["q_laser"], label="laser deposition"); ax[1, 1].plot(x, p["q_rad"], label="bremsstrahlung"); ax[1, 1].plot(x, p["q_wall"], label="radial/wall closure")
    ax[1, 1].set_yscale("symlog", linthresh=1e3); ax[1, 1].set_ylabel("volumetric power [W/m$^3$]"); ax[1, 1].legend(fontsize=8); ax[1, 1].set_title("Energy source/sink terms")
    for a in ax.ravel():
        a.set_xlabel("axial position x [mm]"); a.grid(alpha=.3)
    fig.suptitle(f"LSP result: {res.classification} (converged={res.primary.converged})")
    return _save(fig, outdir, f"{prefix}_profiles.png")


def plot_convergence(res, outdir, prefix="lsp"):
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for k, r in res.branches.items():
        if r.history:
            ax.semilogy([h["it"] for h in r.history], [max(h["res"], 1e-18) for h in r.history], label=f"seed {k}: {'conv' if r.converged else 'NOT conv'}")
    ax.set_xlabel("iteration"); ax.set_ylabel("relative max-norm residual of energy equation"); ax.grid(alpha=.3); ax.legend()
    ax.set_title("Residual history (convergence judged on this, not on peak T)")
    return _save(fig, outdir, f"{prefix}_convergence.png")


def plot_power_chain(ir, outdir):
    labels = [c[0].strip() for c in ir.chain[:8]]
    vals = [c[1] for c in ir.chain[:8]]
    fig, ax = plt.subplots(figsize=(10, 4.8))
    colors = ["tab:blue"] * 6 + ["tab:red", "tab:gray"]
    ax.bar(range(len(vals)), np.array(vals) / 1e3, color=colors)
    ax.set_xticks(range(len(vals))); ax.set_xticklabels([l.replace(" (", "\n(") for l in labels], rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("power [kW]"); ax.set_title("Power chain: delivery (blue) then LSP absorption/residual (red/grey)"); ax.grid(alpha=.3, axis="y")
    return _save(fig, outdir, "integrated_power_chain.png")


def plot_uniform_column_balance(diag, outdir, prefix="lsp"):
    """Deposition vs unavoidable demands for a uniform column at T (ignition/maintenance diagnostic)."""
    T = diag["T"]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.semilogy(T, np.maximum(diag["P_dep"], 1e-30), "r", lw=2, label="laser deposition  P_in(1-exp(-alpha L))")
    ax.semilogy(T, np.maximum(diag["P_loss"], 1e-30), "k", lw=2, label="total demand (enthalpy + radiation + wall)")
    ax.semilogy(T, np.maximum(diag["P_enthalpy"], 1e-30), "b--", label="gas enthalpy flux  mdot (h - h_in)")
    ax.semilogy(T, np.maximum(diag["P_rad"], 1e-30), "g:", label="bremsstrahlung")
    ax.semilogy(T, np.maximum(diag["P_wall"], 1e-30), "m:", label="radial/wall closure")
    ax.set_ylim(1e-3, None); ax.set_xlabel("uniform column temperature T [K]"); ax.set_ylabel("power [W]"); ax.grid(alpha=.3, which="both")
    ax.legend(fontsize=7); ax.set_title("Uniform-column power balance (ignition/maintenance diagnostic)")
    return _save(fig, outdir, f"{prefix}_uniform_column_balance.png")


def make_all(res, outdir, prefix="lsp", diag=None):
    """All LSP figures; returns the list of VERIFIED absolute paths (raises OutputError on any failure)."""
    os.makedirs(outdir, exist_ok=True)
    out = [plot_profiles(res, outdir, prefix), plot_convergence(res, outdir, prefix)]
    if diag is not None:
        out.append(plot_uniform_column_balance(diag, outdir, prefix))
    return out
