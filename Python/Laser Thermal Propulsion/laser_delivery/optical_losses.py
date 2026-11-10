"""
optical_losses.py -- series optical train and cumulative power balance.

Each component receives the ACTUAL output power of the preceding component
(never the original transmitter power). Conservation is enforced per stage.
"""
from __future__ import annotations

from config import SystemInputs
from optical_components import FresnelWindow, LumpedComponent, StageResult


class PowerBalanceError(AssertionError):
    pass


def assert_conserved(stage: StageResult, rtol: float = 1e-12) -> None:
    if abs(stage.residual) > rtol * stage.p_in:
        raise PowerBalanceError(f"{stage.name}: energy not conserved, residual {stage.residual:.3e} W")
    if stage.p_out < 0 or stage.total_loss < -1e-12 * stage.p_in or stage.p_out > stage.p_in * (1 + 1e-12):
        raise PowerBalanceError(f"{stage.name}: non-physical power channel")


def build_train(inp: SystemInputs) -> list:
    """Mirror -> beam-conditioning assembly -> window (representative catalogue components, register v1)."""
    import config as cfg
    ppm = inp.mirror_abs_scatter_ppm
    mirror_known = {}
    if ppm is not None:
        mirror_known = {"absorption + scattering COMBINED (typical measured, upper bound)": ppm * 1e-6}
        m_basis = "typical measured reflectance + combined abs+scatter (PROVISIONAL manufacturer measurement)"
        m_unres = "remainder (1 - R - combined abs+scatter), e.g. transmission through coating/substrate: unattributed"
    else:
        m_basis = "GUARANTEED manufacturer specification (> 99.8 %), used as conservative baseline"
        m_unres = "absorption, scattering and any transmission combined; individual split UNRESOLVED"
    mirror = LumpedComponent(
        "Receiving mirror", inp.mirror_R, m_basis, m_unres,
        ("representative part: catalogue dielectric mirror specified at 45 deg AOI; its applicability to a "
         "large spacecraft receiving mirror, and angle/polarisation dependence, are not established",),
        mirror_known)
    cond = LumpedComponent(
        "Beam-conditioning assembly", inp.cond_T,
        "manufacturer specification (> 96 %), 0.960 used as conservative assembly-level total",
        "reflection, absorption, scatter combined; per-surface split UNRESOLVED",
        (f"representative: 8x Galilean expander; design-derived {cfg.PARAMS['expander_num_lenses'].value} lenses / "
         f"{cfg.PARAMS['expander_num_surfaces'].value} surfaces (not a published prescription)",
         "loss NOT divided among surfaces; assembly total used once (no double counting)",
         "beam-size change (x8) and clear-aperture matching are not modelled in this power budget"))
    if inp.window_mode == "catalogue":
        window = LumpedComponent(
            "Fused-silica window", inp.window_T,
            "manufacturer specification (at least 99.6 %)",
            "total loss combined; bulk absorption and scattering UNRESOLVED",
            ("manufacturer derives T from single-sided reflectance (< 0.2 %): that reflectance is the basis of the "
             "total and is NOT added as an extra loss",))
    else:
        window = FresnelWindow("Fused-silica window", inp.wavelength_m,
                               inp.window_thickness_m, inp.window_alpha_per_m)
    return [mirror, cond, window]


def run_train(p_in: float, train: list) -> list:
    stages, p = [], p_in
    for comp in train:
        st = comp.propagate(p)
        assert_conserved(st)
        stages.append(st)
        p = st.p_out
    return stages


def cumulative_table(stages: list, p_reference: float) -> list:
    """Rows: name, p_in, p_out, stage_loss, stage_loss_pct, cum_throughput_pct (vs p_reference)."""
    rows = []
    for st in stages:
        rows.append((st.name, st.p_in, st.p_out, st.total_loss,
                     100.0 * st.total_loss / st.p_in, 100.0 * st.p_out / p_reference))
    return rows
