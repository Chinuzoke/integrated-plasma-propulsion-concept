"""report.py -- readable engineering report, power-flow table and CSV export."""
from __future__ import annotations

import csv

import config as cfg
import receiver as rx
import scenarios as sc
import uncertainty as un
import validation as val

LINE = "=" * 96


def _w(x):  # watts, readable
    return f"{x:12.3f}"


def power_flow_table(bud) -> str:
    hdr = f"{'Stage':<30}{'P_in [W]':>13}{'P_out [W]':>13}{'Loss [W]':>12}{'Loss/stage in [%]':>19}{'Cum. eff. [%]':>14}"
    lines = [hdr, "-" * len(hdr)]
    lines.append(f"{'Transmitter':<30}{_w(bud.inputs.power_w):>13}{_w(bud.inputs.power_w):>13}{0.0:>12.3f}{0.0:>19.4f}{100.0:>14.4f}")
    for n, pin, pout, loss, pct, pct_tx, cum in bud.table_rows():
        lines.append(f"{n:<30}{_w(pin):>13}{_w(pout):>13}{loss:>12.3f}{pct:>19.4f}{cum:>14.4f}")
    return "\n".join(lines)


def channel_detail(bud) -> str:
    out = []
    for s in bud.stages:
        out.append(f"[{s.name}]  basis: {s.basis or 'n/a'}")
        for k, v in s.known_losses.items():
            out.append(f"    modelled channel   {k:<44}{v:14.4f} W")
        if s.combined_unresolved_loss > 0:
            out.append(f"    COMBINED loss      {'(split unavailable)':<44}{s.combined_unresolved_loss:14.4f} W")
            out.append(f"        -> {s.unresolved_note}")
        for n in s.notes:
            out.append(f"    note: {n}")
    return "\n".join(out)


def build_report(bud, checks, elasticities, tornado_rows, mc=None) -> str:
    inp, b, c = bud.inputs, bud.beam, bud.capture
    L = [LINE, "LASER-DELIVERY MODEL FOR LSP THERMAL PROPULSION -- ENGINEERING REPORT", LINE, ""]
    L.append("SCOPE: transmitter -> vacuum propagation -> receiver aperture -> mirror -> beam-conditioning ->")
    L.append("       fused-silica window -> LSP chamber entrance. No plasma, nozzle, thrust or trajectory modelling.")
    L.append("")
    L.append("PEAK-POWER REFERENCE SCENARIO: 30 kW is a peak-power value from the cited study; CW/average power and duty")
    L.append("cycle are UNRESOLVED. Results below are NOT a claim of sustained 30 kW propulsion operation, and the")
    L.append("selected components have NOT been validated as a complete propulsion laser system.")
    L.append("")
    L.append("1. INPUTS (SI) and provenance")
    for k, p in cfg.PARAMS.items():
        v = "UNRESOLVED" if p.value is None else p.value
        tag = p.status.name + (" / PROVISIONAL" if p.provisional else "") + ("" if p.used else " / recorded, not used")
        L.append(f"   {k:<32}{str(v):>22} {p.unit:<7} [{tag}]")
        if p.note:
            L.append(f"      - {p.note}")
    L.append("")
    L.append("2. BEAM AT THE RECEIVER PLANE (Gaussian, M^2-corrected)")
    L.append(f"   Rayleigh range            z_R      = {b.rayleigh_range_m:12.2f} m   (z/z_R = {b.z_over_zR:.2f}; far field: {b.far_field})")
    L.append(f"   Far-field half-divergence theta    = {b.divergence_half_angle_rad*1e6:12.3f} urad")
    L.append(f"   Beam radius (1/e^2)       w        = {b.beam_radius_m:12.3f} m")
    L.append(f"   Beam diameter             2w       = {b.beam_diameter_m:12.3f} m")
    L.append(f"   On-axis intensity         I0       = {b.on_axis_intensity_w_m2:12.3f} W/m^2")
    L.append("")
    L.append("3. POINTING AND RECEIVER CAPTURE")
    L.append(f"   Pointing error                     = {inp.pointing_arcsec:g} arcsec = {inp.pointing_rad:.4e} rad")
    L.append(f"   Lateral displacement      d        = {c.displacement_m:12.4f} m")
    L.append(f"   Aperture radius           a        = {c.aperture_radius_m:12.3f} m   (a/w = {c.aperture_radius_m/c.beam_radius_m:.4f})")
    L.append(f"   Captured power                     = {c.p_captured_w:12.3f} W")
    L.append(f"   Capture efficiency                 = {100*c.efficiency:12.4f} %")
    L.append(f"   Power missed by receiver           = {c.p_missed_w:12.3f} W")
    eta_c = rx.capture_fraction_centered_analytic(c.aperture_radius_m, c.beam_radius_m)
    L.append(f"   (centred-beam analytic efficiency  = {100*eta_c:12.4f} %; pointing costs {100*(eta_c-c.efficiency):.4f} percentage points)")
    L.append(f"   Numerical evidence: |eta(64 nodes) - eta(32 nodes)| = {c.quadrature_self_error:.2e};  "
             f"|radial - independent 2-D polar| = {c.method_discrepancy:.2e}")
    eta_j = rx.mean_capture_with_jitter(c.aperture_radius_m, c.beam_radius_m, inp.pointing_rad * inp.distance_m)
    L.append(f"   INFORMATIONAL alternative: if {inp.pointing_arcsec:g} arcsec were RMS jitter PER AXIS (not a fixed offset),")
    L.append(f"   mean capture efficiency = {100*eta_j:.4f} %  (fixed-offset assumption is the baseline).")
    L.append("")
    L.append("4. POWER-FLOW TABLE")
    L.append(power_flow_table(bud))
    L.append("")
    L.append("5. CHANNEL DETAIL (what is modelled vs what is unavailable)")
    L.append(channel_detail(bud))
    L.append("")
    L.append("6. FINAL POWER-BUDGET SUMMARY")
    L.append(f"   Transmitter PEAK optical power     = {inp.power_w:12.3f} W (reference scenario)")
    L.append(f"   Delivered to LSP chamber entrance  = {bud.delivered_w:12.3f} W (peak-power scenario)")
    L.append(f"   Total power lost                   = {bud.total_loss_w:12.3f} W")
    L.append(f"   Overall laser-delivery efficiency  = {100*bud.efficiency:12.4f} %")
    if bud.average_delivered_w is not None:
        L.append(f"   Average delivered power (duty cycle {inp.duty_cycle:g}) = {bud.average_delivered_w:12.3f} W")
    else:
        L.append("   Average delivered power            : not computed (duty cycle UNRESOLVED)")
    opt_loss = sum(s.total_loss for s in bud.stages[2:])
    L.append(f"   of which: missed by receiver {100*c.p_missed_w/inp.power_w:.2f} % of P_tx; optical train {100*opt_loss/inp.power_w:.3f} % of P_tx")
    L.append("   NOTE: baseline optical values are guaranteed manufacturer specifications for representative catalogue parts,")
    L.append("   used as conservative figures. They are not measured performance of an assembled spacecraft system.")
    L.append("")
    L.append("6b. ALTERNATIVE SCENARIOS (reference scenario baseline_v1 is unchanged)")
    for n, d, b2 in sc.run_scenarios(inp):
        L.append(f"   {n:<28} delivered {b2.delivered_w:10.3f} W  ({100*b2.efficiency:.4f} %)   {d}")
    L.append("   POINTING sensitivity cases (pointing has NO established source; values are not specifications):")
    L.append(f"      {'arcsec':>8}{'d [m]':>10}{'capture [%]':>14}{'delivered [W]':>16}")
    for a_, d_, e_, p_ in sc.pointing_sweep(inp):
        L.append(f"      {a_:8.3f}{d_:10.3f}{100*e_:14.4f}{p_:16.3f}")
    L.append("   M^2 sensitivity cases (M^2 has NO established source; 11.8 is an assumed case):")
    L.append(f"      {'M^2':>8}{'w [m]':>10}{'capture [%]':>14}{'delivered [W]':>16}")
    for m_, w_, e_, p_ in sc.m2_sweep(inp):
        L.append(f"      {m_:8.2f}{w_:10.3f}{100*e_:14.4f}{p_:16.3f}")
    L.append("")
    L.append("7. SENSITIVITY (see uncertainty.py: ranges are ILLUSTRATIVE, not measured uncertainties)")
    L.append("   Local elasticity d ln P_del / d ln x at baseline:")
    for k, v in sorted(elasticities.items(), key=lambda kv: -abs(kv[1])):
        L.append(f"      {un.PARAM_LABELS[k]:<34}{v:+8.3f}")
    L.append("   One-at-a-time sweep over illustrative ranges (largest swing first):")
    for n, lo, hi, plo, phi, sw in tornado_rows:
        L.append(f"      {un.PARAM_LABELS[n]:<34} [{lo:.6g}, {hi:.6g}] -> P_del {plo:9.1f} .. {phi:9.1f} W  (swing {sw:8.1f} W)")
    if mc:
        L.append(f"   Monte Carlo (N={mc['n']}, seed={mc['seed']}, uniform over illustrative ranges): "
                 f"mean {mc['mean']:.0f} W, std {mc['std']:.0f} W, 5-95 % [{mc['p05']:.0f}, {mc['p95']:.0f}] W")
        L.append("   Spearman rank correlation with P_del: " + ", ".join(f"{un.PARAM_LABELS[k]} {v:+.2f}" for k, v in list(mc['spearman'].items())[:5]))
    L.append("   Parameter uncertainty is partly explored above; MODEL-FORM uncertainty (Gaussian profile for M^2=11.8,")
    L.append("   scalar M^2, fixed-offset pointing, no jitter/aberration/obscuration) is NOT quantified.")
    L.append("")
    L.append("8. VERIFICATION STATUS")
    L.append(f"   {val.summary(checks)}")
    for cat in ("math", "consistency", "literature", "inputs"):
        sub = [x for x in checks if x.category == cat]
        L.append(f"   [{cat}] {sum(x.passed for x in sub)}/{len(sub)} passed")
    for x in checks:
        if not x.passed:
            L.append(f"   FAILED: {x.name} {x.detail}")
    L.append("   Verified mathematically : equation implementation, quadrature convergence, energy conservation.")
    L.append("   Checked vs literature   : only a spot-check of fused-silica index (limited; see references.py).")
    L.append("   NOT validated           : no experimental/flight data. Gaussian-profile assumption, catalogue optical")
    L.append("                             values, pointing model and all unresolved inputs require validation.")
    L.append("")
    L.append("9. UNRESOLVED INPUTS: " + ", ".join(cfg.unresolved_names()))
    L.append(LINE)
    return "\n".join(L)


def write_budget_csv(bud, path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["stage", "P_in_W", "P_out_W", "known_losses_W", "combined_unresolved_loss_W",
                    "total_loss_W", "loss_pct_of_stage_input", "loss_pct_of_transmitter", "cumulative_efficiency_pct"])
        ptx = bud.inputs.power_w
        for s in bud.stages:
            w.writerow([s.name, s.p_in, s.p_out, sum(s.known_losses.values()), s.combined_unresolved_loss,
                        s.total_loss, 100 * s.total_loss / s.p_in, 100 * s.total_loss / ptx, 100 * s.p_out / ptx])
