"""
power_budget.py -- complete series power accounting, transmitter -> LSP chamber entrance.

Stages: transmitter -> vacuum propagation -> receiver aperture capture ->
receiving mirror -> beam-conditioning assembly -> fused-silica window ->
(delivered to LSP chamber entrance). Ends there.
"""
from __future__ import annotations

from dataclasses import dataclass

import beam_propagation as bp
import pointing_error as pe
import receiver as rx
from config import SystemInputs
from optical_components import StageResult
from optical_losses import assert_conserved, build_train, run_train


@dataclass(frozen=True)
class PowerBudget:
    inputs: SystemInputs
    beam: bp.BeamState
    capture: rx.CaptureResult
    # StageResult list, all stages in order (no transmitter row)
    stages: list

    @property
    def delivered_w(self) -> float:
        return self.stages[-1].p_out

    @property
    def average_delivered_w(self):
        """Duty-cycle-weighted delivered power, or None when the duty cycle is UNRESOLVED."""
        d = self.inputs.duty_cycle
        return None if d is None else d * self.delivered_w

    @property
    def total_loss_w(self) -> float:
        return self.inputs.power_w - self.delivered_w

    @property
    def efficiency(self) -> float:
        return self.delivered_w / self.inputs.power_w

    def table_rows(self):
        """name, p_in, p_out, stage_loss, % of stage input lost, % of TX lost here, cumulative eff %."""
        ptx = self.inputs.power_w
        return [(s.name, s.p_in, s.p_out, s.total_loss, 100 * s.total_loss / s.p_in,
                 100 * s.total_loss / ptx, 100 * s.p_out / ptx) for s in self.stages]


def compute_budget(inp: SystemInputs) -> PowerBudget:
    beam = bp.beam_state(inp.power_w, inp.wavelength_m,
                         inp.waist_m, inp.m2, inp.distance_m)
    d = pe.lateral_displacement(inp.distance_m, inp.pointing_rad)

    # Stage 1: vacuum propagation (zero loss by default; hook for future media)
    p_prop_loss = inp.power_w * inp.prop_loss_fraction
    prop = StageResult("Free-space propagation", inp.power_w, inp.power_w - p_prop_loss,
                       {"propagation attenuation": p_prop_loss}, 0.0, "",
                       "vacuum: zero by assumption (no medium model)")
    assert_conserved(prop)

    # Stage 2: capture by the receiver aperture
    cap = rx.capture_power(prop.p_out, beam.beam_radius_m,
                           inp.aperture_radius_m, d)
    capstage = StageResult("Receiver aperture capture", prop.p_out, cap.p_captured_w,
                           {"missed (outside aperture)": cap.p_missed_w}, 0.0, "",
                           "numerical 2-D Gaussian integral over the circular aperture",
                           (f"displacement {d:.4f} m, beam radius {beam.beam_radius_m:.3f} m",))
    assert_conserved(capstage)

    # Stages 3-5: optical train, each fed by the previous output
    optical = run_train(cap.p_captured_w, build_train(inp))
    stages = [prop, capstage] + optical
    for a, b in zip(stages[:-1], stages[1:]):
        if a.p_out != b.p_in:
            raise AssertionError("stages are not chained in series")
    return PowerBudget(inp, beam, cap, stages)
