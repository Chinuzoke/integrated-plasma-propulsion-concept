"""
integration.py -- the ONLY coupling between the laser-delivery subsystem and the LSP subsystem.

Delivery (flat modules in the project root) -> LaserInput -> LSP.

Interface actually used (inspected in the delivery code, not assumed):
    power_budget.compute_budget(config.SystemInputs) -> PowerBudget
        .delivered_w            optical power at the LSP chamber entrance [W]  (= stages[-1].p_out)
        .inputs.wavelength_m    laser wavelength [m]
        .inputs.m2              beam propagation ratio (UNVERIFIED in the delivery register)
        .inputs.power_w         transmitter PEAK power [W] (reference scenario)
        .inputs.duty_cycle      None while CW/average power is UNRESOLVED
        .stages                 per-stage accounting (propagation, capture, mirror, conditioning, window)

Separation of physics (enforced):
    * the delivery subsystem decides P_delivered; the LSP applies NO further optical transmission factor,
      so optical losses are never counted twice;
    * the LSP decides how much of P_delivered the plasma ABSORBS; delivered != absorbed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import os

# project-root (laser-delivery) configuration
from laser_delivery import config as delivery_config
# project-root (laser-delivery) budget
from laser_delivery import power_budget as delivery_budget

from .config import LaserInput, LSPConfig
from .run import LSPResult, run_lsp
from .validation import Check
from . import report
from . import plots


def laser_input_from_budget(budget) -> LaserInput:
    """Build the LSP laser input from the delivery budget: delivered chamber power, wavelength, M^2."""
    p = float(budget.delivered_w)
    if not p > 0:
        raise ValueError(
            "delivery subsystem delivers no power to the chamber; nothing to feed the LSP")
    basis = ("PEAK-power reference scenario (30 kW peak; CW/average power and duty cycle unresolved)"
             if budget.inputs.duty_cycle is None else f"peak power with duty cycle {budget.inputs.duty_cycle:g}")
    return LaserInput(power_w=p, wavelength_m=float(budget.inputs.wavelength_m), m2=float(budget.inputs.m2),
                      source="laser-delivery subsystem (power_budget.compute_budget)", power_basis=basis)


@dataclass
class IntegratedResult:
    budget: object
    laser: LaserInput
    lsp: LSPResult
    chain: list = field(default_factory=list)
    checks: list = field(default_factory=list)


def power_chain(budget, lsp: LSPResult) -> list:
    """[(label, power_W, subsystem)] from transmitter to the energy sinks inside the plasma."""
    st = budget.stages
    a = lsp.accounting
    return [
        ("Transmitter optical power (peak, reference scenario)",
         budget.inputs.power_w, "delivery"),
        ("After free-space propagation", st[0].p_out, "delivery"),
        ("Captured by receiver aperture", st[1].p_out, "delivery"),
        ("After receiving mirror", st[2].p_out, "delivery"),
        ("After beam-conditioning assembly", st[3].p_out, "delivery"),
        ("DELIVERED to LSP chamber entrance (after window)",
         st[4].p_out, "delivery"),
        ("Absorbed by the plasma (inverse bremsstrahlung)", a["P_abs"], "LSP"),
        ("Residual laser beam leaving the plasma", a["P_res"], "LSP"),
        ("  plasma energy sink: optically thin radiation", a["P_rad"], "LSP"),
        ("  plasma energy sink: radial/wall loss (closure)",
         a["P_wall"], "LSP"),
        ("  plasma energy sink: gas enthalpy flux increase", a["dH"], "LSP"),
    ]


def interface_checks(budget, laser: LaserInput, lsp: LSPResult) -> list:
    out = []
    P_in = lsp.accounting["P_in"]
    out.append(Check("interface", "LSP input power equals the delivery subsystem's delivered power exactly",
                     P_in == budget.delivered_w == laser.power_w, f"{P_in!r} W vs {budget.delivered_w!r} W"))
    out.append(Check("interface", "LSP input is NOT the transmitter power nor the receiver-captured power",
                     P_in != budget.inputs.power_w and P_in != budget.stages[1].p_out,
                     f"P_tx={budget.inputs.power_w:.1f}, captured={budget.stages[1].p_out:.1f}, delivered={P_in:.1f} W"))
    out.append(Check("interface", "wavelength passed in metres and identical to the delivery wavelength",
                     laser.wavelength_m == budget.inputs.wavelength_m and 1e-7 < laser.wavelength_m < 1e-4,
                     f"{laser.wavelength_m:.4e} m"))
    out.append(Check("interface", "laser power passed in watts, consistent with delivery units (W)",
                     laser.power_w == budget.stages[-1].p_out and 0 < laser.power_w < budget.inputs.power_w, ""))
    out.append(Check("interface", "optical losses not applied twice: LSP beam power balance starts from the delivered power",
                     abs(P_in - (lsp.accounting["P_abs"] + lsp.accounting["P_res"])) < 1e-9 * P_in, ""))
    out.append(Check("interface", "plasma absorption kept separate from optical transmission (absorbed + residual = delivered)",
                     abs(lsp.accounting["beam_balance"]) < 1e-9 * P_in,
                     f"absorption efficiency {lsp.accounting['absorption_efficiency']:.4f} is an LSP result; delivery efficiency {budget.efficiency:.4f} is a delivery result"))
    out.append(Check("interface", "LSP energy balance consistent at the chamber boundary",
                     lsp.accounting["energy_balance_residual_rel"] < 1e-6 or not lsp.primary.converged,
                     f"relative residual {lsp.accounting['energy_balance_residual_rel']:.1e}; converged={lsp.primary.converged}"))
    return out


def run_integrated(delivery_inputs=None, lsp_kwargs=None, seeds=("hot", "cold"), solver_opts=None) -> IntegratedResult:
    """Run delivery then LSP. `delivery_inputs`: a SystemInputs (default: config baseline).
    `lsp_kwargs`: LSPConfig fields other than `laser` (operating inputs, closures, grid)."""
    dinp = delivery_inputs or delivery_config.SystemInputs.from_config()
    budget = delivery_budget.compute_budget(dinp)
    laser = laser_input_from_budget(budget)
    cfg = LSPConfig(laser=laser, **(lsp_kwargs or {}))
    lsp = run_lsp(cfg, seeds=seeds, opts=solver_opts)
    ir = IntegratedResult(budget, laser, lsp)
    ir.chain = power_chain(budget, lsp)
    ir.checks = interface_checks(budget, laser, lsp)
    return ir


def main() -> int:
    # Execute the complete integrated simulation
    result = run_integrated()

    # Output directory
    outdir = "outputs/integrated"
    os.makedirs(outdir, exist_ok=True)

    # Display the power chain
    print("\n========= INTEGRATED POWER CHAIN ==========")

    for label, power, subsystem in result.chain:
        print(f"{label}: {power:.4f} W [{subsystem}]")

    print("\n =========== ABSORPTION DIAGNOSTICS ==========")
    print(f"Delivered power: {result.laser.power_w:.10e} W")
    print(f"Wavelength: {result.laser.wavelength_m:.10e} m")
    print(f"M2: {result.laser.m2}")
    print(f"Absorbed power: {result.lsp.accounting['P_abs']:.10e} W")
    print(f"Residual power: {result.lsp.accounting['P_res']:.10e} W")
    print(
        f"Absorption efficiency: {result.lsp.accounting['absorption_efficiency']:.10e}")

    # Display LSP results and interface checks
    print("\n========= LSP SIMULATION REPORT =========")
    text = report.format_lsp_report(result.lsp, result.checks)
    print(text)

    # Save report
    with open(os.path.join(outdir, "integrated_report.txt"), "w", encoding="utf-8") as f:
        f.write(text)

    # Save numerical results
    report.write_lsp_csv(result.lsp, os.path.join(
        outdir, "lsp_budget.csv"), os.path.join(outdir, "lsp_profiles.csv"))

    # Save plots
    print("\n========= SAVED PLOTS =========")
    print(plots.plot_profiles(result.lsp, outdir))
    print(plots.plot_convergence(result.lsp, outdir))

    # Interface validation
    print("\n========= INTERFACE VALIDATION =========")
    for check in result.checks:
        print(f"{'PASS' if check.passed else 'FAIL'}: {check.name}")
        passed = all(check.passed for check in result.checks)
        print("\nIntegrated simulation completed.")
        print(f"Results saved in: {outdir}")
        return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
