"""
main_integrated.py -- run BOTH subsystems: laser delivery, then LSP fed by the actual delivered power.

    python main_integrated.py [--pressure-pa P] [--mdot kg/s] [--radius-m R] [--length-m L] [--cells N]
                              [--set delivery_field=value ...] [--demo-ignited] [--scan] [--outdir DIR] [--no-plots]

The delivery-only entry point (main.py) is unchanged. LSP operating inputs default to the legacy
(ORIGINAL_MODEL_DEFAULT, unsourced, illustrative) values; --demo-ignited selects an EXPLORATORY operating point
(R = 0.5 mm, 10 MPa, 1 mg/s) at which this model has an ignited branch for the delivered 1064 nm power. It is not
a design recommendation.
"""
from __future__ import annotations

import argparse
import os

import config as delivery_config
import validation as delivery_validation
from lsp import integration, plots, report, scan, validation as lsp_validation
from lsp.references import write_csv as write_lsp_refs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pressure-pa", type=float); ap.add_argument("--mdot", type=float)
    ap.add_argument("--radius-m", type=float); ap.add_argument("--length-m", type=float)
    ap.add_argument("--cells", type=int, default=200)
    ap.add_argument("--wall-closure", choices=["kirchhoff", "adiabatic", "legacy_newton"], default="kirchhoff")
    ap.add_argument("--set", action="append", default=[], metavar="FIELD=VALUE",
                    help="override a laser-delivery SystemInputs field, e.g. --set aperture_radius_m=10")
    ap.add_argument("--demo-ignited", action="store_true")
    ap.add_argument("--scan", action="store_true", help="exploratory (R, p, mdot) scan at the delivered power")
    ap.add_argument("--outdir", default="outputs/integrated"); ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args(argv)

    dinp = delivery_config.SystemInputs.from_config()
    for item in a.set:
        k, v = item.split("=")
        dinp = dinp.with_(**{k: float(v) if v.replace(".", "", 1).replace("e-", "", 1).replace("e", "", 1).lstrip("-").isdigit() else v})
    kw = {"n_cells": a.cells, "wall_closure": a.wall_closure}
    if a.demo_ignited:
        kw.update(pressure_pa=1.0e7, mdot_kg_s=1e-6, plasma_radius_m=5e-4)
    for src, dst in (("pressure_pa", "pressure_pa"), ("mdot", "mdot_kg_s"), ("radius_m", "plasma_radius_m"), ("length_m", "length_m")):
        if getattr(a, src) is not None:
            kw[dst] = getattr(a, src)

    ir = integration.run_integrated(dinp, kw)
    d_checks = delivery_validation.run_all(dinp)
    l_checks = lsp_validation.run_all()
    text = report.format_integrated_report(ir, l_checks)
    text += f"\n\nLASER-DELIVERY VERIFICATION: {sum(c.passed for c in d_checks)}/{len(d_checks)} checks passed (delivery subsystem unchanged)\n"
    if a.scan:
        rows = scan.exploratory_scan(ir.laser, (1e6, 3e6, 1e7), (1e-7, 1e-6), (2e-3, 1e-3, 5e-4))
        text += "\nEXPLORATORY SCAN at the delivered power (hot-seed solve; NOT a design study; 1-D collimated model):\n" + scan.format_scan(rows) + "\n"
    print(text)
    os.makedirs(a.outdir, exist_ok=True)
    open(os.path.join(a.outdir, "integrated_report.txt"), "w", encoding="utf-8").write(text)
    report.write_lsp_csv(ir.lsp, os.path.join(a.outdir, "lsp_budget.csv"), os.path.join(a.outdir, "lsp_profiles.csv"))
    report.write_register_csv(ir.lsp.cfg, os.path.join(a.outdir, "lsp_parameter_register.csv"))
    write_lsp_refs(os.path.join(a.outdir, "lsp_references.csv"))
    if not a.no_plots:
        for pth in (plots.plot_profiles(ir.lsp, a.outdir, "integrated"), plots.plot_convergence(ir.lsp, a.outdir, "integrated"),
                    plots.plot_power_chain(ir, a.outdir)):
            print("saved", pth)
    ok = all(c.passed for c in d_checks) and all(c.passed for c in l_checks) and all(c.passed for c in ir.checks)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
