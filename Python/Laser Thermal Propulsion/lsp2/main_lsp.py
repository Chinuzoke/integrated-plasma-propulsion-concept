"""
main_lsp.py -- run the LSP subsystem ON ITS OWN.  python -m lsp.main_lsp --power-w W --wavelength-m M [options]

The laser input must be given explicitly (there is no default laser power). `--legacy-case` reproduces the earlier
model's laser input (3000 W at 10.6 um) and operating inputs for comparison purposes only.
For the integrated run use main_integrated.py at the project root.
"""
from __future__ import annotations

import argparse
import os

from . import outputs, report, validation
from .config import LaserInput, LSPConfig
from .run import run_lsp


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--power-w", type=float)
    ap.add_argument("--wavelength-m", type=float)
    ap.add_argument("--legacy-case", action="store_true", help="3000 W, 10.6 um, legacy operating inputs (for comparison only)")
    ap.add_argument("--pressure-pa", type=float); ap.add_argument("--mdot", type=float)
    ap.add_argument("--radius-m", type=float); ap.add_argument("--length-m", type=float)
    ap.add_argument("--cells", type=int, default=200)
    ap.add_argument("--wall-closure", choices=["kirchhoff", "adiabatic", "legacy_newton"], default="kirchhoff")
    ap.add_argument("--outdir", default=None, help="default: <project root>/outputs/lsp (independent of the current directory)"); ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--skip-validation", action="store_true")
    a = ap.parse_args(argv)
    if a.legacy_case:
        laser = LaserInput(3000.0, 10.6e-6, None, "legacy default case (comparison only)", "legacy")
    elif a.power_w and a.wavelength_m:
        laser = LaserInput(a.power_w, a.wavelength_m, None, "explicit command-line input (standalone)", "stated by caller")
    else:
        ap.error("give --power-w and --wavelength-m (or --legacy-case): there is no default laser power")
    kw = {"n_cells": a.cells, "wall_closure": a.wall_closure}
    for src, dst in (("pressure_pa", "pressure_pa"), ("mdot", "mdot_kg_s"), ("radius_m", "plasma_radius_m"), ("length_m", "length_m")):
        if getattr(a, src) is not None:
            kw[dst] = getattr(a, src)
    cfg = LSPConfig(laser=laser, **kw)
    res = run_lsp(cfg)
    checks = None if a.skip_validation else validation.run_all()
    text = report.format_lsp_report(res, checks)
    outdir = outputs.resolve_outdir(a.outdir, "lsp")
    written = []
    try:
        written.append(outputs.write_text(os.path.join(outdir, "lsp_report.txt"), text))
        report.write_lsp_csv(res, os.path.join(outdir, "lsp_budget.csv"), os.path.join(outdir, "lsp_profiles.csv"))
        report.write_register_csv(cfg, os.path.join(outdir, "lsp_parameter_register.csv"))
        written += [os.path.join(outdir, n) for n in ("lsp_budget.csv", "lsp_profiles.csv", "lsp_parameter_register.csv")]
        if not a.no_plots:
            from . import plots
            written += plots.make_all(res, outdir, "lsp", res.diagnostics["uniform_column"])
    except outputs.OutputError as e:
        print(text); print(f"OUTPUT FAILURE: {e}"); return 2
    rows = outputs.verification_table(written)
    text += report.output_verification_section(rows)
    outputs.write_text(os.path.join(outdir, "lsp_report.txt"), text)
    print(text)
    if not all(r["ok"] for r in rows):
        return 2
    return 0 if (checks is None or all(c.passed for c in checks)) else 1


if __name__ == "__main__":
    raise SystemExit(main())
