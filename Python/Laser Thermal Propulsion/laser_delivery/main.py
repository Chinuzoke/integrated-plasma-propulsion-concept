"""
main.py -- run the full laser-delivery model: budget, validation, sensitivity, report, plots.

Usage:  python main.py [--outdir outputs] [--no-plots] [--mc N] [--window-mode catalogue|uncoated_fresnel]
All physical inputs live in config.py (PARAMS); none are set here.
"""
from __future__ import annotations

import argparse
import os

import config as cfg
import power_budget as pb
import references
import report
import uncertainty as un
import validation as val


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--mc", type=int, default=2000, help="Monte Carlo samples (0 to skip)")
    ap.add_argument("--window-mode", choices=["catalogue", "uncoated_fresnel"], default=None)
    args = ap.parse_args(argv)

    inp = cfg.SystemInputs.from_config()
    if args.window_mode:
        inp = inp.with_(window_mode=args.window_mode)
    os.makedirs(args.outdir, exist_ok=True)

    bud = pb.compute_budget(inp)
    checks = val.run_all(inp)
    el = un.local_elasticities(inp)
    tor = un.tornado(inp)
    mc = un.monte_carlo(inp, n=args.mc) if args.mc > 0 else None

    text = report.build_report(bud, checks, el, tor, mc)
    print(text)
    with open(os.path.join(args.outdir, "report.txt"), "w", encoding="utf-8") as f:
        f.write(text)
    report.write_budget_csv(bud, os.path.join(args.outdir, "power_budget.csv"))
    references.write_csv(os.path.join(args.outdir, "references.csv"))
    if not args.no_plots:
        import plots
        for p in plots.make_all(bud, args.outdir):
            print("saved", p)
    return 0 if all(c.passed for c in checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
