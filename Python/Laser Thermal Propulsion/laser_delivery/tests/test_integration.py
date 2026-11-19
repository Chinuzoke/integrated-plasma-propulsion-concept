import os, re, subprocess, sys, tempfile, unittest
import tests._path  # noqa
import config as dcfg
import power_budget as dpb
from lsp import integration as ig, config as lc, run as lr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEMO = dict(pressure_pa=1.0e7, mdot_kg_s=1e-6, plasma_radius_m=5e-4, n_cells=100)


class TestInterface(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ir = ig.run_integrated(lsp_kwargs=DEMO)

    def test_lsp_receives_actual_delivery_output(self):
        b = dpb.compute_budget(dcfg.SystemInputs.from_config())
        self.assertEqual(self.ir.laser.power_w, b.delivered_w)
        self.assertEqual(self.ir.lsp.cfg.laser.power_w, b.delivered_w)
        self.assertIn("laser-delivery subsystem", self.ir.laser.source)

    def test_units_and_wavelength_consistency(self):
        self.assertEqual(self.ir.laser.wavelength_m, self.ir.budget.inputs.wavelength_m)
        self.assertEqual(self.ir.laser.m2, self.ir.budget.inputs.m2)
        self.assertTrue(1e-7 < self.ir.laser.wavelength_m < 1e-4)          # metres, not nm
        self.assertTrue(0 < self.ir.laser.power_w < 1e5)                    # watts, not kW

    def test_not_transmitter_or_captured_power(self):
        a = self.ir.lsp.accounting["P_in"]
        self.assertNotEqual(a, self.ir.budget.inputs.power_w); self.assertNotEqual(a, self.ir.budget.stages[1].p_out)
        self.assertLess(a, self.ir.budget.stages[1].p_out)

    def test_optical_losses_not_counted_twice(self):
        standalone = lr.run_lsp(lc.LSPConfig(laser=lc.LaserInput(self.ir.budget.delivered_w, self.ir.budget.inputs.wavelength_m, None, "explicit", "n/a"), **DEMO))
        for k in ("P_in", "P_abs", "P_res"):
            self.assertAlmostEqual(standalone.accounting[k], self.ir.lsp.accounting[k], places=6)
        self.assertEqual(self.ir.lsp.accounting["P_in"], self.ir.budget.stages[-1].p_out)

    def test_absorption_is_not_transmission(self):
        a = self.ir.lsp.accounting
        self.assertTrue(0 < a["absorption_efficiency"] < 1)
        self.assertNotAlmostEqual(a["absorption_efficiency"], self.ir.budget.efficiency, places=2)
        higher_p = lr.run_lsp(lc.LSPConfig(laser=self.ir.laser, **{**DEMO, "pressure_pa": 3e6}), seeds=("hot",))
        self.assertNotAlmostEqual(higher_p.accounting["absorption_efficiency"], a["absorption_efficiency"], places=2)
        self.assertEqual(higher_p.cfg.laser.power_w, self.ir.laser.power_w)        # LSP choices never change delivery

    def test_energy_balance_consistent(self):
        a = self.ir.lsp.accounting
        self.assertLess(a["energy_balance_residual_rel"], 1e-9)
        self.assertAlmostEqual(a["P_in"], a["P_abs"] + a["P_res"], places=9)
        self.assertTrue(all(c.passed for c in self.ir.checks), msg=str([c for c in self.ir.checks if not c.passed]))

    def test_power_chain_is_series(self):
        ch = {c[0]: c[1] for c in self.ir.chain}
        vals = [c[1] for c in self.ir.chain[:6]]
        self.assertTrue(all(x >= y for x, y in zip(vals, vals[1:])))
        self.assertEqual(vals[-1], self.ir.lsp.accounting["P_in"])


class TestPropagationOfDeliveryChanges(unittest.TestCase):
    def test_aperture_change_propagates(self):
        base = dcfg.SystemInputs.from_config()
        small = ig.run_integrated(base, DEMO)
        large = ig.run_integrated(base.with_(aperture_radius_m=10.0), DEMO)
        self.assertGreater(large.budget.delivered_w, small.budget.delivered_w)
        self.assertEqual(large.lsp.accounting["P_in"], large.budget.delivered_w)
        self.assertGreater(large.lsp.accounting["P_abs"], small.lsp.accounting["P_abs"])

    def test_mirror_reflectance_scales_lsp_input_exactly(self):
        base = dcfg.SystemInputs.from_config()
        a = ig.run_integrated(base, DEMO); b = ig.run_integrated(base.with_(mirror_R=0.99), DEMO)
        self.assertAlmostEqual(b.lsp.accounting["P_in"] / a.lsp.accounting["P_in"], 0.99 / 0.998, places=12)

    def test_wavelength_change_propagates(self):
        base = dcfg.SystemInputs.from_config()
        r = ig.run_integrated(base.with_(wavelength_m=1.07e-6), DEMO)
        self.assertEqual(r.lsp.cfg.laser.wavelength_m, 1.07e-6)

    def test_no_delivered_power_is_an_error(self):
        class B:  # a budget that delivers nothing
            delivered_w = 0.0
            class inputs: wavelength_m = 1e-6; m2 = 1.0; duty_cycle = None
        with self.assertRaises(ValueError): ig.laser_input_from_budget(B)


class TestIndependenceAndEntryPoint(unittest.TestCase):
    def _py(self, code):
        return subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True)

    def test_lsp_runs_without_delivery_modules(self):
        r = self._py("import sys, lsp.run, lsp.validation; assert 'power_budget' not in sys.modules and 'optical_losses' not in sys.modules; print('ok')")
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_delivery_runs_without_lsp(self):
        r = self._py("import sys, config, power_budget; power_budget.compute_budget(config.SystemInputs.from_config()); assert 'lsp' not in sys.modules; print('ok')")
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_no_hardcoded_delivery_power_in_lsp_calculation_modules(self):
        for f in ("run.py", "solver.py", "energy.py", "physics.py", "config.py", "integration.py", "tables.py", "absorption.py"):
            with open(os.path.join(ROOT, "lsp", f), encoding="utf-8") as fh:
                src = fh.read()
            code = re.sub(r"(\"\"\"[\s\S]*?\"\"\"|#.*)", "", src)
            self.assertNotRegex(code, r"3630|3804|30000|30_000", msg=f)

    def test_standalone_lsp_cli_requires_explicit_laser(self):
        r = subprocess.run([sys.executable, "-m", "lsp.main_lsp", "--skip-validation", "--no-plots"], cwd=ROOT, capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0); self.assertIn("no default laser power", r.stderr)

    def test_main_entry_points_run_both_subsystems(self):
        with tempfile.TemporaryDirectory() as d:
            r = subprocess.run([sys.executable, "main_integrated.py", "--no-plots", "--demo-ignited", "--cells", "100", "--outdir", d],
                               cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stdout[-1500:] + r.stderr[-1500:])
            self.assertIn("INTEGRATED LASER-DELIVERY + LSP REPORT", r.stdout)
            self.assertIn("LASER-DELIVERY VERIFICATION", r.stdout)
            self.assertTrue(os.path.exists(os.path.join(d, "integrated_report.txt")))
        r = subprocess.run([sys.executable, "main.py", "--no-plots", "--mc", "0", "--outdir", tempfile.mkdtemp()], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr[-800:])      # original delivery entry point still works unchanged


if __name__ == "__main__":
    unittest.main()
