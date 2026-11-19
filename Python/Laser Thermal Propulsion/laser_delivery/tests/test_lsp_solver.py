import unittest
import tests._path  # noqa
import numpy as np
from lsp import config as lc, energy as en, physics as ph, solver as sv, validation as v, run as lr


def cfg(n=120, **kw):
    return v._ignited_cfg(n=n, **kw)


class TestEnergyEquation(unittest.TestCase):
    def test_jacobian_matches_fd(self):
        c = cfg(n=30, wavelength=10.6e-6, p=1e5, R=2e-3)
        phys, geo = ph.build_physics(c), sv.make_geometry(c)
        T = np.linspace(700, 17000, geo.n)
        J = en.evaluate(T, phys, geo, jacobian=True)["J"]
        fd = np.zeros_like(J)
        for j in range(geo.n):
            d = 1e-3 * T[j]; Tp, Tm = T.copy(), T.copy(); Tp[j] += d; Tm[j] -= d
            fd[:, j] = (en.evaluate(Tp, phys, geo)["G"] - en.evaluate(Tm, phys, geo)["G"]) / (2 * d)
        self.assertLess(np.max(np.abs(J - fd)) / np.max(np.abs(fd)), 1e-3)

    def test_beer_lambert_exact_and_conservation(self):
        class UA(ph.ConstantTestPhysics):
            def alpha(self, T): return np.full_like(np.asarray(T, float), 25.0)
        geo = en.Geometry(0.02, 1e-5, 80, 1e-5, 1e8, 300.0)
        f = en.evaluate(np.full(80, 900.0), UA(), geo)
        self.assertAlmostEqual(f["I"][-1] / (1e8 * np.exp(-25 * 0.02)), 1.0, places=12)
        self.assertAlmostEqual(np.sum(f["q_laser"]) * geo.dx, 1e8 - f["I"][-1], delta=1e-4)

    def test_analytic_convection_diffusion(self):
        errs = []
        for n in (50, 100, 200):
            res, g, G = v._solve_linear(n)
            self.assertTrue(res.converged)
            Ta = v._analytic_conv_diff(g.x, g.length, G, 520.0, 0.5, 6.24e6, 300.0)
            errs.append(np.max(np.abs(res.T - Ta)) / (Ta.max() - 300.0))
        self.assertLess(errs[-1], 3e-3)
        self.assertGreater(errs[0] / errs[1], 1.7); self.assertGreater(errs[1] / errs[2], 1.7)   # first order

    def test_pure_convection_exact(self):
        res, g, G = v._solve_linear(80, k=0.0, G=0.24)
        Tex = 300.0 + 6.24e6 * (g.x + 0.5 * g.dx) / (G * 520.0)
        self.assertLess(np.max(np.abs(res.T - Tex)) / (Tex.max() - 300.0), 1e-10)

    def test_no_heating_stays_cold(self):
        c = cfg(wavelength=1.064e-6, p=1e5, R=2e-3, mdot=5e-5)
        res = sv.solve(ph.build_physics(c), sv.make_geometry(c), sv.seed(sv.make_geometry(c), "cold"))
        self.assertTrue(res.converged); self.assertAlmostEqual(res.T.max(), c.T_inlet_K, places=6)
        self.assertEqual(res.balance["P_abs"], 0.0)
        self.assertAlmostEqual(res.balance["P_res"], c.laser.power_w, places=9)


class TestSolverBehaviour(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c = cfg(n=120)
        cls.r = lr.run_lsp(cls.c, seeds=("hot", "cold"))

    def test_ignited_branch_residual_converged(self):
        h = self.r.branches["hot"]
        self.assertTrue(h.converged); self.assertGreater(h.T.max(), 10000)
        self.assertLess(h.residual, 1e-9); self.assertLess(h.balance_residual, 1e-9)

    def test_bistability(self):
        self.assertLess(self.r.branches["cold"].T.max(), 1000); self.assertGreater(self.r.branches["hot"].T.max(), 10000)

    def test_energy_and_beam_balance(self):
        a = self.r.accounting
        self.assertAlmostEqual(a["P_in"], a["P_abs"] + a["P_res"], places=9)
        self.assertLess(abs(a["energy_balance_residual_W"]), 1e-6)
        self.assertAlmostEqual(a["P_abs"], a["P_rad"] + a["P_wall"] + a["dH"], delta=1e-6)

    def test_bounds(self):
        p = self.r.profiles
        self.assertGreaterEqual(p["T"].min(), self.c.T_inlet_K * (1 - 1e-6)); self.assertTrue(np.all(p["n_e"] >= 0))
        self.assertTrue(np.all(p["ion_fraction"] <= 2.0))
        self.assertTrue(np.all(np.diff(p["I"]) <= 0))          # beam never gains power

    def test_ignited_plasma_flagged_as_inlet_attached(self):
        self.assertTrue(any("ATTACHED TO THE INLET" in w for w in self.r.warnings))

    def test_not_converged_is_reported(self):
        res = sv.solve(ph.build_physics(self.c), sv.make_geometry(self.c), sv.seed(sv.make_geometry(self.c), "hot"),
                       sv.SolverOptions(max_iter=3))
        self.assertFalse(res.converged); self.assertIn("max iterations", res.reason)

    def test_out_of_range_initial_state_fails_clearly(self):
        res = sv.solve(ph.build_physics(self.c), sv.make_geometry(self.c), np.full(self.c.n_cells, 2.0e5))
        self.assertFalse(res.converged); self.assertIn("outside property range", res.reason)

    def test_grid_independence(self):
        pa = []
        for n in (60, 120, 240):
            c = cfg(n=n); r = sv.solve(ph.build_physics(c), sv.make_geometry(c), sv.seed(sv.make_geometry(c), "hot"))
            self.assertTrue(r.converged); pa.append(r.balance["P_abs"])
        self.assertLess(abs(pa[2] - pa[1]), abs(pa[1] - pa[0]) + 1e-9)
        self.assertLess(abs(pa[2] - pa[1]) / pa[2], 0.05)

    def test_more_power_absorbs_more(self):
        p1 = lr.run_lsp(cfg(n=100, power=3000.0), seeds=("hot",)).accounting["P_abs"]
        p2 = lr.run_lsp(cfg(n=100, power=4500.0), seeds=("hot",)).accounting["P_abs"]
        self.assertGreater(p2, p1)

    def test_shorter_wavelength_absorbs_less_in_same_state(self):
        # a 10.6 um beam is absorbed far more strongly than 1.064 um at the same plasma state
        r106 = lr.run_lsp(cfg(n=100, wavelength=10.6e-6), seeds=("hot",))
        self.assertGreater(r106.accounting["absorption_efficiency"], self.r.accounting["absorption_efficiency"])


class TestLegacyAudit(unittest.TestCase):
    def test_legacy_audit_facts(self):
        checks = v._legacy_audit_checks()
        self.assertTrue(all(c.passed for c in checks), msg=str([c for c in checks if not c.passed]))


class TestConfig(unittest.TestCase):
    def test_laser_input_required_and_validated(self):
        with self.assertRaises(TypeError): lc.LSPConfig()
        with self.assertRaises(ValueError): lc.LaserInput(-1.0, 1e-6)
        with self.assertRaises(ValueError): lc.LaserInput(10.0, 0.0)

    def test_provenance_register_flags_unresolved(self):
        reg = lc.parameter_register(cfg())
        self.assertTrue(any(r[3] is lc.Prov.UNRESOLVED for r in reg))
        self.assertTrue(all(r[3] is lc.Prov.ORIGINAL_MODEL_DEFAULT for r in reg if r[0] in ("pressure_pa", "mdot_kg_s", "plasma_radius_m", "length_m")))

    def test_beam_consistency_uses_m2(self):
        c = lc.LSPConfig(laser=lc.LaserInput(1000.0, 1.064e-6, 11.8))
        bc = c.beam_consistency(); self.assertTrue(bc["checked"]); self.assertAlmostEqual(bc["rayleigh_range_m"], np.pi * 4e-6 / (11.8 * 1.064e-6), places=9)


if __name__ == "__main__":
    unittest.main()
