import math, unittest
import tests._path  # noqa
import config as cfg
import optical_components as oc
import optical_losses as ol
import power_budget as pb
import uncertainty as un


class TestOptics(unittest.TestCase):
    def test_lumped_conservation(self):
        st = oc.LumpedComponent("m", 0.998, "b", "n").propagate(1000.0)
        self.assertAlmostEqual(st.p_out, 998.0)
        self.assertAlmostEqual(st.residual, 0.0, places=12)
        self.assertEqual(st.known_losses, {})          # no invented split
        self.assertAlmostEqual(st.combined_unresolved_loss, 2.0)

    def test_fused_silica_index_1064(self):
        self.assertAlmostEqual(oc.fused_silica_index(1064e-9), 1.4496, places=3)
        self.assertAlmostEqual(oc.fused_silica_index(587.56e-9), 1.4585, places=3)

    def test_uncoated_window_physics(self):
        n = oc.fused_silica_index(1064e-9)
        R = oc.fresnel_reflectance_normal(n)
        self.assertAlmostEqual(R, 0.0337, places=3)
        T, Rt, A = oc.plate_fractions(n, None, None)
        self.assertAlmostEqual(T, (1 - R) / (1 + R), places=12)
        self.assertAlmostEqual(T + Rt + A, 1.0, places=12)
        self.assertLess(T, 0.996)  # uncoated cannot reach the AR-coated catalogue figure

    def test_window_with_absorption(self):
        n = oc.fused_silica_index(1064e-9)
        T0, _, _ = oc.plate_fractions(n, None, None)
        T1, R1, A1 = oc.plate_fractions(n, 0.02, 5.0)
        self.assertLess(T1, T0); self.assertGreater(A1, 0); self.assertAlmostEqual(T1 + R1 + A1, 1, places=12)
        st = oc.FresnelWindow("w", 1064e-9, 0.02, 5.0).propagate(100.0)
        self.assertAlmostEqual(st.residual, 0.0, places=10)

    def test_invalid_components(self):
        with self.assertRaises(ValueError): oc.LumpedComponent("x", 1.01, "", "")
        with self.assertRaises(ValueError): oc.LumpedComponent("x", 0.0, "", "")
        with self.assertRaises(ValueError): oc.fresnel_reflectance_normal(0.5)
        with self.assertRaises(ValueError): oc.fused_silica_index(5e-6)


class TestBudget(unittest.TestCase):
    inp = cfg.SystemInputs.from_config()

    def test_hand_check(self):
        b = pb.compute_budget(self.inp)
        self.assertAlmostEqual(b.capture.efficiency, 0.12681, places=4)
        self.assertAlmostEqual(b.delivered_w, 3630.25, delta=0.05)
        self.assertAlmostEqual(b.efficiency, 0.121008, places=5)

    def test_series_chaining_and_conservation(self):
        b = pb.compute_budget(self.inp)
        for s in b.stages: self.assertLess(abs(s.residual), 1e-9)
        for x, y in zip(b.stages, b.stages[1:]): self.assertEqual(x.p_out, y.p_in)
        self.assertAlmostEqual(b.inputs.power_w, b.delivered_w + b.total_loss_w, places=9)

    def test_optical_train_product(self):
        b = pb.compute_budget(self.inp)
        self.assertAlmostEqual(b.delivered_w / b.capture.p_captured_w, 0.998 * 0.96 * 0.996, places=12)

    def test_bounds(self):
        b = pb.compute_budget(self.inp)
        self.assertTrue(0 < b.efficiency <= 1)
        for s in b.stages: self.assertTrue(0 <= s.p_out <= s.p_in)

    def test_pointing_and_aperture_trends(self):
        p0 = pb.compute_budget(self.inp.with_(pointing_arcsec=0.0)).delivered_w
        p1 = pb.compute_budget(self.inp.with_(pointing_arcsec=1.0)).delivered_w
        p2 = pb.compute_budget(self.inp.with_(pointing_arcsec=10.0)).delivered_w
        self.assertGreater(p0, p1); self.assertGreater(p1, p2)
        self.assertGreater(pb.compute_budget(self.inp.with_(aperture_radius_m=10.0)).delivered_w, p0)

    def test_propagation_loss_hook(self):
        b = pb.compute_budget(self.inp.with_(prop_loss_fraction=0.1))
        self.assertAlmostEqual(b.stages[0].p_out, 27000.0)

    def test_uncoated_window_mode_lower(self):
        a = pb.compute_budget(self.inp).delivered_w
        b = pb.compute_budget(self.inp.with_(window_mode="uncoated_fresnel")).delivered_w
        self.assertLess(b, a)


class TestSensitivity(unittest.TestCase):
    def test_elasticities(self):
        el = un.local_elasticities(cfg.SystemInputs.from_config())
        for k in ("power_w", "mirror_R", "cond_T", "window_T"):
            self.assertAlmostEqual(el[k], 1.0, places=3)
        # far-field scaling P ~ a^2 w0^2 /(M^2 lambda z)^2 -> |S| ~ 2 (slightly less: capture saturates)
        for k, s in (("aperture_radius_m", 1), ("waist_m", 1), ("m2", -1), ("wavelength_m", -1), ("distance_m", -1)):
            self.assertAlmostEqual(el[k], 1.87 * s, delta=0.05)
        self.assertLess(abs(el["pointing_arcsec"]), 0.02)

    def test_tornado_sorted(self):
        rows = un.tornado(cfg.SystemInputs.from_config())
        sw = [r[5] for r in rows]
        self.assertEqual(sw, sorted(sw, reverse=True))

    def test_monte_carlo_reproducible(self):
        i = cfg.SystemInputs.from_config()
        self.assertEqual(un.monte_carlo(i, 50, 1)["mean"], un.monte_carlo(i, 50, 1)["mean"])


class TestValidationModule(unittest.TestCase):
    def test_all_checks_pass(self):
        import validation as v
        bad = [c for c in v.run_all() if not c.passed]
        self.assertEqual(bad, [], msg=str(bad))
