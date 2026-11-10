import math, unittest
import tests._path  # noqa
import receiver as rx


class TestReceiver(unittest.TestCase):
    def test_centred_matches_analytic(self):
        for a, w in [(5, 19.18), (5, 3), (1, 10), (30, 5)]:
            self.assertAlmostEqual(rx.capture_fraction_radial(a, w), 1 - math.exp(-2 * a * a / w / w), places=12)

    def test_displaced_three_methods(self):
        for a, w, d in [(5, 19.18, 0.58), (5, 3, 6), (5, 8, 4), (2, 10, 15)]:
            r, p, ref = (rx.capture_fraction_radial(a, w, d), rx.capture_fraction_polar(a, w, d),
                         rx.capture_fraction_reference(a, w, d))
            self.assertAlmostEqual(r, ref, places=10)
            self.assertAlmostEqual(p, ref, places=10)

    def test_displaced_differs_from_centred_formula(self):
        self.assertLess(rx.capture_fraction_radial(5, 8, 6), rx.capture_fraction_centered_analytic(5, 8))

    def test_convergence(self):
        ref, rows = rx.convergence_study(5, 19.18, 0.58, node_list=(2, 4, 8, 16, 32))
        errs = [r[2] for r in rows]
        self.assertLess(errs[-1], 1e-12)
        self.assertLess(errs[2], errs[0])

    def test_monotone_pointing(self):
        e = [rx.capture_fraction_radial(5, 19.18, d) for d in range(0, 60, 3)]
        self.assertTrue(all(x >= y for x, y in zip(e, e[1:])))

    def test_monotone_aperture(self):
        e = [rx.capture_fraction_radial(a, 19.18, 0.58) for a in (0.5, 1, 2, 5, 10, 30, 100)]
        self.assertTrue(all(x <= y for x, y in zip(e, e[1:])))
        self.assertAlmostEqual(e[-1], 1.0, places=9)

    def test_grid_independence(self):
        base = rx.capture_fraction_radial(5, 19.18, 0.58, 4, 64)
        for p, n in [(1, 24), (2, 32), (10, 100)]:
            self.assertAlmostEqual(rx.capture_fraction_radial(5, 19.18, 0.58, p, n), base, places=10)

    def test_power_split_conserved(self):
        c = rx.capture_power(3000.0, 19.18, 5.0, 0.58)
        self.assertAlmostEqual(c.p_captured_w + c.p_missed_w, 3000.0, places=9)
        self.assertLess(c.quadrature_self_error, 1e-10)
        self.assertLess(c.method_discrepancy, 1e-10)

    def test_jitter_limits(self):
        self.assertAlmostEqual(rx.mean_capture_with_jitter(5, 19.18, 0.0), rx.capture_fraction_centered_analytic(5, 19.18), places=12)
        self.assertLess(rx.mean_capture_with_jitter(5, 3.0, 2.0), rx.capture_fraction_centered_analytic(5, 3.0))

    def test_invalid(self):
        with self.assertRaises(ValueError): rx.capture_fraction_radial(-5, 10)
        with self.assertRaises(ValueError): rx.capture_fraction_radial(5, 0)
        with self.assertRaises(ValueError): rx.capture_fraction_radial(5, 10, -1)
