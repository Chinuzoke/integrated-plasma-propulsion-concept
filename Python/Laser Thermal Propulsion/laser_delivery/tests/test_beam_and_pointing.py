import math, unittest
import tests._path  # noqa
import numpy as np
import beam_propagation as bp
import pointing_error as pe
import config as cfg


class TestBeam(unittest.TestCase):
    w0, lam, m2 = 0.25, 1064e-9, 11.8

    def test_rayleigh_known_value(self):
        self.assertAlmostEqual(bp.rayleigh_range(self.w0, self.lam, self.m2), math.pi * 0.0625 / (11.8 * 1064e-9), places=6)
        self.assertAlmostEqual(bp.rayleigh_range(self.w0, self.lam, self.m2), 15638.9, delta=0.1)

    def test_divergence_and_bpp(self):
        th = bp.divergence_half_angle(self.w0, self.lam, self.m2)
        self.assertAlmostEqual(th * 1e6, 15.986, places=3)
        self.assertAlmostEqual(self.w0 * th, self.m2 * self.lam / math.pi, places=18)

    def test_radius_at_rayleigh_range(self):
        zr = bp.rayleigh_range(self.w0, self.lam, self.m2)
        self.assertAlmostEqual(bp.beam_radius(zr, self.w0, self.lam, self.m2) / self.w0, math.sqrt(2), places=12)

    def test_baseline_radius(self):
        self.assertAlmostEqual(bp.beam_radius(1.2e6, self.w0, self.lam, self.m2), 19.1846, places=3)
        self.assertAlmostEqual(bp.beam_diameter(1.2e6, self.w0, self.lam, self.m2), 38.369, places=2)

    def test_m2_not_ignored(self):
        self.assertGreater(bp.beam_radius(1.2e6, self.w0, self.lam, 11.8), 10 * bp.beam_radius(1.2e6, self.w0, self.lam, 1.0) * 0.9)

    def test_intensity_integrates_to_power(self):
        w = 3.0
        r = np.linspace(0, 8 * w, 200001)
        P = 2 * np.pi * np.trapezoid(bp.intensity(r, 7.0, w) * r, r)
        self.assertAlmostEqual(P, 7.0, places=6)

    def test_invalid(self):
        with self.assertRaises(ValueError): bp.rayleigh_range(-1, self.lam, 1)
        with self.assertRaises(ValueError): bp.beam_radius(-5.0, self.w0, self.lam, 1)
        with self.assertRaises(ValueError): bp.on_axis_intensity(-1.0, 1.0)


class TestPointing(unittest.TestCase):
    def test_arcsec(self):
        self.assertAlmostEqual(pe.arcsec_to_rad(1.0), 4.848136811e-6, places=14)

    def test_displacement_baseline(self):
        self.assertAlmostEqual(pe.displacement_from_arcsec(1.2e6, 0.1), 0.58178, places=4)

    def test_zero(self):
        self.assertEqual(pe.lateral_displacement(1.2e6, 0.0), 0.0)

    def test_invalid(self):
        with self.assertRaises(ValueError): pe.lateral_displacement(-1.0, 1e-6)
        with self.assertRaises(ValueError): pe.arcsec_to_rad(-0.1)


class TestConfig(unittest.TestCase):
    def test_unresolved_refused(self):
        with self.assertRaises(ValueError): cfg.value("window_thickness")

    def test_invalid_inputs(self):
        base = cfg.SystemInputs.from_config()
        for bad in (dict(power_w=-1), dict(aperture_radius_m=-1), dict(cond_T=1.2), dict(m2=0.9), dict(mirror_R=0)):
            with self.assertRaises(ValueError): base.with_(**bad)

    def test_baseline_values_unchanged(self):
        b = cfg.SystemInputs.from_config()
        self.assertEqual((b.power_w, b.wavelength_m, b.waist_m, b.m2, b.distance_m, b.aperture_radius_m, b.pointing_arcsec),
                         (30000.0, 1064e-9, 0.25, 11.8, 1.2e6, 5.0, 0.1))
