import math, unittest
import tests._path  # noqa
from laser_source import LaserSource


class TestLaserSource(unittest.TestCase):
    def test_baseline(self):
        s = LaserSource.from_config()
        self.assertAlmostEqual(s.beam_parameter_product_m_rad, 11.8 * 1064e-9 / math.pi, places=18)
        self.assertAlmostEqual(s.waist_radius_m * s.divergence_half_angle_rad, s.beam_parameter_product_m_rad, places=18)

    def test_invalid(self):
        for bad in (dict(power_w=-1.0), dict(wavelength_m=0.0), dict(waist_radius_m=-0.1), dict(m_squared=0.9), dict(power_w=float("nan"))):
            kw = dict(power_w=1.0, wavelength_m=1064e-9, waist_radius_m=0.25, m_squared=1.0); kw.update(bad)
            with self.assertRaises(ValueError): LaserSource(**kw)
