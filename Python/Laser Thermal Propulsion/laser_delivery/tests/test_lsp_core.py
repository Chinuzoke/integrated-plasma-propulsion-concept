import math, unittest
import tests._path  # noqa
import numpy as np
from lsp import absorption as ab, argon_thermo as at, constants as C, radiation as rd, tables, transport as tr, validation as v


class TestThermo(unittest.TestCase):
    def test_saha_consistency(self):
        T = np.array([4000., 9000., 15000., 25000., 60000.])
        c = at.composition(T, 2e5)
        self.assertLess(np.max(np.abs((c["n_a"] + c["n_1"] + c["n_2"] + c["n_e"]) * C.k_B * T / 2e5 - 1)), 1e-12)
        self.assertLess(np.max(np.abs((c["n_1"] + 2 * c["n_2"]) / c["n_e"] - 1)), 1e-9)

    def test_matches_legacy_brentq_value(self):
        # independent legacy result recorded in the audit: n_e(18 kK, 1 bar)
        import sys, os
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lsp", "legacy_original"))
        import importlib
        pp = importlib.import_module("plasma_properties"); lc = importlib.import_module("constants")
        sys.path.pop(0)
        for T in (8000., 14000., 20000.):
            self.assertAlmostEqual(at.composition(T, 1e5)["n_e"][0] / pp.saha_composition(T, 1e5, lc.GAS_ARGON)["n_e"], 1.0, places=9)
        for m in ("plasma_properties", "constants"):
            sys.modules.pop(m, None)

    def test_no_density_floor(self):
        self.assertLess(at.composition(500.0, 1e5)["n_e"][0], 1e-30)

    def test_enthalpy_limits(self):
        c = at.composition(np.array([300.]), 1e5)
        self.assertAlmostEqual(at.enthalpy(np.array([300.]), c)[0] / (2.5 * C.k_B * 300 / C.AR_ATOM_MASS), 1.0, places=12)
        hot = at.composition(np.array([4e4]), 1e5)
        self.assertGreater(at.enthalpy(np.array([4e4]), hot)[0], 4 * 2.5 * C.k_B * 4e4 / C.AR_ATOM_MASS / 4 + 1e7)

    def test_invalid(self):
        with self.assertRaises(ValueError): at.composition(-5.0, 1e5)
        with self.assertRaises(ValueError): at.composition(5000.0, -1.0)


class TestRadiationAbsorption(unittest.TestCase):
    def test_constants_vs_cgs(self):
        self.assertAlmostEqual(v._cgs_from_si_Cff() / 3.692e8, 1.0, places=3)
        self.assertAlmostEqual(v._cgs_from_si_Cbr() / 1.4e-27, 1.0, delta=0.03)

    def test_dimensions(self):
        self.assertTrue(np.allclose(v.dimension_of_Cff(), v._dim(m=5, s=-3, K=0.5)))

    def test_alpha_scalings(self):
        T, ne, z = 15000.0, 1e23, 1e23
        a = ab.alpha_ff(T, ne, z, 1.064e-6)
        self.assertAlmostEqual(ab.alpha_ff(T, 2 * ne, z, 1.064e-6) / a, 2.0, places=12)       # linear in n_e
        self.assertAlmostEqual(ab.alpha_ff(T, ne, 3 * z, 1.064e-6) / a, 3.0, places=12)       # linear in sum Z^2 n_i
        self.assertAlmostEqual(ab.alpha_ff(T, ne, z, 1.064e-6, gaunt=1.5) / a, 1.5, places=12)

    def test_wavelength_dependence(self):
        c = at.composition(np.array([18000.]), 1e5)
        a1 = ab.alpha_ff(18000., c["n_e"][0], ab.zsq_n_ion(c)[0], 1.064e-6)
        a2 = ab.alpha_ff(18000., c["n_e"][0], ab.zsq_n_ion(c)[0], 10.6e-6)
        self.assertGreater(a2 / a1, 100.0)

    def test_laser_density_cutoff(self):
        self.assertTrue(ab.plasma_frequency_ok(1e25, 1.064e-6))
        self.assertFalse(ab.plasma_frequency_ok(1e28, 1.064e-6))

    def test_unresolved_en_ib(self):
        with self.assertRaises(NotImplementedError): ab.electron_neutral_alpha()

    def test_brems_scales_sqrtT(self):
        self.assertAlmostEqual(rd.bremsstrahlung_loss(40000., 1e23, 1e23) / rd.bremsstrahlung_loss(10000., 1e23, 1e23), 2.0, places=12)

    def test_nec_table_no_extrapolation(self):
        t = rd.NECTable([1e4, 2e4], [1e8, 2e8], 1e5, "test")
        self.assertAlmostEqual(t.loss(1.5e4), 4 * math.pi * 1.5e8)
        with self.assertRaises(ValueError): t.loss(3e4)


class TestTransportTables(unittest.TestCase):
    def test_anchors(self):
        self.assertAlmostEqual(float(tr.viscosity_ar(300.0)) / 2.27e-5, 1.0, delta=0.03)
        self.assertAlmostEqual(float(tr.k_heavy(300.0)) / 0.0177, 1.0, delta=0.05)

    def test_extrapolation_flag(self):
        self.assertTrue(bool(tr.tstar_extrapolated(15000.0))); self.assertFalse(bool(tr.tstar_extrapolated(5000.0)))

    def test_tabulated_conductivity(self):
        k = tr.TabulatedConductivity([300, 1000], [0.02, 0.05], "test")
        self.assertAlmostEqual(float(k(650)), 0.035)
        with self.assertRaises(ValueError): k(2000)

    def test_spline_vs_direct(self):
        tab = tables.PropertyTable(1e6, 1.064e-6)
        T = np.array([700., 5000., 12000., 22000.])
        c = at.composition(T, 1e6)
        self.assertLess(np.max(np.abs(tab.h(T) / at.enthalpy(T, c) - 1)), 1e-6)
        self.assertLess(np.max(np.abs(tab.n_e(T) / c["n_e"] - 1)), 1e-6)

    def test_range_error(self):
        tab = tables.PropertyTable(1e5, 1.064e-6)
        with self.assertRaises(tables.PropertyRangeError): tab.h(np.array([1e7]))
        with self.assertRaises(tables.PropertyRangeError): tab.alpha(np.array([100.]))

    def test_kirchhoff_potential_constant_k(self):
        tab = tables.PropertyTable(1e5, 1.064e-6, k_override=lambda T: np.full_like(T, 0.3))
        self.assertAlmostEqual(float(tab.S(np.array([5000.]), 500.0)[0]), 0.3 * 4500.0, delta=0.3 * 4500 * 1e-5)


if __name__ == "__main__":
    unittest.main()
