import unittest
import tests._path  # noqa
import config as cfg
import references as ref
import scenarios as sc
import optical_components as oc
import power_budget as pb


class TestRegister(unittest.TestCase):
    def test_every_param_has_matching_record(self):
        reg = {r.parameter: r for r in ref.build_register()}
        for k, p in cfg.PARAMS.items():
            self.assertIn(k, reg)
            self.assertEqual(reg[k].classification, p.status.name)
            self.assertEqual(reg[k].provisional, "yes" if p.provisional else "no")

    def test_classifications(self):
        P, S = cfg.PARAMS, cfg.Status
        self.assertIs(P["transmitter_power"].status, S.PUBLISHED_LITERATURE); self.assertTrue(P["transmitter_power"].provisional)
        self.assertIs(P["m_squared"].status, S.UNVERIFIED)
        self.assertIs(P["pointing_error"].status, S.UNVERIFIED)
        self.assertIs(P["distance"].status, S.ASSUMED_MODELLING_CHOICE)
        self.assertIs(P["receiver_radius"].status, S.ASSUMED_MODELLING_CHOICE)
        self.assertIs(P["mirror_reflectance"].status, S.MANUFACTURER_SPECIFICATION)
        self.assertIs(P["mirror_typical_reflectance"].status, S.PROVISIONAL_MANUFACTURER_MEASUREMENT)
        self.assertIs(P["expander_num_surfaces"].status, S.DESIGN_DERIVED)
        self.assertIs(P["conditioning_transmission"].status, S.MANUFACTURER_SPECIFICATION)

    def test_unresolved_items(self):
        for k in ("cw_average_optical_power", "duty_cycle", "mirror_absorption_fraction", "mirror_scatter_fraction",
                  "expander_surface_losses", "window_absorption_coefficient", "window_scattering", "window_thickness"):
            self.assertIn(k, cfg.unresolved_names())
            with self.assertRaises(ValueError): cfg.value(k)

    def test_no_unresolved_value_used(self):
        self.assertIsNone(cfg.SystemInputs.from_config().duty_cycle)

    def test_recorded_only_flags(self):
        for k in ("expander_expansion_ratio", "mirror_aoi_deg", "window_single_sided_reflectance", "expander_num_surfaces"):
            self.assertFalse(cfg.PARAMS[k].used)


class TestScenarios(unittest.TestCase):
    base = cfg.SystemInputs.from_config()

    def test_reference_reproducible(self):
        self.assertAlmostEqual(pb.compute_budget(self.base).delivered_w, 3630.2530, places=3)
        names = [n for n, _, _ in sc.run_scenarios(self.base)]
        self.assertEqual(names[0], "baseline_v1")

    def test_typical_mirror(self):
        b = {n: x for n, _, x in sc.run_scenarios(self.base)}
        self.assertGreater(b["typical_mirror"].delivered_w, b["baseline_v1"].delivered_w)
        st = b["typical_mirror"].stages[2]
        self.assertAlmostEqual(list(st.known_losses.values())[0] / st.p_in, 150e-6, places=12)
        self.assertAlmostEqual(st.combined_unresolved_loss / st.p_in, 50e-6, places=10)
        self.assertLess(b["uncoated_window_comparison"].delivered_w, b["baseline_v1"].delivered_w)

    def test_sweeps_monotone(self):
        p = [r[3] for r in sc.pointing_sweep(self.base)]
        self.assertTrue(all(x >= y for x, y in zip(p, p[1:])))
        m = [r[3] for r in sc.m2_sweep(self.base)]
        self.assertTrue(all(x >= y for x, y in zip(m, m[1:])))

    def test_duty_cycle_average(self):
        b = pb.compute_budget(self.base.with_(duty_cycle=0.25))
        self.assertAlmostEqual(b.average_delivered_w, 0.25 * b.delivered_w, places=9)
        self.assertIsNone(pb.compute_budget(self.base).average_delivered_w)

    def test_known_channels_validation(self):
        with self.assertRaises(ValueError): oc.LumpedComponent("x", 0.998, "", "", (), {"a": 0.01})
        with self.assertRaises(ValueError): self.base.with_(mirror_abs_scatter_ppm=5000.0)
        with self.assertRaises(ValueError): self.base.with_(duty_cycle=1.5)
