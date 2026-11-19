import sys, numpy as np
sys.path.insert(0, "src")
from radiation.line_absorption import *

def test_all_ar2_rows_reproduce_ritz_wavelength():
    rows = [r for r in load_lines() if r["ion"] == "Ar II"]
    assert len(rows) == 12
    for r in rows:
        lam = air_wavelength_from_levels(float(r["E_i_cm1"]), float(r["E_k_cm1"]))
        assert abs(lam - float(r["ritz_air_nm"])) < 6e-5

def test_f_and_a_forms_agree():
    lam, gl, gk, f = 1.08e-6, 6.0, 4.0, 0.3
    a = (gl / gk) * f * 2 * np.pi * Q_E**2 / (EPS0 * M_E * C * lam**2)
    assert np.isclose(integrated_cross_section(a_ki=a, g_l=gl, g_k=gk, lam=lam), integrated_cross_section(f=f), rtol=1e-12)
    assert np.isclose(f_from_a(a, gl, gk, lam), f, rtol=1e-12)

def test_classical_integrated_cross_section_value():
    assert np.isclose(integrated_cross_section(f=1.0), 2.654e-6, rtol=1e-3)   # 0.02654 cm^2 Hz

def test_population_rises_with_temperature_and_degeneracy():
    a = lower_level_density_ion(14000, 1e6, 150147.7021, "5/2")
    assert lower_level_density_ion(18000, 1e6, 150147.7021, "5/2") > a
    assert np.isclose(degeneracy("5/2"), 6.0) and np.isclose(degeneracy("1/2"), 2.0)

def test_no_a_values_in_table():
    assert all("A_ki" not in r for r in load_lines())

def test_doppler_width_scale():
    assert 0.01 < doppler_fwhm_nm(15000, 1081.29) < 0.02
