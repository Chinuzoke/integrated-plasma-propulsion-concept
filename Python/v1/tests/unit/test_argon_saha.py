import sys, numpy as np
sys.path.insert(0, "src")
from thermodynamics.argon_saha import *

def test_low_T_unionised_and_high_T_ionised():
    assert ionisation_fraction(1800.0, 1e6) < 1e-20
    assert ionisation_fraction(60000.0, 1e5) > 0.99

def test_x_monotonic_in_T_and_decreases_with_p():
    ts = np.linspace(5000, 25000, 50)
    x = ionisation_fraction(ts, 1e6)
    assert np.all(np.diff(x) > 0)
    assert ionisation_fraction(15000.0, 1e6) < ionisation_fraction(15000.0, 1e5)

def test_closed_form_satisfies_saha_directly():
    t, p = 16000.0, 5e5
    d = number_densities(t, p)
    assert np.isclose(d["ne"] * d["ni"] / d["n0"], saha_constant(t), rtol=1e-10)

def test_pressure_and_charge_neutrality():
    t, p = 14000.0, 5e5
    d = number_densities(t, p)
    assert np.isclose((d["n0"] + d["ni"] + d["ne"]) * K_B * t, p, rtol=1e-12)
    assert np.isclose(d["ne"], d["ni"])

def test_enthalpy_limits_and_ionisation_contribution():
    t, p = 1800.0, 1e6
    assert np.isclose(enthalpy(t, p), 2.5 * K_B * t / M_AR, rtol=1e-12)
    t = 20000.0
    x = ionisation_fraction(t, 1e6)
    assert enthalpy(t, 1e6) > 2.5 * K_B * t / M_AR * (1 + x) and x * E_ION / M_AR > 0

def test_cp_peaks_near_ionisation_and_exceeds_frozen():
    ts = np.linspace(8000, 24000, 161)
    cp = np.array([cp_equilibrium(t, 1e5) for t in ts])
    assert cp.max() > 3 * 2.5 * K_B / M_AR

def test_dimension_scales():
    assert np.isclose(plasma_frequency(1e20) / (2*np.pi), 8.98e10, rtol=0.01)   # f_p = 8.98 sqrt(ne) Hz
    assert np.isclose(debye_length(1e4, 1e20), 7430*np.sqrt(1e4*K_B/Q_E/1e20), rtol=0.01)   # lambda_D = 7430 sqrt(T_eV/ne) m
