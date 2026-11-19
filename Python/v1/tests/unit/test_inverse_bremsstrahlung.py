import sys, numpy as np
sys.path.insert(0, "src")
from radiation.inverse_bremsstrahlung import *
from thermodynamics.argon_saha import plasma_frequency

def test_matches_cgs_textbook_coefficient():
    # cgs constants: e=4.80320e-10 esu, m=9.10938e-28 g, h=6.62607e-27, c=2.99792e10, k=1.38065e-16
    e, m, h, c, k = 4.80320e-10, 9.10938e-28, 6.62607e-27, 2.99792458e10, 1.380649e-16
    a_cgs = 4 * e**6 / (3 * m * h * c) * np.sqrt(2 * np.pi / (3 * k * m))
    assert abs(a_cgs / 3.7e8 - 1) < 0.02                       # textbook value ~3.7e8
    # SI implementation must reproduce the cgs number for a test state (convert units)
    t, ne, ni, nu = 12000.0, 1e17, 1e17, 2.776e14               # cm^-3
    a_ref_cm = a_cgs * ne * ni / (np.sqrt(t) * nu**3) * (1 - np.exp(-h * nu / (k * t)))
    lam = C / nu
    a_si = alpha_ib(t, ne * 1e6, ni * 1e6, lam)                 # m^-1
    assert np.isclose(a_si / 100.0, a_ref_cm, rtol=2e-3)

def test_scales_with_density_squared_and_charge_squared():
    a = alpha_ib(12000, 1e23, 1e23, 1.08e-6)
    assert np.isclose(alpha_ib(12000, 2e23, 2e23, 1.08e-6) / a, 4.0)
    assert np.isclose(alpha_ib(12000, 1e23, 1e23, 1.08e-6, z=2) / a, 4.0)

def test_longer_wavelength_absorbs_more_per_free_electron_pair():
    assert alpha_ib(12000, 1e23, 1e23, 10.6e-6) > 100 * alpha_ib(12000, 1e23, 1e23, 1.08e-6)

def test_low_frequency_limit_nu_to_minus_two():
    # for h nu << kT: (1 - e^-x) ~ x, so alpha ~ nu^-2
    r = alpha_ib(1e6, 1e20, 1e20, 20e-6) / alpha_ib(1e6, 1e20, 1e20, 10e-6)
    assert np.isclose(r, 4.0, rtol=0.01)

def test_laser_frequency_above_plasma_frequency():
    nu = C / 1.08e-6
    assert 2 * np.pi * nu > plasma_frequency(1.6e24)            # no cutoff at the highest ne reached at 1 MPa

def test_fractional_absorption_limits():
    assert fractional_absorption(0.0, 1.0) == 0.0
    assert np.isclose(fractional_absorption(1e12, 1.0), 1.0)
