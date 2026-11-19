import sys, numpy as np
sys.path.insert(0, "src")
from radiation.coulomb_approx import *

def dE_h(n1, n2):  # hydrogen transition energy [hartree]
    return 0.5 * (1 / n1**2 - 1 / n2**2)

def test_radial_function_normalised():
    r, u = radial_function(3.0, 2)
    assert np.isclose(np.sum(u**2) * 0.002, 1.0, rtol=1e-6)

def test_hydrogen_2p_3d_f_value():           # NIST: f(2p-3d) = 0.6958
    assert abs(oscillator_strength(2, 1, 3, 2, dE_h(2, 3)) / 0.6958 - 1) < 0.05

def test_hydrogen_2p_4d_f_value():           # NIST: f(2p-4d) = 0.1218
    assert abs(oscillator_strength(2, 1, 4, 2, dE_h(2, 4)) / 0.1218 - 1) < 0.05

def test_hydrogen_3d_4f_f_value():           # NIST: f(3d-4f) = 1.0071
    assert abs(oscillator_strength(3, 2, 4, 3, dE_h(3, 4)) / 1.0071 - 1) < 0.05
