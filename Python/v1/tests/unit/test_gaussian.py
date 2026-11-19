import sys, numpy as np
sys.path.insert(0, "src")
from laser_delivery.gaussian import *

def test_rayleigh_known_value():
    assert np.isclose(rayleigh_range(1e-3, 1.064e-6), np.pi*1e-6/1.064e-6)

def test_radius_at_zR_is_sqrt2_w0():
    w0, lam = 0.01, 1.08e-6
    assert np.isclose(beam_radius(rayleigh_range(w0, lam), w0, lam), np.sqrt(2)*w0)

def test_capture_limits():
    assert np.isclose(capture_fraction_centred(1e9, 1.0), 1.0)
    assert np.isclose(capture_fraction_centred(1.0, 1.0), 1-np.exp(-2))

def test_offset_matches_centred_at_zero():
    assert np.isclose(capture_fraction_offset(0.7, 1.0, 0.0), capture_fraction_centred(0.7, 1.0), rtol=1e-4)

def test_offset_reduces_capture():
    assert capture_fraction_offset(0.7, 1.0, 0.8) < capture_fraction_centred(0.7, 1.0)

def test_energy_balance_closes():
    b = power_budget(1000, {"lens": 0.99, "window": 0.995}, 0.005, 1.08e-6, 1.1, 10.0, 0.05, offset=0.01)
    assert abs(b["balance_residual"]) < 1e-9
