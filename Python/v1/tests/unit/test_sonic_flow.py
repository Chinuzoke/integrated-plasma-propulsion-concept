import sys, numpy as np
sys.path.insert(0, "src")
from nozzle.sonic_flow import *
A = np.pi / 4 * (0.5e-3) ** 2

def test_eq1_inverts_choked_flow():
    t0, p0 = 1500.0, 8e5
    md = choked_mdot(A, p0, t0)
    assert np.isclose(temperature_from_h(h_av(A, p0, md)), t0, rtol=1e-12)

def test_vacuum_velocity_is_sqrt_2h():
    t = 1200.0
    assert np.isclose(exit_velocity_vacuum(t), np.sqrt(2 * GAMMA / (GAMMA - 1) * R_AR * t))

def test_cold_velocity_matches_paper_560():
    assert abs(exit_velocity_vacuum(300.0) - 560.0) < 3.0

def test_eq6_and_eq7_are_consistent():
    md, P, vc = 2.5e-4, 1000.0, 560.0
    ve = 1363.0
    eta = eta_thrust(md * ve, md * vc, md, P)
    assert np.isclose(isp_simplified(eta, 1.0, P, md, vc), ve / G0, rtol=1e-12)

def test_throat_area_matches_paper_0p2mm2():
    assert abs(A * 1e6 - 0.2) < 0.01

def test_plenum_section_area_matches_paper():
    assert abs(np.pi / 4 * (18**2 - 17**2) - 27.5) < 0.1
