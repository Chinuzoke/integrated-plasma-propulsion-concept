import sys, numpy as np
sys.path.insert(0, "src")
from radiation.line_sensitivity import *

def test_f_reproduces_hibbert_hansen_tabulated_f_values():
    # table 5 (f_l, f_v): 1072.05 nm -> 0.0004, 0.0005 ; 1081.29 nm -> 0.0459, 0.0386 ; 1082.95 nm -> 0.0051, 0.0040
    expect = {1072.05413: (0.0004, 0.0005), 1081.28968: (0.0459, 0.0386), 1082.94584: (0.0051, 0.0040)}
    for lam, (al, av, gl, gk) in VERIFIED_A.items():
        fl, fv = f_abs(lam, al, gl, gk), f_abs(lam, av, gl, gk)
        assert abs(fl - expect[lam][0]) < max(0.06 * expect[lam][0], 6e-5)
        assert abs(fv - expect[lam][1]) < max(0.06 * expect[lam][1], 6e-5)

def test_voigt_profile_normalised():
    x = np.arange(-60, 60, 0.001)
    assert np.isclose(np.sum(voigt_nm(x, 0.015, 0.5)) * 0.001, 1.0, rtol=2e-2)

def test_no_lines_gives_ib_only():
    tot, ibo, lno = absorbed_fraction(18000, 1e6, 1e-3, [], 1080.0, 3.0, 0.1)
    assert np.isclose(tot, ibo) and abs(lno) < 1e-12

def test_more_oscillator_strength_absorbs_more():
    ln = verified_lines()
    ln2 = [dict(l, f=l["f"] * 5) for l in ln]
    a = absorbed_fraction(20000, 1e6, 1e-3, ln, 1080.0, 3.0, 0.15)[0]
    b = absorbed_fraction(20000, 1e6, 1e-3, ln2, 1080.0, 3.0, 0.15)[0]
    assert b > a

def test_thin_limit_matches_integrated_strength():
    # broad Stark width, tiny f: absorbed fraction by line ~ rho(lam0) * S L  (optically thin)
    ln = [dict(verified_lines()[1], f=1e-4)]
    _, _, lno = absorbed_fraction(20000, 1e6, 1e-3, ln, 1081.28968, 3.0, 0.02, include_ib=False)
    from radiation.line_absorption import lower_level_density_ion, line_integral_alpha_lambda
    s = line_integral_alpha_lambda(lower_level_density_ion(20000, 1e6, ln[0]["e_lower_cm"], ln[0]["j_lower"]), 1081.28968e-9, 1e-4) * 1e-3
    rho0 = 1 / (3.0 * FWHM2SIG * np.sqrt(2 * np.pi))
    assert np.isclose(lno, rho0 * s, rtol=0.05)

def test_neutral_lines_and_laser_shape_options_run():
    ar1 = neutral_ar_lines(20000, 1e6, (0.017, 0.017), 1.0)
    g = absorbed_fraction(20000, 1e6, 1e-3, ar1, 1080.0, 3.0, 1.0, include_ib=False)[2]
    l = absorbed_fraction(20000, 1e6, 1e-3, ar1, 1080.0, 3.0, 1.0, include_ib=False, laser_shape="lorentzian")[2]
    assert 0.0 < g < 1.0 and 0.0 < l < 1.0 and l > g
