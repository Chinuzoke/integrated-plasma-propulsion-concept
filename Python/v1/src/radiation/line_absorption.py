"""Atomic-line absorption SCREENING for the 1070-1090 nm band (NIST ASD line table, Ar I / Ar II). SI units unless named.
No transition probability (A_ki) or oscillator strength is available for any listed line, so absolute absorption cannot be computed.
This module provides (i) lower-level populations from Boltzmann equilibrium on the Saha Ar+ density, (ii) the line-integrated absorption
S = integral(alpha d lambda) as a function of the unknown oscillator strength f, and (iii) an upper bound using f <= 1 (a bound, not an estimate).
Assumptions: NIST ion level energies are measured from the ion ground level; LTE (Te = T_heavy); ion partition function U = 6 (ground term only);
line broadening and laser spectral width not modelled (both unresolved)."""
import csv, numpy as np
from thermodynamics.argon_saha import K_B, H, M_E, Q_E, EPS0, number_densities
C = 299792458.0
HC_CM = H * C * 100.0                      # J per cm^-1
U_ION = 6.0

def load_lines(path="data/processed/ar_lines_1070_1090_nm.csv"):
    with open(path) as f:
        return list(csv.DictReader(f))

def degeneracy(j_text):
    num, den = j_text.split("/")
    return float(num) / float(den) * 2 + 1          # 2J + 1

def air_wavelength_from_levels(e_i_cm, e_k_cm):
    """Air wavelength [nm] from two level energies [cm^-1]; Edlen (1966) refractive index (sigma in um^-1)."""
    lam_vac = 1e7 / (e_k_cm - e_i_cm)
    s2 = (1e3 / lam_vac) ** 2
    n = 1 + 1e-8 * (8342.13 + 2406030 / (130 - s2) + 15997 / (38.9 - s2))
    return lam_vac / n

def lower_level_density_ion(t, p, e_lower_cm, j_lower):
    """n_l = n_Ar+ (g_l/U) exp(-E_l/kT)  [m^-3]."""
    ni = number_densities(t, p)["ni"]
    return ni * degeneracy(j_lower) / U_ION * np.exp(-e_lower_cm * HC_CM / (K_B * t))

def integrated_cross_section(f=None, a_ki=None, g_l=None, g_k=None, lam=None):
    """integral sigma d nu [m^2 Hz].  From f: e^2 f/(4 eps0 m_e c).  From A_ki: lam^2 g_k A/(8 pi g_l) (lam in m, c-free form)."""
    if f is not None:
        return Q_E**2 * f / (4 * EPS0 * M_E * C)
    return lam**2 * (g_k / g_l) * a_ki / (8 * np.pi)

def line_integral_alpha_lambda(n_lower, lam_m, f):
    """S = integral(alpha d lambda) [m^-1 nm] for oscillator strength f, ignoring stimulated emission (thin, lower level only)."""
    return n_lower * integrated_cross_section(f=f) * lam_m**2 / C * 1e9

def doppler_fwhm_nm(t, lam_nm, mass_kg=39.948 * 1.66053906660e-27):
    return lam_nm / C * np.sqrt(8 * K_B * t * np.log(2) / mass_kg)

def f_from_a(a_ki, g_l, g_k, lam):
    """Absorption oscillator strength from A_ki: f = 2 pi eps0 m_e c lam^2 ... = (g_k/g_l) A eps0 m_e c lam^2 /(2 pi e^2)."""
    return (g_k / g_l) * a_ki * EPS0 * M_E * C * lam**2 / (2 * np.pi * Q_E**2)
