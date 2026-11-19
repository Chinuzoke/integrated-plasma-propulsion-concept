"""Electron-ion inverse-bremsstrahlung (free-free) absorption, Kramers classical form, SI units.
alpha = (4/3) * (e^2/(4 pi eps0))^3 / (m_e h c) * sqrt(2 pi/(3 k m_e)) * Z^2 n_e n_i / (sqrt(T) nu^3) * (1 - exp(-h nu/kT)) * g_ff
(cgs original: Rybicki and Lightman, Radiative Processes in Astrophysics, Eq. 5.19b; textbook, source rank 4).
Limits: gaunt factor g_ff = 1 (order-unity uncertainty, NOT tuned); electron-neutral absorption and atomic-line absorption NOT included;
Te = T_heavy; no plasma-frequency cutoff correction (omega >> omega_p here, see tests)."""
import numpy as np
from thermodynamics.argon_saha import K_B, H, M_E, Q_E, EPS0, number_densities
C = 299792458.0

def kramers_coefficient():
    """Constant A in alpha = A Z^2 n_e n_i T^-1/2 nu^-3 (1 - e^{-h nu/kT}) g [SI: m^5 K^0.5 s^-3... see test]."""
    e2 = Q_E**2 / (4 * np.pi * EPS0)
    return (4.0 / 3.0) * e2**3 / (M_E * H * C) * np.sqrt(2 * np.pi / (3 * K_B * M_E))

def alpha_ib(t, ne, ni, wavelength, z=1.0, gaunt=1.0):
    """Absorption coefficient [1/m] for intensity I(x) = I0 exp(-alpha x)."""
    nu = C / wavelength
    return kramers_coefficient() * z**2 * ne * ni / (np.sqrt(t) * nu**3) * (1 - np.exp(-H * nu / (K_B * t))) * gaunt

def alpha_ib_equilibrium(t, p, wavelength):
    d = number_densities(t, p)
    return alpha_ib(t, d["ne"], d["ni"], wavelength)

def fractional_absorption(alpha, length):
    """Uniform slab of length [m]: 1 - exp(-alpha L)."""
    return 1.0 - np.exp(-alpha * length)
