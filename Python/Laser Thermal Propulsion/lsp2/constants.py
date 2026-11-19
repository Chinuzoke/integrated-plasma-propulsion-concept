"""
constants.py -- fundamental constants (CODATA 2018 exact/recommended values, SI) and argon data.

Provenance
----------
* Universal constants: CODATA 2018 (exact SI-defining values; m_e, eps0 recommended values), as in the
  legacy constants.py (retained unchanged in value).
* Argon molar mass 39.948 g/mol (standard atomic weight) and ionisation energies 15.7596 eV, 27.6297 eV
  (first and second): taken from the legacy file; they agree with NIST Atomic Spectra Database values as
  recalled by the author, but were NOT re-fetched in this session -> verify against NIST ASD.
* Statistical weights (1, 6, 9): ground-TERM weights (Ar I 1S0 -> 1; Ar II 2P -> 2+4 = 6; Ar III 3P -> 9).
  Excited electronic levels are neglected in the partition functions (error grows with T; a few % in
  the Saha ratio near 2e4 K). Ionisation-potential lowering (Debye) is neglected (~0.05 eV at
  n_e ~ 1e23 m^-3, ~3 % in the Saha ratio at 2e4 K).
"""
import math

k_B = 1.380649e-23          # J/K
e_charge = 1.602176634e-19  # C
h_planck = 6.62607015e-34   # J s
c_light = 2.99792458e8      # m/s
m_e = 9.1093837015e-31      # kg
m_u = 1.66053906660e-27     # kg
N_A = 6.02214076e23         # 1/mol
eps0 = 8.8541878128e-12     # F/m
sigma_SB = 5.670374419e-8   # W/(m^2 K^4)
R_gas = 8.314462618         # J/(mol K)

AR_MOLAR_MASS = 39.948e-3                   # kg/mol
AR_ATOM_MASS = AR_MOLAR_MASS / N_A          # kg per atom (electron-mass difference of Ar+ neglected)
AR_E1 = 15.7596 * e_charge                  # J, Ar -> Ar+
AR_E2 = 27.6297 * e_charge                  # J, Ar+ -> Ar2+
AR_G = (1.0, 6.0, 9.0)                      # ground-term statistical weights (Ar, Ar+, Ar2+)

# Lennard-Jones (12-6) parameters for Ar-Ar, Hirschfelder-Curtiss-Bird tabulation (recalled; verify).
AR_LJ_SIGMA_ANGSTROM = 3.542
AR_LJ_EPS_OVER_K = 93.3                     # K
