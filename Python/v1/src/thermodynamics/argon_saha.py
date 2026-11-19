"""Equilibrium Ar / Ar+ / e- composition and thermodynamics (single ionisation, LTE, Te = Th). SI units.
Limits: ionisation-potential lowering, excited-state partition-function growth, Ar2+ and molecular ions neglected.
Valid as a baseline roughly below ~20,000 K; Ar2+ onset (27.63 eV) is not modelled."""
import numpy as np

K_B = 1.380649e-23; H = 6.62607015e-34; M_E = 9.1093837015e-31
Q_E = 1.602176634e-19; EPS0 = 8.8541878128e-12; AMU = 1.66053906660e-27
M_AR = 39.948 * AMU                  # kg (IUPAC standard atomic weight; reference data)
E_ION = 15.76 * Q_E                  # J (S02 Ch.4.1 citing NIST; parameter register P06)
G_ION, G_NEUTRAL = 6.0, 1.0          # ground-term degeneracies (Ar+ 2P total, Ar 1S0); assumption

def saha_constant(t):
    """S(T) = n_e n_i / n_0 [m^-3] = 2 (U_i/U_0) (2 pi m_e k T / h^2)^(3/2) exp(-E_i/kT)."""
    return 2.0 * G_ION / G_NEUTRAL * (2 * np.pi * M_E * K_B * t / H**2) ** 1.5 * np.exp(-E_ION / (K_B * t))

def ionisation_fraction(t, p):
    """x = n_e/(n_0+n_i) from  x^2/(1-x^2) = S k T / p   =>  x = sqrt(K/(1+K)),  K = S k T / p."""
    k = saha_constant(t) * K_B * t / p
    return np.sqrt(k / (1.0 + k))

def number_densities(t, p):
    x = ionisation_fraction(t, p)
    n_atoms = p / (K_B * t * (1.0 + x))      # heavy nuclei per m^3 (n0 + ni)
    return {"x": x, "n_atoms": n_atoms, "n0": n_atoms * (1 - x), "ni": n_atoms * x, "ne": n_atoms * x}

def density(t, p):
    """Mass density [kg/m^3]; electron mass neglected (m_e/m_Ar = 1.4e-5)."""
    return number_densities(t, p)["n_atoms"] * M_AR

def enthalpy(t, p):
    """h = [ (5/2)(1+x) k T + x E_i ] / m_Ar   [J/kg]  (translational enthalpy of all species + ionisation energy)."""
    x = ionisation_fraction(t, p)
    return (2.5 * (1.0 + x) * K_B * t + x * E_ION) / M_AR

def cp_equilibrium(t, p, dt=1.0):
    """Equilibrium specific heat [J/(kg K)] by central difference of h at fixed p."""
    return (enthalpy(t + dt, p) - enthalpy(t - dt, p)) / (2 * dt)

def debye_length(t, ne):
    return np.sqrt(EPS0 * K_B * t / (ne * Q_E**2))

def plasma_frequency(ne):
    return np.sqrt(ne * Q_E**2 / (EPS0 * M_E))

def mcwhirter_min_ne(t, delta_e_ev):
    """McWhirter necessary (not sufficient) LTE criterion: n_e >= 1.6e12 T^(1/2) dE^3 [cm^-3] -> returned in m^-3."""
    return 1.6e12 * np.sqrt(t) * delta_e_ev**3 * 1e6
