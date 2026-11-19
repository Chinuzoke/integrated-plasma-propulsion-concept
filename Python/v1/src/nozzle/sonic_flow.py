"""Sonic-flow estimation method of Tsukadaira et al. (2026), Eqs. (1)-(7).
Benchmark reproduction of the published ESTIMATION method, not a thrust measurement.
SI units. Ideal monatomic gas, constant gamma, frozen vacuum expansion (as in the paper).
Provisional deviation: the paper obtains T_pl from h_av with NASA CEA; here h = cp*T (valid while argon is un-ionised)."""
import numpy as np

R_AR = 8.314462618 / 0.039948      # J/(kg K): CODATA R / IUPAC molar mass of Ar (reference data, not from the paper)
GAMMA = 5.0 / 3.0                  # paper uses 1.667
G0 = 9.80665                       # m/s^2

def h_av(a_throat, p_pl, mdot, gamma=GAMMA):
    """Eq.(1): mean specific enthalpy [J/kg] from choked-flow mass flux. Assumes discharge coefficient = 1."""
    return (a_throat * p_pl / mdot) ** 2 * gamma**2 / (gamma - 1) * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))

def temperature_from_h(h, gamma=GAMMA, r=R_AR):
    """h = cp T with cp = gamma R/(gamma-1)  [K]  (replaces CEA; provisional)."""
    return h * (gamma - 1) / (gamma * r)

def choked_mdot(a_throat, p0, t0, gamma=GAMMA, r=R_AR):
    """Ideal choked mass flow [kg/s]."""
    return a_throat * p0 * np.sqrt(gamma / (r * t0)) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))

def exit_velocity_vacuum(t_pl, gamma=GAMMA, r=R_AR):
    """Eq.(3) with p_e = 0: v_e = sqrt(2 gamma/(gamma-1) R T_pl) = sqrt(2 h)."""
    return np.sqrt(2 * gamma / (gamma - 1) * r * t_pl)

def thrust(mdot, ve):            return mdot * ve            # Eq.(4)
def isp(ve):                     return ve / G0               # Eq.(5)

def eta_thrust(f, f_cold, mdot, p_laser):
    """Eq.(6): (F^2 - F_cold^2) / (2 mdot P_laser)."""
    return (f**2 - f_cold**2) / (2 * mdot * p_laser)

def isp_simplified(eta_abs, eta_nozzle, p_laser, mdot, v_cold):
    """Eq.(7)."""
    return np.sqrt(2 * eta_abs * eta_nozzle * p_laser / mdot + v_cold**2) / G0
