"""Sensitivity of line + inverse-bremsstrahlung absorption of a finite-bandwidth 1080 nm laser (SI; wavelengths in nm where stated).
Transmission through a uniform slab:  T(lam) = exp(-alpha_IB L - sum_j S_j L phi_j(lam)),  absorbed = 1 - integral(rho_laser T dlam).
phi_j: normalised Voigt (Doppler Gaussian + Lorentzian Stark) in nm^-1.  rho_laser: normalised Gaussian of given FWHM and centre.
Stark width and f-values are INPUT PARAMETERS here, not predictions."""
import numpy as np
from scipy.special import voigt_profile
from radiation.line_absorption import (lower_level_density_ion, line_integral_alpha_lambda, doppler_fwhm_nm,
                                        f_from_a, load_lines)
from radiation.inverse_bremsstrahlung import alpha_ib_equilibrium
C = 299792458.0
FWHM2SIG = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

# Hibbert and Hansen (1994) transition probabilities (length, velocity) [s^-1], verified via f = g_k A lam^2 / (6.6702e15 g_l) reproduction
VERIFIED_A = {  # NIST Ritz air wavelength [nm]: (A_length, A_velocity, g_lower, g_upper)
    1072.05413: (3.746e4, 4.762e4, 6, 4),
    1081.28968: (3.930e6, 3.301e6, 6, 4),
    1082.94584: (5.771e5, 4.559e5, 4, 2),
}

def f_abs(lam_nm, a, gl, gk):
    return f_from_a(a, gl, gk, lam_nm * 1e-9)

def voigt_nm(x_nm, fwhm_doppler, fwhm_lorentz):
    sig = fwhm_doppler * FWHM2SIG
    gam = fwhm_lorentz / 2.0
    return voigt_profile(x_nm, sig, gam)

def absorbed_fraction(t, p, length, lines, laser_centre_nm, laser_fwhm_nm, stark_fwhm_nm,
                      include_ib=True, grid_half_width=25.0, step=0.002, laser_shape="gaussian"):
    """lines: list of dicts {lam_nm, e_lower_cm, j_lower, f}. Returns (total, ib_only, lines_only)."""
    lam = np.arange(1080.0 - grid_half_width, 1080.0 + grid_half_width, step)
    if laser_shape == "gaussian":
        rho = np.exp(-0.5 * ((lam - laser_centre_nm) / (laser_fwhm_nm * FWHM2SIG)) ** 2)
    else:                                              # 'lorentzian': heavy-tailed bounding scenario (ASSUMED)
        rho = 1.0 / (1.0 + ((lam - laser_centre_nm) / (laser_fwhm_nm / 2.0)) ** 2)
    rho /= rho.sum() * step
    tau = np.zeros_like(lam)
    for ln in lines:
        n_l = ln["n_lower"] if "n_lower" in ln else lower_level_density_ion(t, p, ln["e_lower_cm"], ln["j_lower"])
        s = line_integral_alpha_lambda(n_l, ln["lam_nm"] * 1e-9, ln["f"])        # m^-1 nm
        tau += s * length * voigt_nm(lam - ln["lam_nm"], doppler_fwhm_nm(t, ln["lam_nm"]), ln.get("stark_nm", stark_fwhm_nm))
    tau_ib = alpha_ib_equilibrium(t, p, 1.08e-6) * length if include_ib else 0.0
    t_tot = np.sum(rho * np.exp(-tau - tau_ib)) * step
    t_ib = np.exp(-tau_ib)
    t_ln = np.sum(rho * np.exp(-tau)) * step
    return 1 - t_tot, 1 - t_ib, 1 - t_ln

def verified_lines(use="mean"):
    rows = {round(float(r["ritz_air_nm"]), 5): r for r in load_lines() if r["ion"] == "Ar II"}
    out = []
    for lam, (al, av, gl, gk) in VERIFIED_A.items():
        r = rows[round(lam, 5)]
        a = {"length": al, "velocity": av, "mean": 0.5 * (al + av)}[use]
        out.append({"lam_nm": lam, "e_lower_cm": float(r["E_i_cm1"]), "j_lower": r["lower_J"], "f": f_abs(lam, a, gl, gk)})
    return out


def neutral_ar_lines(t, p, f_values, stark_nm):
    """Ar I 1073.387 and 1075.916 nm (3d-5f). Lower-level density: n0 (g/U0) exp(-E/kT) with U0 = 1 (ground state only; ASSUMED).
    f_values: (f_1073, f_1076) are INPUTS (Coulomb-approximation bounds times an assumed line share)."""
    from thermodynamics.argon_saha import number_densities, K_B
    from radiation.line_absorption import HC_CM
    n0 = number_densities(t, p)["n0"]
    out = []
    for lam, e, j, f in ((1073.387, 114821.944, 3, f_values[0]), (1075.916, 113425.969, 2, f_values[1])):
        out.append({"lam_nm": lam, "n_lower": n0 * (2 * j + 1) * np.exp(-e * HC_CM / (K_B * t)), "f": f, "stark_nm": stark_nm})
    return out
