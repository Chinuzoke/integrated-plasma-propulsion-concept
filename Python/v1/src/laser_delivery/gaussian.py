"""Gaussian beam propagation and receiver capture (SI units, 1/e^2 radius convention)."""
import numpy as np

def rayleigh_range(w0, wavelength, m2=1.0):
    """zR = pi w0^2 / (M^2 lambda) [m]."""
    return np.pi * w0**2 / (m2 * wavelength)

def beam_radius(z, w0, wavelength, m2=1.0):
    """w(z) = w0 sqrt(1 + (z/zR)^2) [m]; radius, not diameter."""
    return w0 * np.sqrt(1.0 + (z / rayleigh_range(w0, wavelength, m2))**2)

def capture_fraction_centred(a, w):
    """Fraction through a centred circular aperture of radius a: 1 - exp(-2 a^2 / w^2)."""
    return 1.0 - np.exp(-2.0 * a**2 / w**2)

def capture_fraction_offset(a, w, offset, n_r=400, n_phi=720):
    """Fraction through an aperture displaced by 'offset' from the beam axis (midpoint-rule integration)."""
    r = (np.arange(n_r) + 0.5) * a / n_r
    phi = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    R, PHI = np.meshgrid(r, phi, indexing="ij")
    rho2 = R**2 + offset**2 + 2 * R * offset * np.cos(PHI)
    irr = (2.0 / (np.pi * w**2)) * np.exp(-2.0 * rho2 / w**2)
    return float(np.sum(irr * R) * (a / n_r) * (2 * np.pi / n_phi))

def power_budget(p_tx, transmissions, w0, wavelength, m2, z, a_rx, offset=0.0):
    """Sequential budget; each loss applied once; balance residual should be ~0."""
    w = beam_radius(z, w0, wavelength, m2)
    f_ctr = capture_fraction_centred(a_rx, w)
    f_cap = f_ctr if offset == 0 else capture_fraction_offset(a_rx, w, offset)
    p = p_tx * f_cap
    out = {"p_tx": p_tx, "w_rx": w, "capture": f_cap,
           "pointing_loss_fraction": 1.0 - f_cap / f_ctr if f_ctr > 0 else 0.0,
           "p_captured": p, "p_spill": p_tx - p, "stages": {}}
    for name, t in transmissions.items():
        out["stages"][name] = p * (1 - t)
        p *= t
    out["p_delivered"] = p
    out["overall_efficiency"] = p / p_tx
    out["balance_residual"] = p_tx - (p + out["p_spill"] + sum(out["stages"].values()))
    return out
