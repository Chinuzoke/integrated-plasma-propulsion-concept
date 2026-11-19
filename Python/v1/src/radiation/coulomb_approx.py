"""Numerical Coulomb approximation (quantum-defect, single active electron) for one-electron radial integrals and oscillator strengths.
Atomic units. Radial equation for u = r R:  u'' = [ l(l+1)/r^2 - 2/r + 1/n*^2 ] u, with n* = sqrt(Ry/binding) (non-integer for non-hydrogenic states).
The decaying solution is integrated INWARD by Numerov from large r through the inner classical turning point and cut at the first local minimum of |u| (set to zero inside).
Accuracy is benchmarked on hydrogen only (tests); for argon it is an order-of-magnitude estimate (tens of percent)."""
import numpy as np

def radial_function(n_star, l, h=0.002):
    r_in = n_star**2 * (1 - np.sqrt(1 - l * (l + 1) / n_star**2))
    r_max = 4.0 * n_star**2 + 60.0
    r_low = max(0.25 * r_in, 5 * h)
    n = int((r_max - r_low) / h) + 1
    r = r_max - h * np.arange(n)                     # decreasing grid
    k2 = lambda x: l * (l + 1) / x**2 - 2.0 / x + 1.0 / n_star**2    # u'' = k2 u
    u = np.zeros(n)
    u[0] = 1e-30
    u[1] = u[0] * np.exp(h / n_star)
    for i in range(1, n - 1):
        f0, f1, f2 = k2(r[i - 1]), k2(r[i]), k2(r[i + 1])
        u[i + 1] = (2 * u[i] * (1 + 5 * h * h * f1 / 12) - u[i - 1] * (1 - h * h * f0 / 12)) / (1 - h * h * f2 / 12)
        if abs(u[i + 1]) > 1e150:                     # rescale to avoid overflow during inward growth
            u[: i + 2] *= 1e-150
    # cut inside the inner turning point at the first local minimum of |u| (inward growth of the unphysical solution), else keep all
    a = np.abs(u)
    idx_in = np.where(r <= r_in)[0]
    cut = n
    if idx_in.size:
        i0 = idx_in[0]
        for i in range(i0, n - 1):
            if a[i + 1] > a[i]:
                cut = i + 1; break
    r = r[:cut][::-1]; u = u[:cut][::-1]
    u /= np.sqrt(np.sum(u**2) * h)
    return r, u

def radial_integral(n1, l1, n2, l2, h=0.002):
    """<n1 l1| r |n2 l2> in a0 (n1, n2 effective quantum numbers)."""
    r1, u1 = radial_function(n1, l1, h)
    r2, u2 = radial_function(n2, l2, h)
    r = np.arange(min(r1[0], r2[0]), min(r1[-1], r2[-1]), h)
    f1 = np.interp(r, r1, u1, left=0.0); f2 = np.interp(r, r2, u2, left=0.0)
    return np.sum(f1 * r * f2) * h

def oscillator_strength(n_lo, l_lo, n_up, l_up, delta_e_hartree, h=0.002):
    """One-electron absorption oscillator strength f = (2/3) dE (l_>/(2 l_lo + 1)) R^2   (dE in hartree, R in a0)."""
    lmax = max(l_lo, l_up)
    R = radial_integral(n_lo, l_lo, n_up, l_up, h)
    return 2.0 / 3.0 * delta_e_hartree * lmax / (2 * l_lo + 1) * R**2
