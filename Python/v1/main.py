"""Prints the key results of every implemented module.   Run from the project root:   python3 main.py
Library files under src/ only define functions and print nothing when run directly; this script calls them."""
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
import numpy as np
from laser_delivery.gaussian import power_budget, rayleigh_range
from nozzle.sonic_flow import h_av, temperature_from_h, exit_velocity_vacuum, isp, choked_mdot, G0
from thermodynamics.argon_saha import number_densities, enthalpy, cp_equilibrium
from radiation.inverse_bremsstrahlung import alpha_ib_equilibrium, fractional_absorption
from radiation.line_sensitivity import absorbed_fraction, verified_lines, neutral_ar_lines

os.makedirs("outputs/tables", exist_ok=True)
print("=== 1. Laser delivery (illustrative inputs, not a published case) ===")
b = power_budget(1000, {"lens": 0.99, "window": 0.995}, 0.005, 1.08e-6, 1.1, 10.0, 0.05, offset=0.01)
print(f"Rayleigh range {rayleigh_range(0.005,1.08e-6,1.1):.1f} m | beam radius at receiver {b['w_rx']*1e3:.2f} mm | capture {b['capture']:.4f} | "
      f"delivered {b['p_delivered']:.1f} W | balance residual {b['balance_residual']:.1e} W")

print("\n=== 2. Sonic-flow benchmark of Tsukadaira et al. (2026): estimate, not a thrust measurement ===")
A = np.pi / 4 * (0.5e-3) ** 2
print("cold choked mass flow at 0.5 MPa, 300 K: %.4f g/s" % (choked_mdot(A, 5e5, 300) * 1e3))
T = temperature_from_h(h_av(A, 1.0e6, 0.25e-3))
print(f"p = 1.0 MPa, mdot = 0.25 g/s -> mean plenum temperature {T:.0f} K, Isp {isp(exit_velocity_vacuum(T)):.0f} s (published max: 139 s)")

print("\n=== 3. Equilibrium argon at 1.0 MPa ===")
print("T [K]    ionisation   n_e [m^-3]   h [MJ/kg]   cp [J/kg K]")
rows = []
for t in (1786, 5000, 10000, 15000, 18000, 20000):
    d = number_densities(t, 1e6); rows.append([t, d["x"], d["ne"], enthalpy(t, 1e6) / 1e6, cp_equilibrium(t, 1e6)])
    print(f"{t:6d}   {d['x']:.3e}   {d['ne']:.3e}   {rows[-1][3]:8.3f}   {rows[-1][4]:9.0f}")
with open("outputs/tables/argon_equilibrium_1MPa.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["T_K", "ionisation_fraction", "ne_m3", "h_MJ_per_kg", "cp_J_per_kgK"]); w.writerows(rows)

print("\n=== 4. Absorption screening, 1 mm uniform plasma, 1.0 MPa, 1080 nm, laser 3 nm wide (CALCULATED, NOT validated) ===")
print("T [K]   IB only   + Ar II lines   + Ar I lines (share 0.1, Stark 10 nm)")
for t in (16000, 18000, 20000):
    ib = fractional_absorption(alpha_ib_equilibrium(t, 1e6, 1.08e-6), 1e-3)
    ar2 = absorbed_fraction(t, 1e6, 1e-3, verified_lines(), 1080.0, 3.0, 0.15)[0]
    ar1 = absorbed_fraction(t, 1e6, 1e-3, verified_lines() + neutral_ar_lines(t, 1e6, (0.0173, 0.0172), 10.0), 1080.0, 3.0, 0.15)[0]
    print(f"{t:6d}   {100*ib:6.2f}%   {100*ar2:8.2f}%   {100*ar1:8.2f}%")
print("\nSee docs/FINAL_DEVELOPMENT_REPORT.md for what is and is not established. Tables written to outputs/tables/.")
