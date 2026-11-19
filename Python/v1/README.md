# Laser-sustained plasma (LSP) thermal propulsion model: literature-grounded foundation (1080 nm, argon)

Status: partial. Phase 1 (literature extraction) is complete for the Tsukadaira et al. (2026) benchmark and the atomic data are researched to the limit of available evidence. Phase 2 is partly implemented. The plasma energy balance, plasma-to-gas heat transfer, quasi-1D flow with heat addition, coupling, CAD, ANSYS runs and the system layer are NOT implemented. See docs/FINAL_DEVELOPMENT_REPORT.md.

## Run
    python3 run_tests.py        # 41 tests, needs numpy and scipy only (pytest not required)

## Layout
- src/laser_delivery: Gaussian beam, aperture capture, power budget
- src/nozzle: published sonic-flow estimation method (benchmark reproduction, not a thrust measurement)
- src/thermodynamics: equilibrium Ar/Ar+/e- (Saha), enthalpy, plasma checks
- src/radiation: inverse bremsstrahlung, line screening, finite-bandwidth sensitivity, Coulomb approximation
- data/processed: transcribed NIST Ar line table
- docs: registers (SOURCE, PARAMETER, ATOMIC_DATA), equations, reports, ANSYS and CAD plans
