"""
references.py -- LSP source and parameter register.

Status of each entry is stated honestly:
  VERIFIED_THIS_SESSION  existence / bibliographic data confirmed by a lookup during development
  RECALLED               cited from the author's recall; confirm before relying on it
  LEGACY_CODE            from the earlier LSP implementation (no external source asserted there)
  USER_SUPPLIED          provided by the user, not re-fetched
Kind: experiment | theory/numerical | textbook | data table | software/standard

NO plasma operating condition (pressure, flow rate, dimensions) is taken from these references into the
model: the model's operating inputs are ORIGINAL_MODEL_DEFAULT (see config.py).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, asdict


@dataclass(frozen=True)
class LRef:
    topic: str
    citation: str
    kind: str
    status: str
    url_or_doi: str
    used_for: str
    limits: str


REGISTER = [
    LRef("1-D LSP theory (origin of 1-D models)",
         "Yu. P. Raizer, 'Subsonic propagation of a light spark and threshold conditions for the maintenance of plasma by radiation', "
         "Sov. Phys. JETP 31, 1148-1154 (1970)", "theory", "VERIFIED_THIS_SESSION (appears in the reference list of Jeng & Keefer 1987)",
         "not retrieved", "historical basis of the 1-D formulation", "Original not read in this session"),
    LRef("Limits of 1-D models; 2-D LTE argon/hydrogen LSP",
         "S.-M. Jeng and D. Keefer, 'A Theoretical Investigation of Laser-Sustained Plasma Thruster', AIAA-87-0383, 25th Aerospace "
         "Sciences Meeting, Reno, 1987", "theory/numerical", "VERIFIED_THIS_SESSION (document read)",
         "https://apps.dtic.mil/sti/pdfs/ADA193765.pdf",
         "states that 1-D models are inherently limited because the LSP depends strongly on beam geometry and radial transport; "
         "LTE assumed at their pressures; optically thin + thick radiation split",
         "Hydrogen propellant, 10.6 um, focused beam: not directly comparable to argon at 1 um"),
    LRef("Argon LSP 2-D model",
         "S.-M. Jeng and D. Keefer, 'Theoretical investigation of laser-sustained argon plasmas', J. Appl. Phys. 60(7), 2272-2279 (1986)",
         "theory/numerical", "VERIFIED_THIS_SESSION (appears in reference list of Jeng & Keefer 1987)", "not recorded",
         "benchmark class for a future 2-D upgrade", "Not read in this session; no values used"),
    LRef("Argon LSP numerical model",
         "Akarapu, Nassar, Copley, Todd, 'Numerical model of a laser-sustained argon plasma', J. Laser Appl. 21, 169-175 (2009)",
         "theory/numerical", "VERIFIED_THIS_SESSION (search listing only)", "https://pure.psu.edu/en/publications/numerical-model-of-a-laser-sustained-argon-plasma/",
         "comparison target for a future 2-D model", "Not read in this session"),
    LRef("Fiber/diode-laser argon LSP experiments (Shizuoka)",
         "Homme, Takano, Arakaki, Matsui, 'Evaluation of generating conditions and characteristics of argon laser sustained plasma using fiber laser', "
         "Applied Plasma Science 31(2) (2024); reported (search abstract): minimum laser power ~460 W at 1.0 MPa, "
         "absorption 12-40 %, temperature 11,900-15,700 K", "experiment", "VERIFIED_THIS_SESSION (abstract-level only)",
         "https://www.jstage.jst.go.jp/article/aps/31/2/31_61/_article/-char/en",
         "QUALITATIVE plausibility check only (MPa-level pressure at ~1 um; T ~1.2e4-1.6e4 K; absorption tens of %)",
         "Focused-beam, multi-dimensional geometry: a 1-D collimated model cannot be quantitatively compared"),
    LRef("Diode-laser argon LSP experiment",
         "'Generation of argon LSP using diode laser', Trans. Japan Soc. Aeronaut. Space Sci. 66(4) (T-22-55); reported (abstract): 4-kW diode laser, "
         "1.55-2.00 MPa", "experiment", "VERIFIED_THIS_SESSION (abstract-level only)",
         "https://www.jstage.jst.go.jp/article/tjsass/66/4/66_T-22-55/_article/-char/en",
         "qualitative: ignition thresholds at MPa pressures for ~1 um class lasers", "Authors/details not recorded"),
    LRef("Argon LSP thruster laboratory model",
         "Duplay, 'Argon Laser-Plasma Thruster', TU Delft repository (2024) (3 kW, 1070 nm fibre laser per the listing)",
         "experiment/design", "VERIFIED_THIS_SESSION (listing only)", "https://repository.tudelft.nl/record/uuid:3a853f2e-2d8c-496b-a016-9e9855e8069c",
         "context only", "Not read"),
    LRef("Free-free absorption and emission (IB, brems)",
         "G. B. Rybicki and A. P. Lightman, Radiative Processes in Astrophysics, Wiley (1979), eqs. 5.15b, 5.19b",
         "textbook", "RECALLED", "ISBN not recorded",
         "alpha_ff and eps_ff (SI forms verified numerically against the legacy cgs constants in validation.py)",
         "Kramers/Gaunt approximation; e-ion only"),
    LRef("Argon thermodynamic/transport properties (LTE)",
         "A. B. Murphy and C. J. Arundell, 'Transport coefficients of argon, nitrogen, oxygen, argon-nitrogen, and argon-oxygen plasmas', "
         "Plasma Chem. Plasma Process. 14, 451-490 (1994)", "data table", "VERIFIED_THIS_SESSION (existence)",
         "https://doi.org/10.1007/BF01570207", "intended source for a sourced conductivity table (TabulatedConductivity)",
         "NO values were extracted or used here"),
    LRef("Thermal plasma fundamentals",
         "M. I. Boulos, P. Fauchais, E. Pfender, Thermal Plasmas: Fundamentals and Applications, Vol. 1, Plenum (1994)",
         "textbook", "RECALLED", "not recorded", "background for LTE, transport, radiation", "no values used"),
    LRef("Viscosity collision integral fit",
         "P. D. Neufeld, A. R. Janzen, R. A. Aziz, J. Chem. Phys. 57, 1100 (1972) (Omega22* fit)", "theory fit", "RECALLED",
         "not recorded", "Ar-Ar viscosity -> heavy-particle conductivity", "valid T* 0.3-100; extrapolated above ~9,300 K"),
    LRef("Lennard-Jones argon parameters", "Hirschfelder, Curtiss, Bird, Molecular Theory of Gases and Liquids (1954): sigma 3.542 A, eps/k 93.3 K",
         "textbook table", "RECALLED", "not recorded", "Ar-Ar viscosity", "mu(300 K) check ~2.27e-5 Pa s passes (validation.py)"),
    LRef("Argon ionisation energies, mass", "NIST Atomic Spectra Database values 15.7596 eV, 27.6297 eV; M = 39.948 g/mol", "data table",
         "LEGACY_CODE (values from legacy constants.py; agree with recalled NIST values)", "https://physics.nist.gov/asd", "Saha equilibrium",
         "not re-fetched"),
    LRef("Beam propagation / M^2", "see delivery register (Siegman 1986; ISO 11146-1)", "textbook/standard", "RECALLED", "n/a",
         "collimation consistency check only", "M^2 = 11.8 is UNVERIFIED in the delivery register"),
    LRef("Delivery-subsystem inputs", "laser power, wavelength, M^2 come from power_budget.compute_budget()", "calculated", "USER_SUPPLIED",
         "n/a", "LaserInput", "peak-power reference scenario; see delivery register"),
    LRef("e-Ar momentum-transfer cross-section", "none", "data", "UNRESOLVED", "n/a", "electron-neutral IB; electron conductivity",
         "required for ignition-threshold accuracy"),
    LRef("Net emission coefficient for argon", "none", "data table", "UNRESOLVED", "n/a", "line + recombination radiation", "brems-only otherwise"),
    LRef("Reactive + electron thermal conductivity", "none", "data table", "UNRESOLVED", "n/a", "conduction closure", "conduction underestimated otherwise"),
]


def write_csv(path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(asdict(REGISTER[0]).keys()))
        wr.writeheader()
        for r in REGISTER:
            wr.writerow(asdict(r))
