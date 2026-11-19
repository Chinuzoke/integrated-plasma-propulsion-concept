"""
config.py -- LSP configuration with provenance for every input.

The LSP subsystem owns plasma/chamber inputs. It does NOT own laser power, wavelength or beam
quality: these arrive through `LaserInput`, built either by `lsp.integration` from the laser-delivery
budget (the integrated path) or explicitly by the caller (standalone tests). There is no default laser power.

Provenance classes
    ORIGINAL_MODEL_DEFAULT  value inherited from the earlier LSP implementation. It is NOT an experimental or
                            literature operating condition, has no source, and is illustrative only.
    MODELLING_ASSUMPTION    a choice made in this reconstruction (adjustable).
    FROM_LASER_DELIVERY     supplied by the laser-delivery subsystem.
    NUMERICAL_SETTING       solver/grid parameter.
    UNRESOLVED              needed for a quantitative prediction but not available (no value assumed).
"""
from __future__ import annotations

from dataclasses import dataclass, replace, field
from enum import Enum


class Prov(Enum):
    ORIGINAL_MODEL_DEFAULT = "inherited from earlier LSP code; unsourced; illustrative"
    MODELLING_ASSUMPTION = "modelling assumption (adjustable)"
    FROM_LASER_DELIVERY = "supplied by the laser-delivery subsystem"
    NUMERICAL_SETTING = "numerical setting"
    UNRESOLVED = "unresolved (no value assumed)"


@dataclass(frozen=True)
class LaserInput:
    """Laser state at the LSP chamber entrance. Built from the delivery budget or given explicitly."""
    power_w: float               # optical power ENTERING the chamber [W]
    wavelength_m: float          # [m]
    m2: float | None = None      # beam propagation ratio (used only for the collimation consistency check)
    source: str = "explicit input"
    power_basis: str = "unspecified"

    def __post_init__(self):
        import math
        if not (isinstance(self.power_w, (int, float)) and math.isfinite(self.power_w) and self.power_w > 0):
            raise ValueError(f"laser power must be finite and > 0, got {self.power_w!r}")
        if not (math.isfinite(self.wavelength_m) and self.wavelength_m > 0):
            raise ValueError(f"wavelength must be finite and > 0, got {self.wavelength_m!r}")
        if self.m2 is not None and not (math.isfinite(self.m2) and self.m2 >= 1.0):
            raise ValueError(f"M^2 must be >= 1, got {self.m2!r}")


# Legacy operating inputs (ORIGINAL_MODEL_DEFAULT), kept for reproducibility of the earlier case.
LEGACY_OPERATING = dict(pressure_pa=1.0e5, mdot_kg_s=5.0e-5, plasma_radius_m=2.0e-3, length_m=0.02,
                        T_inlet_K=500.0, T_wall_K=500.0)


@dataclass(frozen=True)
class LSPConfig:
    laser: LaserInput
    # ---- operating inputs: ORIGINAL_MODEL_DEFAULT (not from experiments or literature) ----
    pressure_pa: float = LEGACY_OPERATING["pressure_pa"]
    mdot_kg_s: float = LEGACY_OPERATING["mdot_kg_s"]
    plasma_radius_m: float = LEGACY_OPERATING["plasma_radius_m"]   # = channel radius = beam radius (1-D assumption)
    length_m: float = LEGACY_OPERATING["length_m"]
    T_inlet_K: float = LEGACY_OPERATING["T_inlet_K"]
    T_wall_K: float = LEGACY_OPERATING["T_wall_K"]
    # ---- modelling choices ----
    physics: str = "reconstructed"          # 'reconstructed' | 'legacy_equivalent'
    wall_closure: str = "kirchhoff"         # 'kirchhoff' | 'adiabatic' | 'legacy_newton'
    gaunt_ff: float = 1.0                   # IB Gaunt factor (ASSUMED, as legacy)
    gaunt_brems: float = 1.2                # brems Gaunt factor (Rybicki & Lightman)
    k_multiplier: float = 1.0               # bracket for the unresolved reactive/electron conductivity
    # legacy (tuned) wall parameters, used only with wall_closure == 'legacy_newton' or physics == 'legacy_equivalent'
    legacy_wall_h: float = 100.0
    legacy_wall_emissivity: float = 3.0e-4
    # ---- numerics ----
    n_cells: int = 200

    def __post_init__(self):
        import math
        for name in ("pressure_pa", "mdot_kg_s", "plasma_radius_m", "length_m", "T_inlet_K", "T_wall_K",
                     "gaunt_ff", "gaunt_brems", "k_multiplier"):
            v = getattr(self, name)
            if not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0):
                raise ValueError(f"{name} must be finite and > 0, got {v!r}")
        if self.n_cells < 10:
            raise ValueError("n_cells must be >= 10")
        if self.physics not in ("reconstructed", "legacy_equivalent"):
            raise ValueError("physics must be 'reconstructed' or 'legacy_equivalent'")
        if self.wall_closure not in ("kirchhoff", "adiabatic", "legacy_newton"):
            raise ValueError("wall_closure must be 'kirchhoff', 'adiabatic' or 'legacy_newton'")

    @property
    def area_m2(self) -> float:
        import math
        return math.pi * self.plasma_radius_m**2

    @property
    def I0_w_m2(self) -> float:
        """Top-hat entrance intensity: delivered power over the plasma-column cross-section."""
        return self.laser.power_w / self.area_m2

    @property
    def mass_flux(self) -> float:
        return self.mdot_kg_s / self.area_m2

    def with_(self, **kw) -> "LSPConfig":
        return replace(self, **kw)

    def beam_consistency(self) -> dict:
        """Check the 1-D collimated-top-hat assumption against the delivery beam quality, if supplied.
        Rayleigh range z_R = pi R^2 / (M^2 lambda) must greatly exceed the domain length."""
        import math
        out = {"checked": False}
        if self.laser.m2 is not None:
            zR = math.pi * self.plasma_radius_m**2 / (self.laser.m2 * self.laser.wavelength_m)
            out = {"checked": True, "rayleigh_range_m": zR, "length_m": self.length_m,
                   "ratio": zR / self.length_m, "collimated_ok": zR >= 10.0 * self.length_m}
        return out


def parameter_register(cfg: LSPConfig) -> list:
    """[(name, value, unit, Prov, note)] -- one row per input, for the report and the CSV register."""
    P = Prov
    L = cfg.laser
    return [
        ("laser.power_w", L.power_w, "W", P.FROM_LASER_DELIVERY if "delivery" in L.source else P.MODELLING_ASSUMPTION,
         f"source: {L.source}; basis: {L.power_basis}"),
        ("laser.wavelength_m", L.wavelength_m, "m", P.FROM_LASER_DELIVERY if "delivery" in L.source else P.MODELLING_ASSUMPTION, L.source),
        ("laser.m2", L.m2 if L.m2 is not None else "n/a", "-", P.FROM_LASER_DELIVERY if "delivery" in L.source else P.MODELLING_ASSUMPTION,
         "UNVERIFIED in the delivery register; used only for the collimation check"),
        ("pressure_pa", cfg.pressure_pa, "Pa", P.ORIGINAL_MODEL_DEFAULT, "illustrative; NOT an experimental condition"),
        ("mdot_kg_s", cfg.mdot_kg_s, "kg/s", P.ORIGINAL_MODEL_DEFAULT, "illustrative; flow rate unspecified by any source here"),
        ("plasma_radius_m", cfg.plasma_radius_m, "m", P.ORIGINAL_MODEL_DEFAULT, "= channel = beam radius (top-hat, 1-D)"),
        ("length_m", cfg.length_m, "m", P.ORIGINAL_MODEL_DEFAULT, "chamber/absorption length unspecified by any source"),
        ("T_inlet_K", cfg.T_inlet_K, "K", P.ORIGINAL_MODEL_DEFAULT, "inlet gas temperature"),
        ("T_wall_K", cfg.T_wall_K, "K", P.ORIGINAL_MODEL_DEFAULT, "containment wall temperature"),
        ("physics", cfg.physics, "-", P.MODELLING_ASSUMPTION, "reconstructed enthalpy/Kirchhoff physics vs legacy-equivalent"),
        ("wall_closure", cfg.wall_closure, "-", P.MODELLING_ASSUMPTION, "radial-loss closure for the 1-D model"),
        ("gaunt_ff", cfg.gaunt_ff, "-", P.MODELLING_ASSUMPTION, "IB Gaunt factor"),
        ("gaunt_brems", cfg.gaunt_brems, "-", P.MODELLING_ASSUMPTION, "brems Gaunt factor (Rybicki & Lightman ~1.2)"),
        ("k_multiplier", cfg.k_multiplier, "-", P.MODELLING_ASSUMPTION, "brackets the omitted reactive+electron conductivity"),
        ("n_cells", cfg.n_cells, "-", P.NUMERICAL_SETTING, "finite-volume cells"),
        ("e-Ar momentum-transfer cross-section", "unresolved", "m^2", P.UNRESOLVED, "needed for e-neutral IB and electron conductivity"),
        ("net emission coefficient table", "unresolved", "W/m^3/sr", P.UNRESOLVED, "line + recombination radiation"),
        ("reactive + electron thermal conductivity", "unresolved", "W/m/K", P.UNRESOLVED, "sourced table needed"),
        ("laser focusing / beam profile in chamber", "unresolved", "-", P.UNRESOLVED, "1-D top-hat assumed; focusing optics not modelled"),
        ("two-temperature (T_e != T_heavy) data", "unresolved", "-", P.UNRESOLVED, "LTE assumed"),
    ]
