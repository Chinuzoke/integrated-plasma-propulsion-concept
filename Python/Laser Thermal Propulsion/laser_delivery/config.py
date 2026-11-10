"""
config.py -- central reference parameters and user-adjustable inputs (register v1).

EDIT INPUTS ONLY HERE (the PARAMS table). No calculation module contains a
hard-coded physical input.

Every parameter carries a classification so that guaranteed specifications,
typical measurements, published design inputs, modelling choices, derived
values, unverified scenario values and unresolved items are never confused.
Nothing is silently upgraded: a parameter only changes class by editing this
table (and references.py, which tests keep in agreement).

    PUBLISHED_LITERATURE                 design input from a published study; see `provisional` flag)
    MANUFACTURER_SPECIFICATION           guaranteed/stated manufacturer specification
    PROVISIONAL_MANUFACTURER_MEASUREMENT typical measured manufacturer value (not guaranteed)
    ASSUMED_MODELLING_CHOICE             selected by the modeller; adjustable
    DESIGN_DERIVED                       follows from a conventional design, not a published prescription
    UNVERIFIED                           no source established; illustrative scenario input only
    UNRESOLVED                           no defensible value; value=None; requesting it raises an error

`used=False` marks parameters that are recorded for traceability but do not
enter the power calculation.

Internal units are SI throughout (m, W, rad, s). Arcseconds are converted at the boundary.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from enum import Enum

# 1 arcsec in rad (exact definition)
ARCSEC_TO_RAD = math.pi / (180.0 * 3600.0)


class Status(Enum):
    PUBLISHED_LITERATURE = "published literature"
    MANUFACTURER_SPECIFICATION = "manufacturer specification"
    PROVISIONAL_MANUFACTURER_MEASUREMENT = "typical manufacturer measurement (not guaranteed)"
    ASSUMED_MODELLING_CHOICE = "assumed modelling choice"
    DESIGN_DERIVED = "design-derived (not a manufacturer prescription)"
    UNVERIFIED = "unverified (illustrative scenario input)"
    UNRESOLVED = "unresolved (no value assumed)"


@dataclass(frozen=True)
class Param:
    value: float | str | None
    unit: str
    status: Status
    note: str = ""
    provisional: bool = False   # classified as provisional in the register
    used: bool = True           # False: recorded only, not used in the calculation


S = Status
# --------------------------------------------------------------------------
# THE PARAMETER TABLE -- edit values here.
# --------------------------------------------------------------------------
PARAMS: dict[str, Param] = {
    # ===== A. transmitter and propagation ====================================
    "transmitter_power": Param(30_000.0, "W", S.PUBLISHED_LITERATURE,
                               "PEAK-power reference scenario (30 kW Yb-doped fibre laser, MacRobbie et al. 2025). The source does not "
                               "establish CW operation or that this is the continuous output. NOT a claim of sustained 30 kW operation.",
                               provisional=True),
    "cw_average_optical_power": Param(None, "W", S.UNRESOLVED, "continuous-wave / average optical power", used=False),
    "duty_cycle": Param(None, "-", S.UNRESOLVED, "duty cycle; if supplied, an average-power line is reported"),
    "wavelength": Param(1064e-9, "m", S.PUBLISHED_LITERATURE, "representative Yb-fibre laser wavelength", provisional=True),
    "beam_waist_radius": Param(0.25, "m", S.PUBLISHED_LITERATURE, "250 mm per the cited study; 1/e^2 intensity radius",
                               provisional=True),
    "m_squared": Param(11.8, "-", S.UNVERIFIED,
                       "NO SOURCE established. Retained only as an assumed sensitivity-case input; assumed a second-moment "
                       "value, equal in x and y. Replace when a source is identified."),
    "distance": Param(1.2e6, "m", S.ASSUMED_MODELLING_CHOICE,
                      "selected 1,200 km test distance; NOT an operating distance from the cited study"),
    "receiver_radius": Param(5.0, "m", S.ASSUMED_MODELLING_CHOICE,
                             "10 m diameter; NOT a receiver specification from the cited study; hard-edged, unobstructed, centred"),
    "pointing_error": Param(0.1, "arcsec", S.UNVERIFIED,
                            "illustrative scenario only (no source); treated as a fixed offset; see pointing sweep and jitter variant"),
    "propagation_medium": Param("vacuum", "-", S.ASSUMED_MODELLING_CHOICE,
                                "initial line-of-sight simplification. Atmospheric extinction/turbulence excluded. This does NOT imply "
                                "that pointing jitter, wavefront error or other disturbances are negligible", used=False),
    "propagation_loss_fraction": Param(0.0, "-", S.ASSUMED_MODELLING_CHOICE,
                                       "vacuum: zero. Hook for future medium attenuation"),
    # ===== B. receiving mirror (representative: Excelitas/LINOS DLHS NIR 1064 nm) =====
    "optical_component_wavelength": Param(1064e-9, "m", S.MANUFACTURER_SPECIFICATION,
                                          "operating wavelength stated for the mirror, expander and window datasheets", used=False),
    "mirror_aoi_deg": Param(45.0, "deg", S.MANUFACTURER_SPECIFICATION,
                            "catalogue part is specified at 45 deg AOI; a spacecraft receiving mirror need not operate there. "
                            "Angle/polarisation dependence is not modelled", used=False),
    "mirror_reflectance": Param(0.998, "-", S.MANUFACTURER_SPECIFICATION,
                                "GUARANTEED > 99.8 %; 0.998 used as the conservative baseline"),
    "mirror_typical_reflectance": Param(0.9998, "-", S.PROVISIONAL_MANUFACTURER_MEASUREMENT,
                                        "typical measured > 99.98 %; alternative scenario only, not the guaranteed baseline", provisional=True),
    "mirror_typical_abs_scatter_ppm": Param(150.0, "ppm", S.PROVISIONAL_MANUFACTURER_MEASUREMENT,
                                            "typical measured COMBINED absorption+scatter < 150 ppm (an upper bound); not separately specified",
                                            provisional=True),
    "mirror_absorption_fraction": Param(None, "-", S.UNRESOLVED, "individual absorption not specified"),
    "mirror_scatter_fraction": Param(None, "-", S.UNRESOLVED, "individual scattering not specified"),
    # ===== C. beam-conditioning assembly (representative: EKSMA Optics 160-0081, 8x Galilean) =====
    "expander_expansion_ratio": Param(8.0, "-", S.MANUFACTURER_SPECIFICATION,
                                      "8x. Recorded only: the power budget has no input-beam diameter, and clear-aperture clipping and "
                                      "beam-size matching to the receiver are not modelled", used=False),
    "expander_material": Param("AR-coated fused silica", "-", S.MANUFACTURER_SPECIFICATION, "", used=False),
    "expander_arrangement": Param("Galilean", "-", S.MANUFACTURER_SPECIFICATION, "", used=False),
    "expander_num_lenses": Param(2, "-", S.DESIGN_DERIVED,
                                 "conventional Galilean two-lens architecture; not a published internal prescription", used=False),
    "expander_num_surfaces": Param(4, "-", S.DESIGN_DERIVED,
                                   "two lenses x two surfaces; not a published internal prescription", used=False),
    "conditioning_transmission": Param(0.96, "-", S.MANUFACTURER_SPECIFICATION,
                                       "stated > 96 % for the assembly; 0.960 used as the conservative baseline; assembly-level total"),
    "expander_surface_losses": Param(None, "-", S.UNRESOLVED,
                                     "per-surface reflection/absorption/scatter unavailable; the 4 % is NOT split among the surfaces"),
    # ===== D. fused-silica / synthetic-quartz window (representative: Shimadzu 1064 nm laser window) =====
    "window_material": Param("synthetic quartz", "-", S.MANUFACTURER_SPECIFICATION, "", used=False),
    "window_incidence_deg": Param(0.0, "deg", S.MANUFACTURER_SPECIFICATION, "normal incidence", used=False),
    "window_transmission": Param(0.996, "-", S.MANUFACTURER_SPECIFICATION,
                                 "at least 99.6 %. Shimadzu calculates it from the single-sided reflectance"),
    "window_single_sided_reflectance": Param(0.002, "-", S.MANUFACTURER_SPECIFICATION,
                                             "< 0.2 % per side. Basis of the transmission figure: NOT added as an extra loss", used=False),
    "window_absorption_coefficient": Param(None, "1/m", S.UNRESOLVED, "bulk absorption"),
    "window_scattering": Param(None, "-", S.UNRESOLVED, "scattering"),
    "window_thickness": Param(None, "m", S.UNRESOLVED, "actual window thickness for this spacecraft"),
    "window_mode": Param("catalogue", "-", S.ASSUMED_MODELLING_CHOICE,
                         "'catalogue' (use the specified total T) or 'uncoated_fresnel' (comparison case only: the representative "
                         "Shimadzu window is AR-coated)"),
}

# Numerical settings (receiver quadrature). validation.py shows results are insensitive to them.
NUMERICS = {"n_panels": 4, "n_nodes": 64, "n_phi": 256}


def value(name: str):
    """Return a parameter value; refuse UNRESOLVED ones instead of guessing."""
    p = PARAMS[name]
    if p.value is None:
        raise ValueError(f"'{name}' is UNRESOLVED ({p.note or 'no source'}); "
                         "supply a sourced value in config.PARAMS first.")
    return p.value


def unresolved_names() -> list[str]:
    return [k for k, p in PARAMS.items() if p.status is Status.UNRESOLVED]


# ---------------------------------------------------------------------------
# Input validation helpers
# ---------------------------------------------------------------------------
def require_positive(name: str, x: float) -> None:
    if not isinstance(x, (int, float)) or not math.isfinite(x) or x <= 0.0:
        raise ValueError(f"{name} must be finite and > 0, got {x!r}")


def require_nonnegative(name: str, x: float) -> None:
    if not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0.0:
        raise ValueError(f"{name} must be finite and >= 0, got {x!r}")


def require_fraction(name: str, x: float, allow_zero: bool = False) -> None:
    lo_ok = x >= 0.0 if allow_zero else x > 0.0
    if not isinstance(x, (int, float)) or not math.isfinite(x) or not lo_ok or x > 1.0:
        rng = "[0, 1]" if allow_zero else "(0, 1]"
        raise ValueError(f"{name} must be in {rng}, got {x!r}")


@dataclass(frozen=True)
class SystemInputs:
    """All inputs of one model evaluation, in SI units. Immutable and validated."""
    power_w: float                      # PEAK optical power of the reference scenario
    wavelength_m: float
    waist_m: float
    m2: float
    distance_m: float
    aperture_radius_m: float
    pointing_arcsec: float
    mirror_R: float
    cond_T: float
    window_T: float
    prop_loss_fraction: float = 0.0
    window_mode: str = "catalogue"
    window_thickness_m: float | None = None
    window_alpha_per_m: float | None = None
    # typical combined absorption+scatter (upper bound), optional
    mirror_abs_scatter_ppm: float | None = None
    # optional; enables an average-power line
    duty_cycle: float | None = None

    def __post_init__(self):
        require_positive("transmitter power [W]", self.power_w)
        require_positive("wavelength [m]", self.wavelength_m)
        require_positive("beam waist radius [m]", self.waist_m)
        require_positive("M^2", self.m2)
        if self.m2 < 1.0:
            raise ValueError(
                f"M^2 must be >= 1 (diffraction limit), got {self.m2}")
        require_positive("propagation distance [m]", self.distance_m)
        require_positive(
            "receiver aperture radius [m]", self.aperture_radius_m)
        require_nonnegative("pointing error [arcsec]", self.pointing_arcsec)
        require_fraction("mirror reflectance", self.mirror_R)
        require_fraction("beam-conditioning transmission", self.cond_T)
        require_fraction("window transmission", self.window_T)
        require_fraction("propagation loss fraction",
                         self.prop_loss_fraction, allow_zero=True)
        if self.window_mode not in ("catalogue", "uncoated_fresnel"):
            raise ValueError(
                f"window_mode must be 'catalogue' or 'uncoated_fresnel', got {self.window_mode!r}")
        if self.window_thickness_m is not None:
            require_positive("window thickness [m]", self.window_thickness_m)
        if self.window_alpha_per_m is not None:
            require_nonnegative(
                "window absorption coefficient [1/m]", self.window_alpha_per_m)
        if self.mirror_abs_scatter_ppm is not None:
            require_nonnegative(
                "mirror combined absorption+scatter [ppm]", self.mirror_abs_scatter_ppm)
            if self.mirror_abs_scatter_ppm * 1e-6 > (1.0 - self.mirror_R) + 1e-15:
                raise ValueError(
                    "mirror absorption+scatter exceeds the total mirror loss 1 - R")
        if self.duty_cycle is not None:
            require_fraction("duty cycle", self.duty_cycle)

    @property
    def pointing_rad(self) -> float:
        return self.pointing_arcsec * ARCSEC_TO_RAD

    @classmethod
    def from_config(cls) -> "SystemInputs":
        """Baseline v1 from config.PARAMS (guaranteed specifications; optional unresolved items stay None)."""
        def opt(k): return PARAMS[k].value  # may legitimately be None
        return cls(
            power_w=value("transmitter_power"),
            wavelength_m=value("wavelength"),
            waist_m=value("beam_waist_radius"),
            m2=value("m_squared"),
            distance_m=value("distance"),
            aperture_radius_m=value("receiver_radius"),
            pointing_arcsec=value("pointing_error"),
            mirror_R=value("mirror_reflectance"),
            cond_T=value("conditioning_transmission"),
            window_T=value("window_transmission"),
            prop_loss_fraction=value("propagation_loss_fraction"),
            window_mode=value("window_mode"),
            window_thickness_m=opt("window_thickness"),
            window_alpha_per_m=opt("window_absorption_coefficient"),
            duty_cycle=opt("duty_cycle"),
        )

    def with_(self, **changes) -> "SystemInputs":
        """Return a modified, re-validated copy."""
        return replace(self, **changes)
