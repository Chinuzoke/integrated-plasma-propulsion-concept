"""
laser_source.py -- transmitter power, wavelength and beam characteristics.

Immutable, validated description of the transmitter. Beam parameter product
BPP = M^2 lambda / pi = w0 * theta is conserved in free propagation.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import beam_propagation as bp
import config as cfg


@dataclass(frozen=True)
class LaserSource:
    power_w: float
    wavelength_m: float
    waist_radius_m: float
    m_squared: float

    def __post_init__(self):
        cfg.require_positive("power [W]", self.power_w)
        cfg.require_positive("wavelength [m]", self.wavelength_m)
        cfg.require_positive("waist radius [m]", self.waist_radius_m)
        cfg.require_positive("M^2", self.m_squared)
        if self.m_squared < 1.0:
            raise ValueError(f"M^2 must be >= 1 (diffraction limit), got {self.m_squared}")

    @property
    def beam_parameter_product_m_rad(self) -> float:
        return bp.beam_parameter_product(self.wavelength_m, self.m_squared)

    @property
    def rayleigh_range_m(self) -> float:
        return bp.rayleigh_range(self.waist_radius_m, self.wavelength_m, self.m_squared)

    @property
    def divergence_half_angle_rad(self) -> float:
        return bp.divergence_half_angle(self.waist_radius_m, self.wavelength_m, self.m_squared)

    @classmethod
    def from_inputs(cls, inp: cfg.SystemInputs) -> "LaserSource":
        return cls(inp.power_w, inp.wavelength_m, inp.waist_m, inp.m2)

    @classmethod
    def from_config(cls) -> "LaserSource":
        return cls.from_inputs(cfg.SystemInputs.from_config())
