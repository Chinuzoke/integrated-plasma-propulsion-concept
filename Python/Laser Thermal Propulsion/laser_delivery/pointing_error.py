"""
pointing_error.py -- angular pointing error -> lateral beam displacement.

A transmitter pointing error `delta` (rad) displaces the beam centroid at the
receiver plane by

    d = z * tan(delta)      [m]      (~ z * delta for delta << 1)

Assumptions: the error is a fixed, deterministic angular offset of the beam
axis about the transmitter; the receiver is perpendicular to the nominal line
of sight and its position is known exactly. Random jitter, receiver-side
tracking error and platform motion during light-travel time (point-ahead)
are NOT included; see receiver.mean_capture_with_jitter for a statistical
variant.
"""
from __future__ import annotations

import math

import config as cfg


def arcsec_to_rad(arcsec: float) -> float:
    cfg.require_nonnegative("angle [arcsec]", arcsec)
    return arcsec * cfg.ARCSEC_TO_RAD


def lateral_displacement(distance_m: float, pointing_rad: float) -> float:
    """d = z tan(delta) [m]."""
    cfg.require_positive("distance", distance_m)
    cfg.require_nonnegative("pointing angle [rad]", pointing_rad)
    if pointing_rad >= math.pi / 2:
        raise ValueError("pointing angle must be < pi/2 rad")
    return distance_m * math.tan(pointing_rad)


def displacement_from_arcsec(distance_m: float, arcsec: float) -> float:
    return lateral_displacement(distance_m, arcsec_to_rad(arcsec))
