"""
references.py -- traceable parameter/source register (register v1).

Values, units, classifications and notes are read from config.PARAMS so the
register cannot drift from the model inputs; this file adds only the SOURCES.
tests/test_register.py enforces that every config parameter has a record.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, asdict

import config as cfg

_MACROBBIE = ("Space-Based Solar Power: Implications for Operational Robustness in Lunar EVAs and Exploration "
              "Architectures", "MacRobbie et al.", "2025", "(journal: Aerotecnica Missili & Spazio)")
_MIRROR = ("Excelitas/LINOS DLHS NIR 1064 nm dielectric-coated plane mirror (product page + datasheet)",
           "Excelitas / LINOS", "n/a",
           "https://prod.excelitas.com/product/dielectric-coated-plane-mirrors-dlhs-nir-1064-nm ; "
           "https://linosoptics.excelitas.com/out/Graphics/en/00125299_0.pdf")
_EXPANDER = ("EKSMA Optics 160-0081 compact 8x Galilean beam expander (product page + brochure)", "EKSMA Optics", "n/a",
             "https://www.eksmaoptics.com/v/160-0081/2869 ; "
             "https://eksmaoptics.com/out/media/EKSMA_Optics_Zoom_Beam_Expander.pdf")
_WINDOW = ("Shimadzu laser windows for 1064 nm (product page)", "Shimadzu", "n/a",
           "https://www.shimadzu.com/opt/products/highpower/o-k25cur0000007rp3.html")
_NONE_ASSUMED = ("none: modelling choice / scenario input",
                 "this model", "n/a", "n/a")
_NONE_UNRES = ("none: unresolved", "n/a", "n/a", "n/a")
_NONE_UNVER = ("none established", "n/a", "n/a", "n/a")
_DERIVED = ("derived from conventional two-lens Galilean architecture (not a manufacturer prescription)",
            "this model", "n/a", "n/a")

SOURCES = {
    "transmitter_power": _MACROBBIE, "wavelength": _MACROBBIE, "beam_waist_radius": _MACROBBIE,
    "cw_average_optical_power": _NONE_UNRES, "duty_cycle": _NONE_UNRES,
    "m_squared": _NONE_UNVER, "pointing_error": _NONE_UNVER,
    "distance": _NONE_ASSUMED, "receiver_radius": _NONE_ASSUMED, "propagation_medium": _NONE_ASSUMED,
    "propagation_loss_fraction": _NONE_ASSUMED, "window_mode": _NONE_ASSUMED,
    "optical_component_wavelength": ("Manufacturer datasheets (mirror, expander, window)", "Excelitas, EKSMA, Shimadzu", "n/a",
                                     _MIRROR[3] + " ; " + _EXPANDER[3] + " ; " + _WINDOW[3]),
    "mirror_aoi_deg": _MIRROR, "mirror_reflectance": _MIRROR, "mirror_typical_reflectance": _MIRROR,
    "mirror_typical_abs_scatter_ppm": _MIRROR,
    "mirror_absorption_fraction": _NONE_UNRES, "mirror_scatter_fraction": _NONE_UNRES,
    "expander_expansion_ratio": _EXPANDER, "expander_material": _EXPANDER, "expander_arrangement": _EXPANDER,
    "expander_num_lenses": _DERIVED, "expander_num_surfaces": _DERIVED, "conditioning_transmission": _EXPANDER,
    "expander_surface_losses": _NONE_UNRES,
    "window_material": _WINDOW, "window_incidence_deg": _WINDOW, "window_transmission": _WINDOW,
    "window_single_sided_reflectance": _WINDOW,
    "window_absorption_coefficient": _NONE_UNRES, "window_scattering": _NONE_UNRES, "window_thickness": _NONE_UNRES,
}

# Literature used by the physics (not model inputs)
SUPPORTING = [
    ("fused_silica_refractive_index(lambda) [Sellmeier]", "Interspecimen Comparison of the Refractive Index of Fused Silica",
     "I. H. Malitson", "1965", "https://doi.org/10.1364/josa.55.001205",
     "Paper and DOI verified by lookup. Sellmeier coefficients taken from the PyXOpto implementation "
     "(secondary source): confirm against the paper. Used only for the uncoated-window comparison case."),
    ("gaussian_beam_M2_formalism", "Lasers",
     "A. E. Siegman", "1986", "ISBN 0-935702-11-5"),
    ("beam_width_and_M2_definitions", "ISO 11146-1 (laser beam widths, divergence angles, beam propagation ratios, Part 1)",
     "ISO", "2005", "not recorded", "Title verified."),
]


@dataclass(frozen=True)
class Ref:
    parameter: str
    value: str
    unit: str
    classification: str
    provisional: str
    used_in_calculation: str
    source_title: str
    author_or_manufacturer: str
    year: str
    url_or_doi: str
    notes: str


def build_register() -> list:
    reg = []
    for name, p in cfg.PARAMS.items():
        title, who, year, url = SOURCES[name]
        reg.append(Ref(name, "UNRESOLVED" if p.value is None else str(p.value), p.unit, p.status.name,
                       "yes" if p.provisional else "no", "yes" if p.used else "recorded only",
                       title, who, year, url, p.note))
    for name, title, who, year, url, note in SUPPORTING:
        reg.append(Ref(name, "(supporting physics)", "-", "PUBLISHED_LITERATURE", "no", "yes (window comparison)"
                       if "refractive" in name else "reference", title, who, year, url, note))
    return reg


def write_csv(path: str) -> None:
    reg = build_register()
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(asdict(reg[0]).keys()))
        wr.writeheader()
        for r in reg:
            wr.writerow(asdict(r))
