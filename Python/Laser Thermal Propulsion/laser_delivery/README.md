# Laser-delivery model for externally laser-heated (LSP) thermal propulsion

Modular Python model of a space-based laser link: transmitter -> vacuum Gaussian
propagation -> pointing error -> receiver aperture -> receiving mirror ->
beam-conditioning assembly -> fused-silica window -> **power at the LSP chamber entrance**.
It deliberately stops there: no plasma, nozzle, thrust or trajectory model, and no
magnetic nozzle is assumed. Vacuum is assumed; a propagation-loss hook exists for a
future atmosphere/medium model.

## Run
```
pip install -r requirements.txt
python main.py                       # report, CSVs and 5 figures in ./outputs
python main.py --window-mode uncoated_fresnel   # physical uncoated-window case
python -m unittest discover -s tests -t . -v    # 52 tests (also runs under pytest)
```
All inputs are edited **only** in `config.py` (`PARAMS`). Each carries a classification (register v1):
`PUBLISHED_LITERATURE`, `MANUFACTURER_SPECIFICATION`, `PROVISIONAL_MANUFACTURER_MEASUREMENT`,
`ASSUMED_MODELLING_CHOICE`, `DESIGN_DERIVED`, `UNVERIFIED` (illustrative), `UNRESOLVED` (no value; requesting it
raises an error). `references.py` builds the source register from `config.PARAMS`, and a test keeps them in agreement.
Alternative cases live in `scenarios.py`; the reference scenario `baseline_v1` is never overwritten.

## Files
| File | Responsibility |
|---|---|
| `config.py` | Parameter table + provenance, validated `SystemInputs`, numerical settings |
| `laser_source.py` | Validated transmitter object: power, wavelength, waist, M^2, BPP, z_R, divergence |
| `beam_propagation.py` | z_R, divergence, w(z), diameter, I0, intensity profile (M^2-corrected) |
| `pointing_error.py` | arcsec -> rad -> lateral displacement d = z tan(delta) |
| `receiver.py` | Aperture capture: production radial quadrature, independent 2-D polar quadrature, analytic non-central chi-square reference, convergence study, jitter-averaged variant |
| `optical_components.py` | Lumped (assembly-level) components; physical Fresnel window; Sellmeier index of fused silica |
| `optical_losses.py` | Series optical train, per-stage energy-conservation enforcement, cumulative table |
| `power_budget.py` | End-to-end series budget; every stage fed by the previous stage's output |
| `uncertainty.py` | Local elasticities, illustrative tornado, seeded Monte Carlo |
| `validation.py` | 34 runtime checks (math / consistency / literature spot-check / input rejection) |
| `references.py` | Source register built from `config.PARAMS` -> `outputs/references.csv` |
| `scenarios.py` | Named alternative cases, pointing and M^2 sensitivity sweeps |
| `report.py`, `plots.py`, `main.py` | Report, CSV, figures, orchestration |
| `tests/` | 52 independent unit tests |

## Equations (SI; derivations in module docstrings)
* BPP = M^2 lambda / pi = w0 theta
* z_R = pi w0^2 / (M^2 lambda); w(z) = w0 sqrt(1 + (z/z_R)^2); theta = M^2 lambda / (pi w0)
* I(r) = (2P / pi w^2) exp(-2r^2/w^2); I0 = 2P / (pi w^2)
* Pointing: d = z tan(delta)
* Capture, centred: eta = 1 - exp(-2a^2/w^2). Displaced: integral of I over the disc, evaluated as
  eta = (4/w^2) int_0^a r exp(-2(r-d)^2/w^2) i0e(4rd/w^2) dr (angular integral done analytically), cross-checked
  by a full 2-D quadrature and by eta = ncx2.cdf((2a/w)^2, 2, (2d/w)^2).
* Uncoated window: R = ((n-1)/(n+1))^2, T = (1-R)^2 tau / (1 - R^2 tau^2), tau = exp(-alpha t)
* Conservation per stage: P_in = P_out + modelled losses + combined-unresolved loss

## Baseline result (hand-checkable)
w = 19.18 m at 1,200 km (z/z_R = 76.7), d = 0.582 m, capture 12.68 % (3.804 kW), optical train 0.9542,
**delivered 3.630 kW (12.10 % overall)** in the **peak-power** reference scenario. The 30 kW is a peak value from the
cited study; CW/average power and duty cycle are unresolved, so nothing here claims sustained 30 kW operation. The
receiver geometry (a/w = 0.26) dominates; the optical train costs 0.58 % of transmitter power. Pointing at 0.1 arcsec
costs only 0.02 percentage points of capture because d << w. Alternative cases: typical-mirror (3,637 W), uncoated-window
comparison, and pointing / M^2 sensitivity sweeps (neither has an established source).

## Parameter classification (register v1; as supplied by the user)
* **Published literature (provisional):** 30 kW peak power, 1064 nm, 0.25 m waist (MacRobbie et al. 2025).
* **Unverified:** M^2 = 11.8 (assumed sensitivity case) and 0.1 arcsec pointing (illustrative).
* **Assumed modelling choices:** 1,200 km distance, 5 m receiver radius, vacuum.
* **Manufacturer specification (guaranteed):** mirror R > 99.8 %, expander T > 96 %, window T >= 99.6 % (used as baseline).
* **Typical manufacturer measurement (provisional):** mirror R > 99.98 %, combined absorption+scatter < 150 ppm (scenario only).
* **Design-derived:** two lenses, four surfaces (not a published prescription).
* **Unresolved:** CW/average power, duty cycle, mirror and window absorption/scatter split, per-surface expander losses,
  window thickness. None is used as data.
* The expander's single 4 % loss is used once as an assembly total; the window's single-sided reflectance is the basis of its
  total and is never added again.

## Key assumptions (all adjustable or flagged)
1. M^2 is a second-moment value, identical in x and y; Gaussian transverse profile (an assumption for a multimode beam).
2. 0.1 arcsec is a fixed offset; a jitter alternative is reported as information only.
3. Receiver aperture is hard-edged, unobstructed, perpendicular to the beam, centred on the nominal axis.
4. Representative catalogue parts (45 deg dielectric mirror, 8x compact expander, 1064 nm window) are stand-ins; their
   applicability to a large spacecraft optical train (aperture size, AOI, beam-size matching, clear-aperture clipping) is
   not established. The expander's x8 magnification is recorded but not used in the power budget.
5. Vacuum; no turbulence, thermal blooming, aberrations, point-ahead or platform motion. Vacuum does not imply that
   pointing jitter or wavefront error are negligible.

## Data sources and status
See `references.py` / `outputs/references.csv`. The classifications and URLs come from the user-supplied register; they were
**not independently re-fetched** here. Only Malitson (1965) was looked up (DOI 10.1364/josa.55.001205; Sellmeier
coefficients from a secondary implementation, to be confirmed). Siegman (1986) and ISO 11146-1 are cited from recall.

## What is and is not verified
* **Verified mathematically:** implementation of the equations; agreement of three independent capture
  evaluations to ~1e-16; quadrature convergence and grid independence; energy conservation at every stage;
  physical bounds; correct rejection of impossible inputs; far-field elasticities consistent with theory.
* **Checked against literature:** only a spot-check of fused-silica index (limited).
* **NOT validated:** no experimental or flight data. The Gaussian-profile assumption for M^2 = 11.8, the pointing
  model, the catalogue optical values and every unresolved input still need validation. "Tests pass" means the
  code does what the equations say, not that the equations describe the real system.

## Sensitivity and uncertainty
Elasticities (d ln P / d ln x) are mathematical. Ranges in `uncertainty.ILLUSTRATIVE_RANGES` are
**illustrative, not measured**; replace them with real tolerances when available. Model-form uncertainty is not quantified.

## Limitations / next steps
Obtain real datasheets (mirror, expander, window); define M^2 and pointing statistics from the source study;
add transmitter-side aperture clipping and jitter/point-ahead if needed; consider a non-Gaussian profile
(e.g. flat-top) sensitivity study; add medium attenuation only if an atmosphere is introduced.
