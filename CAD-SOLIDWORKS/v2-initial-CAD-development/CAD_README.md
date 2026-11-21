# CAD package: Model III-inspired baseline (preliminary)

**Status:** preliminary reference geometry. Not pressure-certified, not flight-qualified, not manufacturing-ready. CAD Kernel (SolidWorks) was used. 

Axes: These are yet to be defined.

## Folder structure
```
CAD-Solidworks/v2-initial-CAD-development/components/<part>_{front}.SLDPRT   (20 parts) 
CAD-Solidworks/v2-initial-CAD-development/Assembly/<part>_{views}.jpg (8 views)
CAD-Solidworks/v2-initial-CAD-development/Drawings/(To be added)
CAD-Solidworks/v2-initial-CAD-development/Renders/(To be added)
CAD-Solidworks/v2-initial-CAD-development/CAD_README.md 
```


## Drawing sheets (reference-sheet layout)
To be done.

## Published versus proposed dimensions
| Parameter | Value | Status | Basis |
|---|---|---|---|
| throat_d | 0.5 mm | published | S01 Table 4 / Sect. 3.2 |
| lsp_to_throat | 4 mm | published | S01 Sect. 4.2 |
| lens_d | 25.4 mm | published | S01 Sect. 3.1 (1-inch lens) |
| lens_efl | 125 mm | published | S01 Table 3 (F12.5 = baseline for thruster tests) |
| beam_d | 10 mm | published | S01 Sect. 3.1 (about 10 mm) |
| electrode_d | 3.2 mm | published | S01 Sect. 3.2 (3.2 mm W + 2% ThO2 rod) |
| plenum_bore_D | 18 mm | published | S01 Fig. 4b (D = 18 mm); annulus area 27.5 mm^2 reproduced |
| insert_od | 17 mm | published | S01 Fig. 4b (d = 17 mm), interpreted as insert outer diameter |
| generation_d | 6 mm | published (Model II); proposed for Model III | S01 Sect. 4.4 gives 6 mm for the generation section of Model II |
| electrode_z | -30 mm | proposed | carries the Model I electrode-to-throat distance (S01 Sect. 3.2.1) to Model III; not published for Model III |
| chamber_end_z | -7 mm | proposed | end of the 6 mm chamber; sets converging section; must be set from CFD |
| confinement_d | 1.5 mm | proposed | narrow confinement channel; diameter not published |
| confinement_z0 | -3 mm | proposed | start of confinement channel (3 mm length with end at -1... see geometry) |
| throat_land | 0.5 mm | proposed | cylindrical throat land length |
| diverge_half_angle | 15 deg | proposed | conical diverging section |
| exit_d | 1 mm | proposed | area ratio 4 against the 0.5 mm throat; nozzle shape not published |
| body_od | 50 mm | proposed | main chamber body |
| body_z0 | -74 mm | proposed | body upstream face |
| body_z1 | 16 mm | proposed | body downstream face |
| window_d | 25 mm | proposed | laser window diameter |
| window_t | 5 mm | proposed | laser window thickness |
| lens_n | 1.5 - | proposed | refractive index assumed to derive curvature and back focal distance; lens glass not published |
| lens_edge_t | 3 mm | proposed | lens edge thickness |
| tip_gap | 2 mm | proposed | electrode tip-to-tip gap; not published |
| flange_od | 80 mm | proposed | chamber mounting flange |
| plate_w | 120 mm | proposed | thruster mounting plate |

Derived: lens back focal distance 122.13 mm and curvature radius 62.5 mm from the assumed refractive index 1.5 (a proposed assumption). Beam cone check (derived): half-angle = (beam_d/2)/EFL = 0.04 rad; at 60 mm from focus the beam diameter is 4.8 mm, inside the 6 mm bore and the 20 mm retaining-ring aperture. Annulus area between the published D18 and D17: 27.5 mm^2, matching the paper.

## Components
Material entries are candidates only; no numerical material property is used anywhere in this package, so no property source applies yet and datasheets must be selected before analysis. Surface finish, thread and tolerance values are not assigned except where stated; fits and GD&T are to be set with the mechanical analysis (datum A = beam axis; faces perpendicular to A). Inspection: dimensional check of critical features, leak test and optical inspection after assembly (to be specified).

### main_chamber_body 
- Function: Holds the plenum bore, gas manifold, ports and the nozzle housing; carries pressure loads.
- Geometry and dimensions: L 90, OD 50, bore D18 (pub), pocket D26, bore D10 at nozzle, cooling relief OD 46 (all proposed except D18)
- Material recommendation: 316L-class stainless steel or a copper alloy where heat flux demands; selection pending thermal analysis
- Manufacturing: turning and milling; cross-drilled ports
- Critical dimensions and datums: datum A = bore axis (beam axis); concentricity of D26, D18, D10 bores to A critical
- Unresolved: clearance 0.5 mm radial (D18/D17) is a published-dimension consequence; locating pilot for the insert is NOT defined

### lsp_sustaining_insert
- Function: Replaceable insert: D6 generation section, converging section, narrow confinement channel, electrode and feed holes.
- Geometry and dimensions: L 61, OD 17 (pub), D6, converging 4 long, D1.5 confinement (proposed), 4 feed holes D2, electrode holes D3.4
- Material recommendation: refractory ceramic (boron nitride or alumina) or high-temperature alloy; published experiments do not state the material
- Manufacturing: turning; ceramic grinding or EDM for holes
- Critical dimensions and datums: coaxial with datum A; D1.5 bore finish critical (flow stability)
- Unresolved: located only by the shoulder at z = -1; radial location and anti-rotation are unresolved

### copper_throat_insert
- Function: Converging section, 0.5 mm throat and diverging section.
- Geometry and dimensions: L 4, OD 10, throat D0.5 (pub), land 0.5, half-angle 15 deg, exit D1.0 (proposed; area ratio 4)
- Material recommendation: copper (published for the throat); alternatives to assess: tungsten-copper, molybdenum alloys, tungsten
- Manufacturing: turning then wire/sinker EDM or micro-drilling for the throat
- Critical dimensions and datums: throat D0.5 and roundness; coaxial with A
- Unresolved: copper erosion, oxidation and softening at the throat unresolved; discharge coefficient unknown

### nozzle_retainer
- Function: Threaded retainer that preloads the throat insert and sets the exit opening.
- Geometry and dimensions: L 13, collar OD 14, sleeve OD 10, bore D6 (proposed)
- Material recommendation: stainless steel
- Manufacturing: turning, thread (to be specified)
- Critical dimensions and datums: bore coaxial with A
- Unresolved: thermal growth against copper unresolved

### window_retaining_ring
- Function: Retains the window stack; aperture D20 (beam at this plane about 5.6 mm).
- Geometry and dimensions: t 4, OD 26, aperture D20 (proposed)
- Material recommendation: stainless steel or aluminium
- Manufacturing: turning, thread to be specified
- Critical dimensions and datums: face flatness against the window
- Unresolved: torque limit must be defined: fused silica must not be loaded by excessive ring force

### laser_window
- Function: Laser induction window.
- Geometry and dimensions: D25, t 5 (proposed); anti-reflection coating for 1080 nm (proposed)
- Material recommendation: fused silica
- Manufacturing: purchased optic
- Critical dimensions and datums: surface quality from the optical design
- Unresolved: edge stress, stress birefringence and thermal lensing unevaluated

### window_seal
- Function: Seals the window face against the seat.
- Geometry and dimensions: ID 18, OD 24, t 1 (proposed)
- Material recommendation: metal C-seal or a cooled seal after the seal temperature is checked; no elastomer near the plasma without that check
- Manufacturing: purchased or stamped
- Critical dimensions and datums: seat face flatness
- Unresolved: seal temperature and chemical environment unknown

### window_seat
- Function: Seat ring supporting the seal and window.
- Geometry and dimensions: t 2, OD 25.9, ID 14 (proposed)
- Material recommendation: soft metal (e.g. copper) or compliant ceramic
- Manufacturing: turning
- Critical dimensions and datums: parallelism to the window
- Unresolved: contact stress on the glass unresolved

### focusing_lens 
- Function: Focuses the collimated 10 mm beam onto the LSP.
- Geometry and dimensions: D25.4 (pub, 1 inch), EFL 125 (pub), R = 62.5 and centre thickness 4.30 derived from an assumed n = 1.5, edge 3.0 (proposed)
- Material recommendation: fused silica or similar with 1080 nm coating (lens glass not published)
- Manufacturing: purchased optic
- Critical dimensions and datums: position along z sets the focus 4 mm from the throat
- Unresolved: back focal distance derived: 122.13 mm; actual lens specification needed

### lens_holder
- Function: Holds the lens (flat face on the internal shoulder).
- Geometry and dimensions: L 14, OD 34, bore D25.8 (proposed)
- Material recommendation: aluminium, anodised
- Manufacturing: turning
- Critical dimensions and datums: bore coaxial with A
- Unresolved: lens retention ring is not defined

### optical_alignment_mount 
- Function: Clamp ring on a post and base for lens/laser-head alignment.
- Geometry and dimensions: ring OD 48, post D14, base 70 x 70 x 8 (proposed)
- Material recommendation: aluminium
- Manufacturing: milling and turning
- Critical dimensions and datums: axis height set by datum A
- Unresolved: real alignment needs translation stages; slit, clamp screw and stage not modelled

### beam_axis_reference
- Function: Virtual datum line of the laser beam axis (not a physical part).
- Geometry and dimensions: length 100, D0.7 display only
- Material recommendation: none
- Manufacturing: none
- Critical dimensions and datums: datum A
- Unresolved: display aid only

### gas_inlet_fitting
- Function: Argon feed to the manifold.
- Geometry and dimensions: tube OD 6, ID 4, L 14 (proposed)
- Material recommendation: brass or stainless fitting
- Manufacturing: purchased standard fitting
- Critical dimensions and datums: -
- Unresolved: thread and seal type not defined

### pressure_port_fitting
- Function: Plenum pressure measurement port.
- Geometry and dimensions: tube OD 6, ID 4, L 14; D4 port at the manifold (proposed)
- Material recommendation: brass or stainless
- Manufacturing: purchased
- Critical dimensions and datums: -
- Unresolved: S01 measured the 27.5 mm^2 section; the port position here is proposed

### diagnostic_port_plug
- Function: Radial viewport to the LSP region (proposed engineering addition).
- Geometry and dimensions: OD 6 plug, flange OD 12, central D2 view bore (proposed)
- Material recommendation: stainless steel with a window material to be chosen
- Manufacturing: turning
- Critical dimensions and datums: -
- Unresolved: window seal and the cross-hole through the insert weaken the design; unresolved

### ignition_electrode
- Function: Arc ignition electrode (cathode at -y, anode at +y; polarity unresolved).
- Geometry and dimensions: D3.2 (pub), L 60, tip gap 2.0 (proposed)
- Material recommendation: published: tungsten with 2% thoria; a non-thoriated tungsten (ceriated or lanthanated) should be assessed
- Manufacturing: ground rod
- Critical dimensions and datums: tip geometry and gap
- Unresolved: ignition performance of the non-thoriated alternative unresolved

### electrode_insulator_feedthrough
- Function: Insulates and seals the electrode through the body wall.
- Geometry and dimensions: tube OD 9, ID 3.4, L 30, cap OD 14 (proposed)
- Material recommendation: alumina or a machinable ceramic
- Manufacturing: grinding
- Critical dimensions and datums: coaxial with the port
- Unresolved: high-voltage rating, seal and temperature unresolved

### cooling_jacket
- Function: Coolant sleeve around the nozzle section (2 mm annular gap).
- Geometry and dimensions: L 16, OD 58, ID 50, two D4 ports (proposed)
- Material recommendation: aluminium or stainless steel
- Manufacturing: turning, drilling
- Critical dimensions and datums: fit to body OD 50
- Unresolved: coolant, seals (away from the plasma), flow rate and heat duty not defined

### coolant_fitting
- Function: Coolant inlet and outlet.
- Geometry and dimensions: tube OD 6, ID 4, L 14 (proposed)
- Material recommendation: brass or stainless
- Manufacturing: purchased
- Critical dimensions and datums: -
- Unresolved: leak path at the jacket interface unresolved

### chamber_mounting_flange
- Function: Flange clamping the body to the mounting plate.
- Geometry and dimensions: OD 80, t 6, bore D50, 6 x D4.5 on D64 (proposed)
- Material recommendation: aluminium or stainless
- Manufacturing: turning, drilling
- Critical dimensions and datums: faces square to A
- Unresolved: joint to the body (fit, adhesive or weld) not defined

### thruster_mounting_plate
- Function: Support plate for the flange.
- Geometry and dimensions: 120 x 120 x 8, bore D50, 6 x D3.3 (proposed)
- Material recommendation: aluminium
- Manufacturing: milling
- Critical dimensions and datums: square to A
- Unresolved: support to the optical table not defined

### fastener_m4
- Function: Socket-head cap screw for the flange joint (6 used).
- Geometry and dimensions: M4 x 12 shank, head D7 x 4 (proposed)
- Material recommendation: A2/A4 stainless, ISO 4762 class
- Manufacturing: purchased
- Critical dimensions and datums: -
- Unresolved: preload and torque not calculated

## Assembly sequence (proposed yet-to-be-done)
1. Insert the sustaining insert; check the shoulder contact z = 1.
2. Fit seat ring, seal, laser window and retaining ring (torque-limited) into the body pocket.
3. Fit the copper throat insert, then the nozzle retainer.
4. Fit electrode insulators and electrodes; set the tip gap.
5. Fit the cooling jacket and fittings; fit gas, pressure and diagnostic ports.
6. Clamp the flange to the plate (6 x M4).
7. Mount and align lens, holder and alignment mount to place the focus 4 mm upstream of the throat.

## Preliminary hand checks (CALCULATED, plenum pressure 1.2 MPa from Fig. 11 of the benchmark paper; not an analysis)
- Thick-wall cylinder (OD 50, bore 18): hoop stress at the bore 1.56 MPa, before stress concentration at the cross-ports (not evaluated).
- Window as a clamped flat plate, unsupported radius 10 mm, t 5 mm: peak stress about 3.6 MPa; a simply supported edge gives roughly 1.55 to 1.65 times that. Axial load on the window about 416 N at a mean seal diameter of 21 mm.
- Not evaluated: thermal stress, thermal-expansion mismatch (copper, steel, glass), nozzle thermal loading, chamber deformation, fatigue and thermal cycling, seal compression, fastener preload, vibration, optical alignment stability, coolant leakage.

## Unresolved mechanical risks
1. lsp sustaining insert internal dimensions are unpublished; confinement diameter, converging angle, nozzle contour and insert length are placeholders.
2. The insert is not located radially (0.5 mm gap) and has no anti-rotation feature.
3. Window purge path is not modelled; the paper states gas was directed onto the window but not how.
4. Window seal type and temperature, and retaining-ring load on fused silica, are undefined.
5. Throat material (copper) wear, thermal growth and joining unresolved; high-temperature alternatives not yet assessed.
6. Thoriated electrode: a non-thoriated alternative is recommended but its ignition performance is unknown; electrode polarity and gap are proposed.
7. The diagnostic viewport breaches the insert and body; sealing unresolved.
8. Optical alignment depends on translation stages that are not modelled.
9. All pressure-containing parts require full mechanical analysis (including ASME-type pressure-vessel assessment where applicable) before manufacture; no compliance is claimed.
