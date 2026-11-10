# Integrated Plasma Propulsion System

**An integrated computational and engineering investigation of plasma-based propulsion, combining laser delivery, plasma modelling, thruster design, multiphysics analysis, and propulsion-system modelling.**

**Project status:** Active Development

## 1. Project Overview

This project aims to develop an integrated computational framework for investigating plasma-based propulsion systems.

The work combines numerical modelling of laser delivery and plasma generation with physical thruster design, engineering analysis, and system-level performance evaluation.

Rather than treating these disciplines as independent exercises, the project seeks to establish a connected engineering workflow in which the outputs of individual subsystems inform and support the analysis of the overall propulsion system.

## 2. Research Objectives

- Develop and integrate laser-delivery and laser-sustained plasma models.
- Investigate the interaction between laser energy deposition and plasma generation.
- Develop a CAD-based thruster architecture.
- Evaluate relevant thermal-fluid and multiphysics behaviour.
- Develop propulsion-performance and system-level models.
- Investigate subsystem coupling, sensitivity, and overall system performance.
- Validate numerical predictions against published research and available experimental data.

## 3. Computational and Engineering Framework

| Subsystem | Software | Purpose |
|---|---|---|
| Laser delivery | Python | Beam propagation, optical losses and delivered power |
| Plasma generation | Python | Laser-sustained plasma modelling and energy coupling |
| Thruster design | SolidWorks | Three-dimensional geometry and component integration |
| Multiphysics analysis | ANSYS | Relevant thermal, fluid and engineering analysis |
| Performance modelling | MATLAB | Numerical analysis and propulsion performance |
| System modelling | Simulink | System-level modelling and subsystem interaction |

## 4. Repository Structure

- `CAD-SolidWorks/` — Thruster geometry, components, assemblies, engineering drawings and renders.
- `Python/` — Laser-delivery calculations, plasma models and supporting numerical modules.
- `ANSYS/` — Multiphysics simulation files and engineering analysis.
- `MATLAB/` — Propulsion-performance calculations and numerical studies.
- `Simulink/` — System-level models and subsystem integration.
- `Documentation/` — Research methodology, equations, assumptions and technical documentation.
- `Results/` — Simulation outputs, figures, comparisons and validation results.

## 5. Current Development Status

**SolidWorks CAD**
- Preliminary three-dimensional thruster design and component arrangement established.
- Further dimensional refinement and engineering development ongoing.

**Python computational modelling**
- Laser-delivery subsystem developed.
- Optical losses and delivered laser power calculations implemented.
- Integration and rebuilding of the laser-sustained plasma model is underway.

**Other engineering and computational subsystems**
- ANSYS, MATLAB and Simulink development to follow as the integrated modelling framework progresses.

## 6. Research Methodology

The project follows an iterative engineering-development approach:

1. Physical system definition and literature investigation.
2. Mathematical formulation of individual subsystems.
3. Numerical implementation and verification.
4. CAD-based physical design.
5. Engineering and multiphysics analysis.
6. Subsystem coupling and system-level modelling.
7. Parameter studies and sensitivity analysis.
8. Validation against published research and experimental evidence.

## 7. Validation and Research Integrity

Published experimental and numerical results will be used wherever appropriate to assess model accuracy.

Assumptions, provisional parameters, numerical limitations and validation status will be documented. Predicted performance will be distinguished from experimentally demonstrated performance.

## 8. Development Philosophy

The project is being developed progressively, with emphasis on physical consistency, numerical reliability, engineering feasibility and traceable research decisions.

Individual software modules are treated as interconnected subsystems of a larger propulsion investigation.

**Note:** This is an ongoing research and engineering-development project. Its architecture, models and conclusions remain subject to refinement and validation.
