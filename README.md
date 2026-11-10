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


## 4. Integrated System Architecture

The following diagram illustrates the intended architecture of the integrated propulsion investigation and the information flow between its computational and engineering subsystems.

The individual models are developed and verified progressively, with subsystem coupling introduced as their governing equations, inputs and outputs are established.

flowchart TD
    A["Operating Conditions<br/>Laser Power • Propellant • Geometry"] 

    subgraph PY["PYTHON — Computational Plasma and Laser Models"]
        B["Laser Delivery Model"]
        C["Optical Losses and Captured Laser Power"]
        D["Laser-Sustained Plasma (LSP) Model"]
        E["Plasma State<br/>Temperature • Density • Ionization<br/>Energy Deposition"]

        B --> C
        C -->|"Captured laser power"| D
        D --> E
    end

    subgraph CAD["SOLIDWORKS — Thruster Engineering Design"]
        F["3D Thruster Geometry"]
        G["Component Arrangement<br/>Chamber • Nozzle • Optical Access"]
        F --> G
    end

    subgraph ANSYS["ANSYS — Multiphysics Engineering Analysis"]
        H["Geometry and Computational Mesh"]
        I["Thermal / Fluid Analysis"]
        J["Engineering Field Results"]
        H --> I --> J
    end

    subgraph MATLAB["MATLAB — Propulsion Performance"]
        K["Propulsion Performance Model"]
        L["Thrust • Specific Impulse<br/>Mass Flow • Efficiency"]
        K --> L
    end

    subgraph SIM["SIMULINK — System-Level Modelling"]
        M["Integrated System Model"]
        N["Operating Conditions and<br/>System Response"]
        M --> N
    end

    A --> B
    G --> H
    E --> I
    E --> K
    J --> K
    L --> M
    N -.->|"Operating-condition updates"| A

    V["Literature and Experimental Data<br/>Verification • Validation • Model Refinement"]
    E -.-> V
    J -.-> V
    L -.-> V


## 5. Repository Structure

- `CAD-SolidWorks/` — Thruster geometry, components, assemblies, engineering drawings and renders.
- `Python/` — Laser-delivery calculations, plasma models and supporting numerical modules.
- `ANSYS/` — Multiphysics simulation files and engineering analysis.
- `MATLAB/` — Propulsion-performance calculations and numerical studies.
- `Simulink/` — System-level models and subsystem integration.
- `Documentation/` — Research methodology, equations, assumptions and technical documentation.
- `Results/` — Simulation outputs, figures, comparisons and validation results.

## 6. Current Development Status

**SolidWorks CAD**
- Preliminary three-dimensional thruster design and component arrangement established.
- Further dimensional refinement and engineering development ongoing.

**Python computational modelling**
- Laser-delivery subsystem developed.
- Optical losses and delivered laser power calculations implemented.
- Integration and rebuilding of the laser-sustained plasma model is underway.

**Other engineering and computational subsystems**
- ANSYS, MATLAB and Simulink development to follow as the integrated modelling framework progresses.

## 7. Research Methodology

The project follows an iterative engineering-development approach:

1. Physical system definition and literature investigation.
2. Mathematical formulation of individual subsystems.
3. Numerical implementation and verification.
4. CAD-based physical design.
5. Engineering and multiphysics analysis.
6. Subsystem coupling and system-level modelling.
7. Parameter studies and sensitivity analysis.
8. Validation against published research and experimental evidence.

## 8. Validation and Research Integrity

Published experimental and numerical results will be used wherever appropriate to assess model accuracy.

Assumptions, provisional parameters, numerical limitations and validation status will be documented. Predicted performance will be distinguished from experimentally demonstrated performance.

## 9. Development Philosophy

The project is being developed progressively, with emphasis on physical consistency, numerical reliability, engineering feasibility and traceable research decisions.

Individual software modules are treated as interconnected subsystems of a larger propulsion investigation.

**Note:** This is an ongoing research and engineering-development project. Its architecture, models and conclusions remain subject to refinement and validation.
