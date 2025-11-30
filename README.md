# SECAR Detection – Geant4 Simulation

This repository contains a **Geant4-based simulation of the detector systems used at SECAR (FRIB/NSCL)** for studies of charged-particle and neutron detection near the target position.

⚠️ **Note:**  
This simulation **does not model the SECAR separator optics**. It focuses exclusively on the **target region and downstream detector response**.

---

## Features

The simulation includes modular components for:
Configuration selection can be determined at **include/DetectorConstruction.hh**

- **Primary Source** 
  - Primary beam generation at the target
  - **Radioactive source simulation at the target location (optional)**

- **Target** *(optional)*
  - Gas targets (He/H)
  - Solid targets

- **Stripper Foil** *(optional)*
  - Thin downstream foil for charge-state manipulation studies

- **JENSA Gas-Jet Chamber** *(optional)*
  - Full geometry and material description

- **Neutron Detectors** *(optional)*
  - Liquid scintillator detector array

- **Charged-Particle Detectors** *(optional)*
  - **Ion Chamber (IC)**
  - **Double-Sided Silicon Strip Detector (DSSD)**
  - Positioned immediately downstream of the target

- **Monitor Detectors** *(optional)*
  - **PIPs detectors inside the JENSA chamber** *(under development)*

All detector and geometry components can be enabled or disabled independently to model different experimental configurations.

---

## Scope

This code is intended for:

- Detector response simulations
- Geometry and layout optimization
- Efficiency and resolution studies
- Background estimation near the target
- Acceptance modeling for target-proximate detectors

It is **not intended for beam transport or separator modeling**.

---

## Requirements

- **Geant4 version 11.3** (tested and currently supported)

---

## Outputs

The simulation produces standard Geant4 event data, including detector hits and energy deposits.
Generated products typically include:

- ROOT files

These outputs are **locally generated and excluded from version control**.

---

## Notes

Large generated outputs, CAD models, and detector response files (e.g., `.root`, `.csv`, `.stl`) are intentionally **not tracked** in this repository.

---

## Maintainer

Developed and maintained for SECAR detector studies at FRIB/NSCL.

