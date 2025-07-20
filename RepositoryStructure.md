# ioc-hkl project's purpose.

Repository Overview: EPICS IOC for HKL Diffractometer Control
-------------------------------------------------------------

This is an EPICS (Experimental Physics and Industrial Control System) IOC (Input/Output Controller) designed for controlling HKL diffractometers used in X-ray and neutron scattering experiments. The repository provides a Python-based interface to the HKL library for diffractometer calculations and control.

### Key Components:

#### 1\. Core Architecture

-   EPICS Base: Uses EPICS 7.0.8 as the control system framework

-   PyDevice: Python device support for EPICS, allowing Python code to interface with EPICS records

-   HKL Library: The underlying C library (v5.0.0.3357) that performs diffractometer calculations

#### 2. Main Python Module (python/hkl.py)

The heart of the system is the hklCalculator class that provides:

-   Geometry Support: Multiple diffractometer geometries (E4CH, E4CV, K4CV, E6C, K6C)

-   Lattice Parameters: Crystal structure parameters (a, b, c, α, β, γ)

-   UB Matrix: Orientation matrix calculations

-   Reflection Management: Adding and managing crystal reflections

-   Forward/Backward Calculations: Converting between motor positions and reciprocal space coordinates

-   Pseudo-axes: Virtual axes for HKL coordinates, Q-vectors, etc.

#### 3. EPICS Database (db/hkl_main.db)

Contains thousands of EPICS records that expose the Python functionality as process variables (PVs):

-   Input/Output Records: For setting and reading values

-   Waveform Records: For arrays and complex data

-   Multi-bit Records: For geometry selection and mode control

#### 4\. Configuration Files

-   initials.template: Template for initial parameter values

-   st.cmd: IOC startup script that sets up the Python environment

-   st_base.cmd: Base configuration with default values

### Supported Diffractometer Geometries:

1.  E4CH/E4CV: Eulerian 4-circle (omega, chi, phi, 2theta)

1.  K4CV: Kappa 4-circle (komega, kappa, kphi, 2theta)

1.  E6C: Eulerian 6-circle (mu, omega, chi, phi, gamma, delta)

1.  K6C: Kappa 6-circle (mu, komega, kappa, kphi, gamma, delta)

### Key Features:

#### Crystal Orientation

-   Two-reflection method for UB matrix calculation

-   Support for multiple reflections for refinement

-   Automatic lattice parameter calculations

#### Reciprocal Space Control

-   HKL coordinate system support

-   Q-vector calculations

-   Multiple pseudo-axes for different coordinate systems

#### Motor Control

-   Forward calculations: HKL → motor positions

-   Backward calculations: motor positions → HKL

-   Multiple solution handling

-   Constraint management

#### Data Management

-   CIF file support for crystal structure data

-   Intensity calculations

-   Error handling and logging

### Usage:

1.  Installation: Follow the README.md instructions to install EPICS, HKL library, and this IOC

1.  Startup: Run ./st.cmd in the iocBoot/iocpydev directory

1.  Control: Use EPICS clients (like caget/caput) to interact with the PVs

1.  Testing: Use commands like pydev("hklApp.test()") in the EPICS shell

### Integration:

This IOC is designed to work with beamline control systems, providing a standardized interface for diffractometer control that can be integrated with other EPICS-based instrumentation and data acquisition systems.

The repository represents a complete solution for controlling diffractometers in synchrotron and neutron facilities, bridging the gap between the powerful HKL calculation library and the EPICS control system framework.

