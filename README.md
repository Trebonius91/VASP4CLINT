# VASP4CLINT
**Utility scripts and programs for VASP calculations of (liquid or solid) interface systems**

## General overview

This repository contains a list of scripts and programs that can be used to set up, manage and evaluate VASP
calcuations. It was initially designed to support calculations within the CRC 1453 CLINT (catalysis at liquid
interfaces), but the scripts and programs should be general enough to use them for arbitrary VASP calculations
of surface or bulk systems.
In addition to the scripts, a manual is delivered (subfolder manual), where a general overview for the setup
of VASP calculations for interface systems is given besides the detailed explanation of all included scripts and programs.

The scripts and programs are divided by the programming language of their implementation.

Currently included are:

## Python scripts:

 - **modify_poscar.py** : Perform several simple operations on a POSCAR file, such as shifting its atoms, multiply the cel, transform from cartesian to internal coordinates and the other way round
 - **build_scalms.py** :  Build unit cells of liquid metal alloys on a simple cubic grid (bulk, slab and cluster possible)
 - **built_adsorbates.py** : Place adsorbate atoms or molecules on or in surface slabs. Positions and rotations can be controlled by input files.

## Bash scripts

 - **ml_long.sh** : Do VASP machine learning force field on-the-fly learning trajectories for arbitrary long times, even if the calculation cluster has a walltime limit.

## Fortran proggrams:

 - **analyze_scalms** : Analyze trajectories of liquid metal surface or bulk simulations
 - **partial_dos** : Extract the partial density of states for certain atoms, elements or orbitals from a DOSCAR file

