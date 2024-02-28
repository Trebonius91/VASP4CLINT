# VASP4CLINT
**Utility scripts and programs for VASP calculations of (liquid or solid) interface systems**

by

Julien Steffen, julien.steffen@fau.de

Andreas Mölkner, andreas.moelkner@fau.de

## General overview

This repository contains a list of scripts and programs that can be used to set up, manage and evaluate VASP
calcuations. It was initially designed to support calculations within the [CRC 1453 CLINT](https://www.sfb1452.research.fau.eu/) (catalysis at liquid
interfaces), but the scripts and programs should be general enough to use them for arbitrary VASP calculations
of surface or bulk systems.
In addition to the scripts, a manual is delivered (subfolder manual), where a general overview for the setup
of VASP calculations for interface systems is given besides the detailed explanation of all included scripts and programs.

The scripts and programs are grouped by the programming language of their implementation and alphabetically within their sections.

Currently included are:

## Python scripts:

 - **built_adsorbates.py** : Place adsorbate atoms or molecules on or in surface slabs. Positions and rotations can be controlled by input files.
 - **build_scalms.py** :  Build unit cells of liquid metal alloys on a simple cubic grid (bulk, slab and cluster possible)
 - **check_vac.py** : Checks if single atoms have departed from a surface during ML (utility for **ml_long.sh**!)
 - **eval_neb.py** : Evaluates a nudged elastic band calculation, no matter if already started or finished
 - **manage_mlff_md.py** : Starts and supervises a ML-FF trajectory of a surface slab, restarts if errors occur
 - **modify_poscar.py** : Perform several simple operations on a POSCAR file, such as shifting its atoms, multiply the cel, transform from cartesian to internal coordinates and the other way round
 - **struc_insert.py** : Insert one POSCAR (e.g. a molecule or crystal) into another POSCAR (e.g. a solvent)

## Bash scripts

 - **ml_long.sh** : Do VASP machine learning force field on-the-fly learning trajectories for arbitrary long times, even if the calculation cluster has a walltime limit.
 - **opt_long.sh** : Do a VASP geometry optimization run for arbitrary long times, even if the calculation cluster has a wallime limit

## Fortran programs:

 - **analyze_scalms** : Analyze trajectories of liquid metal surface or bulk simulations
 - **cut_unitcell** : Cut an arbitrary shaped surface slab unit cell from a given larger surface slab
 - **eval_bader** : Evaluates and visualizes Bader charge calculations
 - **mlff_select** : Selects basis functions for a VASP machine learning force field from given ML_AB files
 - **modify_xdatcar** : Modify XDATCAR trajectory files: shifts, multiplications, writing of xyz files
 - **split_freq** : Splits a large frequency calculation into arbitrary many parts and combine the results after finishing
 - **partial_dos** : Extract the partial density of states for certain atoms, elements or orbitals from a DOSCAR file
  
