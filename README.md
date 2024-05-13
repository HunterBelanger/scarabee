# Scarabée

Scarabée is a simplified lattice physics code for light water reactor (LWR)
neutronics calculations. It currently has the following features:

* Resonance self-shielding according to Carlvik's 2-term rational approximation
* 1D Collision Probabilities annular pin cell solver for k-eigenvalue problems
* 2D Method of Characteristics solver for fixed-source and k-eigenvalue problems

This project closely follows the methods outlined in *Methods of Steady-State
Reactor Physics in Nuclear Design* by Stamm'ler and Abbate, and *Lattice
Physics Computations* by Knott and Yamamoto, from the *Handbook of Nuclear
Engineering*, edited by D. Cacuci.

Scarabée uses a custom HDF5 formated nuclear data library which is easy and
intuitive to understand. The `data` directory contains a helper library and
scripts to generate a multi-group nuclear data library with the FRENDY nuclear
data processing code. FRENDY is free/open-source software, and can be downloaded
[here](https://rpg.jaea.go.jp/main/en/program_frendy/). Currently, only FRENDY
is supported, due to its advanced multi-group processing features.
