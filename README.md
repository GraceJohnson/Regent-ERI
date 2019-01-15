# Two-Electron Repulsion Integrals in Regent

This code uses the Regent task-based programming language to compute two-electron repulsion integrals (ERIs)
for s functions given an input geometry and basis set, a computational bottleneck in electronic structure 
calculations. Code written by K. Grace Johnson for CS315B at Stanford University, Fall 2018. 

## Contents

This repository includes the following contents:

  * `eri.rg`: Source code for calculating ERIs
  * `eri_config.rg`: utility file for parsing input arguments
  * `geom`: folder with several input XYZ files (specifying nuclear coordinats of atoms)
  * `basis`: folder with two basis set examples 
  * `submit.sh`: example submission script (for XStream cluster)
  * `ref`: folder with correct output for h2.xyz with both basis sets and 4x4x4 with STO-3G basis set (output with -o flag)
