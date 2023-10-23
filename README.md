# Mg(OH)2-precipitation_1DModel-Tmixer
Matlab code for the 1D solution of precipitation 
Test case for Mg(OH)2
Population Balance Model (PNM), solved using the Quadrature Method of Moments (QMOM)

@article{Raponi2023_TMixer,
author = {Raponi, Antonello and Romano, Salvatore and Battaglia, Giuseppe and Buffo, Antonio and Vanni, Marco and Cipollina, Andrea and Marchisio, Daniele},
title = {Computational Modeling of Magnesium Hydroxide Precipitation and Kinetics Parameters Identification},
journal = {Crystal Growth \& Design},
volume = {23},
number = {7},
pages = {4748-4759},
year = {2023},
doi = {10.1021/acs.cgd.2c01179},

URL = { https://doi.org/10.1021/acs.cgd.2c01179 },
eprint = { https://doi.org/10.1021/acs.cgd.2c01179 },

abstract = { Magnesium is a critical raw material and its recovery as Mg(OH)2 from saltwork brines can be realized via precipitation. 
The effective design, optimization, and scale-up of such a process require the development of a computational model accounting for the effect of fluid dynamics, 
homogeneous and heterogeneous nucleation, molecular growth, and aggregation. 
The unknown kinetics parameters are inferred and validated in this work by using experimental data produced with a T2mm-mixer and a T3mm-mixer, guaranteeing fast and efficient mixing. 
The flow field in the T-mixers is fully characterized by using the k-ε turbulence model implemented in the computational fluid dynamics (CFD) code OpenFOAM. 
The model is based on a simplified plug flow reactor model, instructed by detailed CFD simulations. It incorporates Bromley’s activity coefficient correction and a micro-mixing model for the calculation of the supersaturation ratio. 
The population balance equation is solved by exploiting the quadrature method of moments, and mass balances are used for updating the reactive ions concentrations, accounting for the precipitated solid. 
To avoid unphysical results, global constrained optimization is used for kinetics parameters identification, exploiting experimentally measured particle size distribution (PSD). 
The inferred kinetics set is validated by comparing PSDs at different operative conditions both in the T2mm-mixer and the T3mm-mixer. 
The developed computational model, including the kinetics parameters estimated for the first time in this work, will be used for the design of a prototype for the industrial precipitation of Mg(OH)2 from saltwork brines 
in an industrial environment.}}
