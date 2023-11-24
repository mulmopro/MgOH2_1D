# Mg(OH)<sub>2</sub> precipitation: 1D Model for T- and Y-Mixer
MatLab code for the 1D solution of precipitation in static mixers (T- and Y-mixer)
Test case for Mg(OH)2
The Population Balance Model (PBM) is solved using the Quadrature Method of Moments (QMOM)

![T-Y](https://github.com/mulmopro/MgOH2_1D/assets/102947817/21c2a471-01ee-4b05-9983-cdf7733d51c3)

# References
If you use these codes for your research, we would be glad if you could cite our works:

## Y-Mixer
@article{Raponi2023_YMixer,
title = {Population balance modelling of magnesium hydroxide precipitation: Full validation on different reactor configurations},
journal = {Chemical Engineering Journal},
volume = {477},
pages = {146540},
year = {2023},
issn = {1385-8947},

doi = {https://doi.org/10.1016/j.cej.2023.146540},
url = {https://www.sciencedirect.com/science/article/pii/S1385894723052713},

author = {Antonello Raponi and Ramona Achermann and Salvatore Romano and Silvio Trespi and Marco Mazzotti and Andrea Cipollina and Antonio Buffo and Marco Vanni and Daniele Marchisio},
keywords = {Magnesium hydroxide, Aggregation, Predictive model, QMOM, Precipitation},
abstract = {In this paper, a predictive mono-dimensional (1D) model for Mg(OH)2 precipitation is proposed and its predictive capability is tested. 
Two different reactor configurations are analyzed and compared, namely a T-mixer and a Y-mixer followed by two consecutive diverging channels and a final coil of constant diameter. 
Both setups were chosen for their high mixing efficiency. The suspension samples were characterized by Dynamic Light Scattering (DLS), thus obtaining particle size distributions (PSD). 
The experimental data collected using the T-mixer was used to identify the kinetics parameters set, while the data obtained through the Y-mixer setup was employed to assess the model predictive capability under different fluid dynamics conditions. Computational Fluid Dynamics (CFD) simulations were conducted to characterize the flow fields and the turbulence, which were integrated into the 1D model. 
Predictions were found to be in good agreement with the experimental data and further improved after introducing a novel correction factor for the aggregation kernel.}
}

## T-Mixer
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
