
# Andrade & Duggan (2020)

This repository contains code for the paper:

[Jair Andrade](https://www.linkedin.com/in/jandraor/) and [Jim
Duggan](https://ie.linkedin.com/in/jduggan). *An evaluation of
Hamiltonian Monte Carlo performance to calibrate incidence data from
disaggregated compartmental SEIR models*

All the analysis in this study can be reproduced by executing the file
**Appendix\_A.rmd**

## Abstract

Hamiltonian Monte Carlo (HMC) is a Markov chain Monte Carlo method to
estimate unknown quantities through sample generation from a target
distribution for which an analytical solution is difficult. The strength
of this method lies in its geometrical foundations, which render it
efficient for traversing high-dimensional spaces. This paper analyses
the performance of HMC in calibrating five variants of the SEIR
age-structured model. Four of these variants are related to restriction
assumptions that modellers devise to handle high-dimensional parameter
spaces. The other one corresponds to the symmetric variant. To provide a
robust answer, we compare the performance of HMC and the Nelder-Mead
algorithm (NMS), a common routine for non-linear optimisation.
Furthermore, the calibration is performed on synthetic data in order to
avoid confounding effects from errors in model selection. The results
show that HMC fits the data accurately, and provide reliable estimates
for the basic reproduction number and the age-dependent transmission
rates. Unlike NMS, the HMC performance is robust in the presence of
underreported incidences and high-dimensional complexity. This study
suggests that stringent assumptions on epidemiological can be lifted in
favour of more realistic representations. The supplementary section
presents the full set of results.
