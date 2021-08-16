
# Andrade & Duggan (2020)

This repository contains code for the
[paper](https://doi.org/10.1016/j.epidem.2020.100415):

[Jair Andrade](https://www.linkedin.com/in/jandraor/) and [Jim
Duggan](https://ie.linkedin.com/in/jduggan). *An evaluation of
Hamiltonian Monte Carlo performance to calibrate age-structured
compartmental SEIR models to incidence data*

The analysis in this study can be reproduced by executing the files:

-   **Appendix\_A.rmd**
-   **Appendix\_B.rmd**
-   **Appendix\_C.rmd**

## Abstract

Hamiltonian Monte Carlo (HMC) is a Markov chain Monte Carlo method to
estimate unknown quantities through sample generation from a target
distribution for which an analytical solution is difficult. The strength
of this method lies in its geometrical foundations, which render it
efficient for traversing high-dimensional spaces. First, this paper
analyses the performance of HMC in calibrating five variants of inputs
to an age-structured SEIR model. Four of these variants are related to
restriction assumptions that modellers devise to handle high-dimensional
parameter spaces. The other one corresponds to the unrestricted
symmetric variant. To provide a robust analysis, we compare HMC’s
performance to that of the Nelder-Mead algorithm (NMS), a common method
for non-linear optimisation. Furthermore, the calibration is performed
on synthetic data in order to avoid confounding effects from errors in
model selection. Then, we explore the variation in the method’s
performance due to changes in the scale of the problem. Finally, we fit
an SEIR model to real data. In all the experiments, the results show
that HMC approximates both the synthetic and real data accurately, and
provides reliable estimates for the basic reproduction number and the
age-dependent transmission rates. HMC’s performance is robust in the
presence of underreported incidences and high-dimensional complexity.
This study suggests that stringent assumptions on age-dependent
transmission rates can be lifted in favour of more realistic
representations. The supplementary section presents the full set of
results.
