# On community structure validation in real networks

Data and `R` code to implement the statistical methodology, and reproduce the simulations and data analysis, presented in:

Signorelli, M., Cutillo, L. (2022). On community structure validation in real networks. *Computational Statistics*, 37, 1165–1183.

You can read and download the paper (with open access) here: https://doi.org/10.1007/s00180-021-01156-6

This repository comprises 3 subfolders:

- `functions` contains functions for the computation of the one-tailed implementation of NEAT (originally implemented as a two-tailed test) and of the Community Structure Validation indices; moreover, it includes functions to simulate from the degree-corrected stochastic blockmodel for edge-valued graphs proposed in Section 4.1, as well as some further functions that are used in the simulations and data analyses;

- `simulations` contains the code used for the simulations in Sections 4.2-4.4 of the manuscript;

- `applications` contains the code used for the two example applications presented in Section 5.
