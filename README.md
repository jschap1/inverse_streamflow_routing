<!--
When you have set up your repo you may need to change those badges
-->

[![MISS_HIT](https://github.com/Remi-Gau/template_matlab_analysis/actions/workflows/miss_hit.yml/badge.svg?branch=main)](https://github.com/Remi-Gau/template_matlab_analysis/actions/workflows/miss_hit.yml)
[![MATLAB: test and coverage](https://github.com/Remi-Gau/template_matlab_analysis/actions/workflows/matlab_test_and_coverage.yaml/badge.svg)](https://github.com/Remi-Gau/template_matlab_analysis/actions/workflows/matlab_test_and_coverage.yaml)
[![Octave: test and coverage](https://github.com/Remi-Gau/template_matlab_analysis/actions/workflows/octave_test_and_coverage.yml/badge.svg?branch=main)](https://github.com/Remi-Gau/template_matlab_analysis/actions/workflows/octave_test_and_coverage.yml)
[![codecov](https://codecov.io/gh/Remi-Gau/template_matlab_analysis/branch/master/graph/badge.svg?token=aFXb7WSAsm)](https://codecov.io/gh/Remi-Gau/template_matlab_analysis)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Documentation Status: latest](https://readthedocs.org/projects/template_matlab_analysis/badge/?version=latest)](https://template_matlab_analysis.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Remi-Gau/template_matlab_analysis/main)

- [Template repository for MATLAB / Octave projects](#template-repository-for-matlab--octave-projects)
   - [How to install and use this template](#how-to-install-and-use-this-template)
      - [Install with Github](#install-with-github)
      - [Install with cookiecutter](#install-with-cookiecutter)
   - [Configuration](#configuration)

This repository was created using the cookiecutter template. There may be some unused files and folders left over from the template.

# Inverse streamflow routing

Inverse streamflow routing (ISR) uses a flow direction map and time series of discharge measurements at points along the river to estimate runoff throughout the river basin.

Requires GDAL and the R packages gstat and CoSMoS.

## Contents

Scripts: Workflows for ISR
Functions: Main and secondary functions for performing ISR and evaluating the results.

## References
* Pan, M., & Wood, E. F. (2013). Inverse streamflow routing. Hydrology and Earth System Sciences, 17(11), 4577–4588. https://doi.org/10.5194/hess-17-4577-2013

* Fisher, C. K., Pan, M., & Wood, E. F. (2020). Spatiotemporal assimilation-interpolation of discharge records through inverse streamflow routing. Hydrology and Earth System Sciences, 24(1), 293–305. https://doi.org/10.5194/hess-24-293-2020

* Yang, Y., Lin, P., Fisher, C. K., Turmon, M., Hobbs, J., Emery, C. M., … Pan, M. (2019). Enhancing SWOT discharge assimilation through spatiotemporal correlations. Remote Sensing of Environment, 234(September), 111450. https://doi.org/10.1016/j.rse.2019.111450

Fisher et al., 2021, HESS
