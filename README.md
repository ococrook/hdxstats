# hdxstats
hdxstats: An R-package for statistical analysis of hydrogen deuterium exchange mass-spectrometry data. `hdxstats` provides flexible data analysis approaches for hdx-ms data. By exploiting functional data analysis techniques and emprirical Bayes methodology, the approaches within control false positives whilest maintaining power.

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)


# Basic ideas and concepts

- Start with HDX-MS data in a `QFeatures` object and perform statistical workflows and visualisation. 

`hdxstats` require data to be stored in a particular format, we provide guidance based on the output of standard software but any .csv should be coercible. The `QFeatures` object ensures reproducibility and standardisation

`hdxstats` uses rigourous statistical testing for determining the differences between hdx-ms on a peptide-per-peptide basis.

`hdxstats` provides a number of visualisation tools, please see the vignettes for examples


# Installation requirements

Users will require a working version of R, currently at least version >4.1. It is recommend to use the last version of RStudio. The package can then be installed using the `devtools` package. The package should take a few minutes to install on a regular desktop or laptop. The package will need to be loaded using `library(hdxstats)`

```{r,}
devtools::install_github("ococrook/hdxstats")
```

# Examples

See vignettes for example analysis, including epitope mapping.

# Conpatible experimental designs

We suggest discussing with a statistican prior to performing experiments. However, the following rules of thumb will be useful:

1) Aim for at least 3 time points e.g. 0, 30, 300 seconds with 1 replicate.
2) Adding more timepoints and replicates will substantially improve statistical power (ability to detect differences). 
3) Two replicates at each time point, will stablise analysis considerably.
4) For more subtle differences, more replicates will be better than more time points.
5) For better mechanistic insights, more timepoints are preferred.

# Documentation

See vignettes and package manual for documentation

# Changes

Major changes ar reported in the `NEWS` file.

# Issues and bug reports

Please report bugs by opening an issue

# feature requests

Feature requests are welcomed, please open an issue so we can dicuss your feature

# Contributions

Contributions are welcome, if you wish to contribute to the package please open an issue and we can discuss your suggest contribution.

# Funding

This development of this package was support by funding from GlaxoSmithKline, a Todd-Bird Junior Research Fellowship and the EPSRC (EP/R511742/1).

