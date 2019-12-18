
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NetCoupler

<!-- badges: start -->

[![Join the chat at
https://gitter.im/NetCoupler/Lobby](https://badges.gitter.im/NetCoupler/Lobby.svg)](https://gitter.im/NetCoupler/Lobby)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build
status](https://travis-ci.org/NetCoupler/NetCoupler.svg?branch=master)](https://travis-ci.org/NetCoupler/NetCoupler)
[![Codecov test
coverage](https://codecov.io/gh/NetCoupler/NetCoupler/branch/master/graph/badge.svg)](https://codecov.io/gh/NetCoupler/NetCoupler?branch=master)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/NetCoupler/NetCoupler?branch=master&svg=true)](https://ci.appveyor.com/project/NetCoupler/NetCoupler)
<!-- badges: end -->

The goal of NetCoupler is to estimate causal links between a set of
-omic ( e.g. metabolomics, lipidomics) or other high-dimensional
metabolic data and either a disease outcome, an exposure, or both. The
NetCoupler-algorithm, formulated by Clemens Wittenbecher and converted
into an R package by Luke Johnston, links conditional dependency
networks with an external outcome or exposure to identify direct effects
between them.

TODO: Add figure demonstrating NetCoupler

See the `vignette("netcoupler.Rmd")` for getting started.

# Installation

So far there is only the development version.

``` r
# install.packages("remotes")
remotes::install_github("NetCoupler/NetCoupler")
```

# Questions? Comments?

  - Do you have an informal question, interested in
    contributing/collaborating, or want to learn more about NetCoupler?
    Check out our [Gitter Chat
    Room](https://gitter.im/NetCoupler/Lobby).

# Contributing and Code of Conduct

Checkout the [guidelines](.github/CONTRIBUTING.md) for details on
contributing. Please note that the ‘NetCoupler’ project is released with
a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to
this project, you agree to abide by its terms.
