
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NetCoupler

<!-- badges: start -->

[![R-CMD-check](https://github.com/NetCoupler/NetCoupler/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NetCoupler/NetCoupler/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/NetCoupler/NetCoupler/branch/main/graph/badge.svg)](https://app.codecov.io/gh/NetCoupler/NetCoupler/branch/main)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![CRAN
checks](https://cranchecks.info/badges/worst/NetCoupler)](https://cran.r-project.org/web/checks/check_results_NetCoupler.html)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

The goal of NetCoupler is to estimate potential causal links between a
set of -omic (e.g. metabolomics, lipidomics) or other high-dimensional
metabolic data as a conditional dependency network and either a disease
outcome, an exposure, or both. These potential causal links are
classified as direct, ambigious, or no effects. This algorithm is
largely meant to be used with -omic style data to generate the networks
and while theoretically non-omic data could be used, we have not tested
it in that context. Given the algorithms nature, it’s primarily designed
to be used for exploration of potential mechanisms and used to
complement other analyses for a research question. It could also be used
to confirm a **pre-specified and explicit** hypothesis, similar to how
structural equation models are used. However, this might be a more niche
use.

![Overview of the NetCoupler algorithm.](man/figures/algorithm.svg)

Why or when might you want to use NetCoupler?

1.  You are interested in asking a research question on how some factor
    might influence another factor and how it might mediate through a
    metabolic network.
2.  If you want to explore how a factor might influence a metabolic
    network or how a metabolic network might influence a factor.
3.  You have an -omic dataset and want another method to explore how it
    relates to your variable of interest.

Basically, if you’re research question or objective has the general form
of:

![Type of questions or objectives that NetCoupler aims to help
answer.](man/figures/aim-question.png)

… So that you can ultimately have an answer that looks like:

![General result that NetCoupler provides that might help answer your
question.](vignettes/aim-output.png)

There are a few vignettes available in this package:

-   [Get
    Started](https://netcoupler.github.io/NetCoupler/articles/NetCoupler.html)
    (`vignette("NetCoupler")`) describes a simple overview of how and
    when to use NetCoupler, as well as a basic explanation of some of
    the components of NetCoupler.
-   [Examples with different
    models](https://netcoupler.github.io/NetCoupler/articles/examples.html)
    (`vignette("examples")`) lists different models we’ve tested that
    work with NetCoupler. If you have tried a model out that isn’t
    listed and seen success, let us know by opening an Issue or
    submitting a Pull Request (see the [contributing
    guidelines](https://github.com/NetCoupler/NetCoupler/blob/main/.github/CONTRIBUTING.md)
    for instructions on doing this).
    <!-- TODO: Add link to description vignette when its done -->

# Installation

To install the official CRAN version, use:

``` r
install.packages("NetCoupler")
```

To install the development version, use:

``` r
# install.packages("remotes")
remotes::install_github("NetCoupler/NetCoupler")
```

# Contributing and Code of Conduct

Checkout the
[guidelines](https://github.com/NetCoupler/NetCoupler/blob/main/.github/CONTRIBUTING.md)
for details on contributing. Please note that the ‘NetCoupler’ project
is released with a [Contributor Code of
Conduct](https://github.com/NetCoupler/NetCoupler/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
