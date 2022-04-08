# NetCoupler (development version)

# NetCoupler 0.1.0

* General preparation for submitting package to CRAN.

## New features

* Added to `nc_estimate_*` function output the full model list as an attribute,
that is really only necessary for those interested in the underlying models used 
for classifying the effects
* Added a continuous outcome variable to simulated data that also links in with
the DAG so that the linkage is more obvious (#97)
* Added function to create an edge table (#117)
* Incorporate tidyselect helpers into functions for selection of variables (#62)
* Added Getting Started vignette and an article on examples of using different models (#70)
* Added argument to `nc_estimate_*_links()` functions to set thresholds for
classifying links (#157)
* Added weights to be included to `as_edge_tbl()` (#142)

## Removed features

* Removed `nc_classify_effects()` and `nc_filter_estimates()`, merged them into
the two main estimation functions instead
* Model summary statistics for `lm` and `glm` models were removed for improving
computing speed (they slowed things down quite a bit)

## Internal changes

* Output all models used for classification as an attribute for the `nc_estimate_*` 
functions output
* Use lavaan instead of dagitty to generate the simulated data
* Use standard GitHub Actions and remove AppVeyor
* Refactored some code within estimation method so it runs faster
* Tidied up the unit tests to run faster
* Removed duplicate or extra roxygen examples and instead referenced a common 
source with `@seealso`
* Removed survival dependency
* Switch to using main instead of master branch

# NetCoupler 0.0.4

## Added features

* For `lm` and `glm` models, model summary statistics are added (#88).
* Add a function to classify the direct effects between outcome or exposure and 
the network (#98).
* Add function to plot network graph: `nc_plot_network()` (#89, #110).
* Added helper functions `nc_adjacency_graph()`,
`nc_adjacency_matrix()`, and `nc_partial_corr_matrix()` 
to help create the weights for the network plot.
(Issue #80, PR #89).
* Removed soft deprecated functions. Using MuMIn over glmulti doesn't change the
results too much, see #60 for details (#83).
* Removed stringr dependency (#65, #83).

## Fixed bugs and other problems

* Fix bug where too many digits caused a problem for `pcor()` (#125, #131).
* Fix bug that didn't properly filter variables nor identify neighbour nodes
in `nc_filter_estimates()` (#109).
* Fix problem with `nc_standardize()` that prevented the ability to use the `.regressed_on`.
argument to extract residuals (#108).
* Input dataset can include missingness. Input data is treated as complete case
for only the variables used in the modelling (#88).

## Internal changes

* Rewrote underlying model estimation algorithm so it doesn't use MuMIn and
so there is one unified function for both outcome and exposure side estimation
(#101)

# NetCoupler 0.0.3.9000

* Add `nc_standardize()` function to standardize the metabolic variables (#73).
* Export tidyselect functions like `matches()` or `starts_with()` (#73).
* Add CONTRIBUTING guidelines (#56).
* Add lifecycle badges to functions, soft deprecating `net_coupler_out()`,
`getExp.coef.permetabolite()`, and `getExp.coef.out()` (#59)
* Add defensive checks to input arguments with assertive.types (#59).
* Add AppVeyor to repo. Started Travis to run on repo (#61).
* Added function for exposure side estimation: `nc_exposure_estimates()`

# NetCoupler 0.0.2.9000 

* Major revision of underlying code for generating the outcome-network link estimation (#55),
resulting in created and streamlined `nc_outcome_estimates()` function. Because of
this streamlining, the code is much faster and with the move to use MuMIn we can 
remove our dependency on rJava via glmulti.
* Tidied up `nc_create_network()` function so that only the graph skeleton is output (#55).
* Started cleaning up, along with leftover files.
* Updated and generated documentation of `nc_create_network()`.
* Added unit tests for `nc_create_network()` and the outcome estimation functions.
Travis and code coverage were added as well.
* Renamed `nc_make_network()` to `nc_create_network()` and moved into own file.
* Modularized `nc_make_network()` code and moved into another file.

# NetCoupler 0.0.1.9000

* Added a `NEWS.md` file to track changes to the package.
* Added package infrastructure
* Added initial code that other projects use
* Added basic introduction vignette
