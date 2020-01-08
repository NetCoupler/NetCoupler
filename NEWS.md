# NetCoupler 0.0.3.9000 (development version)

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
