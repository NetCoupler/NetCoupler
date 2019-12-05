#' @import stringr
#' @importFrom dplyr select mutate everything matches
#' @importFrom purrr map map_dfr map2 imap imap_dfr
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom tibble tibble
## usethis namespace: end
NULL

#' Simulated dataset with an underlying Directed Graph structure for the metabolites.
#'
#' @format The simulated dataset, output as a [tibble][tibble::tibble-package]
#'   contains two outcome variables (`survival_time` and `case_status`), one
#'   generic `exposure`, and 12 `metabolite_*` variables.
"simulated_data"
