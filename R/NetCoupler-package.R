#' @importFrom checkmate assert_data_frame assert_character assert_number
#'   assert_list assert_function assert_logical
#' @importFrom dplyr select mutate group_by if_else
#' @importFrom rlang .data
#' @importFrom purrr map map_dfr map2 imap imap_dfr
#' @import ggraph
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom tibble tibble
#' @importFrom lifecycle deprecate_soft
## usethis namespace: end
NULL

#' Simulated dataset with an underlying Directed Graph structure for the metabolites.
#'
#' @format The simulated dataset, output as a [tibble][tibble::tibble-package]
#'   contains two outcome variables (`survival_time` and `case_status`), one
#'   generic `exposure`, 12 `metabolite_*` variables, and an `age` variable.
"simulated_data"
