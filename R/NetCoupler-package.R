#' @importFrom assertive.types assert_is_data.frame assert_is_s4
#'   assert_is_a_string assert_is_character assert_is_function assert_is_a_number
#'   assert_is_list assert_is_logical
#' @importFrom dplyr select mutate select_at group_by summarize_at
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
