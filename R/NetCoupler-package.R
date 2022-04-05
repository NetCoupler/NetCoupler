#' @importFrom checkmate assert_data_frame assert_character assert_number
#'   assert_list assert_function assert_logical
#' @importFrom dplyr select mutate group_by if_else
#' @importFrom rlang .data
#' @importFrom purrr map map_dfr map2 imap imap_dfr
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom tibble tibble
#' @importFrom lifecycle deprecate_soft
## usethis namespace: end
NULL

#' Simulated dataset with an underlying Directed Graph structure for the metabolites.
#'
#' @format The simulated dataset is a [tibble][tibble::tibble-package] with the
#'   following variables:
#'
#'   - Two outcome variables (`outcome_continuous` and `outcome_binary`) along
#'   with survival time (`outcome_event_time`) that is used for the
#'   `outcome_binary` variable
#'   - A generic `exposure` variable as continuous
#'   - 12 `metabolite_*` variables
#'   - An `age` variable used as a confounder
#'
"simulated_data"
