#' Estimating pathways to an outcome from an estimated DAG.
#'
#' \lifecycle{experimental}
#'
#' This algorithm estimates direct effect of a predefined exposure
#' (network-variable) on time-to-event for all causal models that agree with the
#' input-network: Cox prop. hazards regression models are used to estimate the
#' efect of all network-variables on survival time adjusted for all possible
#' combinations of direct neighbors (adjacency set) -> Output is a multiset of
#' possible causal effects.
#'
#' @param .data A data.frame that contains the data with the metabolic variables
#'   and the outcome.
#' @param .graph The estimated graph skeleton obtained from [nc_create_network()].
#' @param .outcome Character. The outcome variable of interest.
#' @param .adjustment_vars Character vector. The variables to adjust for in the model.
#' @param .model_function A function object. The function to use for the modeling.
#' @param ... Options to pass to the model function.
#'
#' @return Outputs a [tibble::tibble()] with all the models computed and their
#'   estimates.
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(survival)
#' metabolite_network <- simulated_data %>%
#'   select(matches("metabolite")) %>%
#'   nc_create_network()
#' simulated_data %>%
#'   nc_outcome_estimates(
#'     .graph = metabolite_network,
#'     .outcome = "survival::Surv(survival_time, case_status)",
#'     .adjustment_vars = "Age",
#'     .model_function = survival::coxph
#'    )
#'
nc_outcome_estimates <- function(.data, .graph, .outcome, .adjustment_vars, .model_function, ...) {
    assert_is_data.frame(.data)
    assert_is_s4(.graph)
    assert_is_a_string(.outcome)
    assert_is_character(.adjustment_vars)
    assert_is_function(.model_function)

    network_edges <- .graph@graph@edgeL

    all_possible_model_formulas <- network_edges %>%
        imap(~ c(.graph@graph@nodes[.x$edges], .y, .adjustment_vars)) %>%
        map2(.outcome, ~ stats::reformulate(.x, response = .y))

    # TODO: Need to consider missing values.
    all_possible_models <- all_possible_model_formulas %>%
        map(~ .model_function(
            formula = .x,
            data = .data,
            na.action = "na.fail"
            # TODO: Need to figure out how to pass other options to function
        )) %>%
        map(~ suppressMessages(MuMIn::dredge(.x)))

    all_top_models_tidied <- all_possible_models %>%
        map(~ MuMIn::get.models(.x, subset = delta <= 5)) %>%
        imap_dfr(~ .tidy_all_model_outputs(.x, .y)) %>%
        mutate(outcome = .outcome)

    return(all_top_models_tidied)
}

.tidy_all_model_outputs <- function(.list, .names) {
    .list %>%
        map_dfr(.tidy_model_output) %>%
        mutate(index_node = .names) %>%
        select(index_node, everything())
}

.tidy_model_output <- function(.object) {
    .object %>%
        broom::tidy(exponentiate = TRUE) %>%
        # Give a unique id for the model.
        mutate(model_id = ids::random_id(1, bytes = 8)) %>%
        select(model_id, everything())
}
