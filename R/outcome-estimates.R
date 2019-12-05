#' Estimating pathways to an outcome from a NetCoupler DAG.
#'
#' This algorithm estimates direct effect of a predefined exposure
#' (network-variable) on time-to-event for all causal models that agree with the
#' input-network: Cox prop. hazards regression models are used to estimate the
#' efect of all network-variables on survival time adjusted for all possible
#' combinations of direct neighbors (adjacency set) -> Output is a multiset of
#' possible causal effects.
#'
#' @param .data Renamed samples x metabolites data matrix
#' @param .graph Estimated DAG skeleton of samples x metabolites data matrix
#' @param adjustment_data Exposure/phenotype data
#'
#' @return Outputs a list with model details and outcome estimates.
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
        imap_dfr(~ .tidy_all_model_outputs(.x, .y))

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
