#' Estimating pathways to an estimated DAG from an exposure.
#'
#' \lifecycle{experimental}
#'
#' This algorithm estimates direct effects of a predefined exposure on each
#' network-variable for all causal models that agree with the input-network:
#' models are adjusted for all possible combinations of direct neighbors
#' (==variables in the adjacency set) -> Output is a multiset of possible
#' effects
#'
#' @inheritParams nc_outcome_estimates
#' @param .exposure Character. The exposure variable of interest.
#'
#' @return Outputs a [tibble][tibble::tibble-package] with all the models computed and their
#'   estimates.
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' metabolite_network <- simulated_data %>%
#'   select(matches("metabolite")) %>%
#'   nc_create_network()
#' simulated_data %>%
#'   nc_exposure_estimates(
#'     .graph = metabolite_network,
#'     .exposure = "exposure",
#'     .adjustment_vars = "age",
#'     .model_function = lm
#'    )
#'
nc_exposure_estimates <- function(.data, .graph, .exposure, .adjustment_vars, .model_function,
                                 .exponentiate = FALSE, ...) {
    assert_is_data.frame(.data)
    assert_is_s4(.graph)
    assert_is_a_string(.exposure)
    assert_is_character(.adjustment_vars)
    assert_is_function(.model_function)

    network_edges <- .graph@graph@edgeL

    all_possible_model_formulas <- network_edges %>%
        imap(~ stats::reformulate(c(.graph@graph@nodes[.x$edges], .adjustment_vars, .exposure),
                                  response = .y))

    # The central variable surrounded by neighbour variables in the network.
    index_node <- names(network_edges)

    # TODO: Need to consider missing values.
    all_possible_models <- all_possible_model_formulas %>%
        map(~ .model_function(
            formula = .x,
            data = .data,
            na.action = "na.fail"
            # TODO: Need to figure out how to pass other options to function
        )) %>%
        map2(index_node,
             ~ suppressMessages(MuMIn::dredge(.x, fixed = c(
                 .adjustment_vars, .exposure
             ))))

    all_top_models_tidied <- all_possible_models %>%
        # TODO: Have argument for threshold? For choosing number of models?
        map(~ MuMIn::get.models(.x, subset = TRUE)) %>%
        imap_dfr(~ .tidy_all_model_outputs(.x, .y, .exponentiate = .exponentiate)) %>%
        mutate(exposure = .exposure) %>%
        select_at(vars("exposure", everything()))

    return(all_top_models_tidied)
}

