#' Classify direct, ambiguous, or no effect between exposure and network nodes.
#'
#' @description
#' \lifecycle{experimental}
#'
#' Using the output of the [nc_exposure_estimates()] or [nc_outcome_estimates()],
#' classify whether the exposure variable has a direct, ambiguous, or no link
#' with the index node metabolic variable, conditional on neighbouring metablic
#' variables and on potential confounders.
#'
#' @param .tbl Multi-model estimates generated from [nc_exposure_estimates()] or
#'   [nc_outcome_estimates()] that contain the model summaries.
#'
#' @return Outputs a [tibble][tibble::tibble-package] with model estimates
#'   between the exposure and the individual index network nodes, along with the
#'   classification of direct effect or not.
#' @export
#' @seealso See the `vignette("description")` for a
#'   detailed description of the algorithm used to classify direct effects.
#'
#' @examples
#'
#' library(dplyr)
#' metabolite_network <- simulated_data %>%
#'   select(matches("metabolite")) %>%
#'   nc_create_network()
#' multimodel_exposure <- simulated_data %>%
#'   mutate(Random = rnorm(nrow(.))) %>%
#'   nc_exposure_estimates(
#'     .graph = metabolite_network,
#'     .exposure = "exposure",
#'     .adjustment_vars = c("age", "Random"),
#'     .model_function = lm
#'    )
#'
#' nc_filter_estimates(multimodel_exposure)
#' nc_classify_effects(multimodel_exposure)
#'
#' library(survival)
#' multimodel_outcome <- simulated_data %>%
#'   mutate(Random = rnorm(nrow(.))) %>%
#'   nc_outcome_estimates(
#'     .graph = metabolite_network,
#'     .outcome = "survival::Surv(survival_time, case_status)",
#'     .adjustment_vars = c("age", "Random"),
#'     .model_function = survival::coxph,
#'     .exponentiate = TRUE
#'    )
#'
#' nc_filter_estimates(multimodel_outcome)
#' nc_classify_effects(multimodel_outcome)
#'
nc_classify_effects <- function(.tbl) {
    # TODO: Use an attribute as a "check"
    assert_is_data.frame(.tbl)

    if (!all(c("index_node", "estimate", "std_error") %in% names(.tbl))) {
        rlang::abort(c("The data frame provided does not contain the proper columns. Make sure to use the data frame generated from either of these:",
                       "nc_exposure_estimates()",
                       "nc_outcome_estimates()"))
    }

    # Need to round since some p-values can be really small (basically zero),
    # and others can be exactly zero. So need to assume both are same.
    rounded_p_values <- .tbl %>%
        dplyr::mutate_at("p_value", ~ round(., 6))

    filtered_estimates <- rounded_p_values %>%
        nc_filter_estimates()

    external_variable <- filtered_estimates %>%
        select(matches("^(exposure|outcome)$")) %>%
        names()

    no_neighbours <- filtered_estimates %>%
        dplyr::filter_at("neighbour_vars", ~ . == "") %>%
        mutate(no_neighbours_adj_p_value = stats::p.adjust(.data$p_value, "fdr")) %>%
        dplyr::rename(no_neighbours_p_value = .data$p_value,
                      no_neighbours_estimate = .data$estimate) %>%
        select_at(c(
            external_variable,
            "index_node",
            "no_neighbours_adj_p_value",
            "no_neighbours_p_value",
            "no_neighbours_estimate"
        ))

    classify_direct_effects <- filtered_estimates %>%
        mutate(adj_p_value = stats::p.adjust(.data$p_value, "fdr")) %>%
        dplyr::full_join(no_neighbours, by = c(external_variable, "index_node")) %>%
        mutate(
            smaller_pvalue = .data$adj_p_value <= .data$no_neighbours_adj_p_value,
            same_direction = sign(.data$estimate) == sign(.data$no_neighbours_estimate),
            smaller_pvalue = dplyr::if_else(.data$neighbour_vars == "", NA, .data$smaller_pvalue),
            same_direction = dplyr::if_else(.data$neighbour_vars == "", NA, .data$same_direction)
        ) %>%
        dplyr::group_by_at(c(external_variable, "index_node")) %>%
        mutate(direct_effect = dplyr::case_when(
            all(.data$smaller_pvalue, na.rm = TRUE) &
                all(.data$same_direction, na.rm = TRUE) &
                # When standard error is smaller than the estimate.
                all(.data$std_error < abs(.data$estimate), na.rm = TRUE) ~ "direct",
            any(.data$smaller_pvalue, na.rm = TRUE) &
                any(.data$same_direction, na.rm = TRUE) ~ "ambiguous",
            TRUE ~ "none"
        )) %>%
        dplyr::filter(.data$neighbour_vars == "") %>%
        select(-dplyr::matches("no_neighbour"),
               -dplyr::all_of(
                   c(
                       "statistic",
                       "adjusted_vars",
                       "neighbour_vars",
                       "smaller_pvalue",
                       "same_direction"
                   )
               ))

    return(classify_direct_effects)
}

#' @describeIn nc_classify_effects Filter out estimates so only exposure
#'   estimates are kept. Adds a column with the neighbouring metabolite
#'   variables listed.
#' @export
nc_filter_estimates <- function(.tbl) {
    .filter_by <- .tbl %>%
        select(matches("^(exposure|outcome)$")) %>%
        names()
    if (.filter_by == "outcome") {
        .tbl %>%
            dplyr::group_by_at("model_id") %>%
            mutate(neighbour_vars = .extract_neighbour_nodes(.data$term,
                                                             .data$adjusted_vars,
                                                             .data$index_node)) %>%
            dplyr::ungroup() %>%
            dplyr::filter(.data$index_node == .data$term) %>%
            select(-dplyr::all_of(c("term", "model_id")))
    } else {
        .tbl %>%
            dplyr::group_by_at("model_id") %>%
            mutate(neighbour_vars = .extract_neighbour_nodes(.data$term, .data$adjusted_vars)) %>%
            dplyr::ungroup() %>%
            dplyr::filter(.data$term == .filter_by) %>%
            select(-dplyr::all_of(c("term", "model_id")))

    }
}

.extract_neighbour_nodes <- function(.term_var, .adj_var, .index_node = NULL) {
    adjusted_variables <- .adj_var %>%
        unique() %>%
        strsplit(", ") %>%
        unlist()

    if (!is.null(.index_node)) {
        .index_node <- unique(.index_node)
        .term_var <- .term_var[!.term_var %in% .index_node]
    }

    neighbour_nodes <- .term_var[!.term_var %in% c("exposure", "outcome", "(Intercept)")]

    neighbour_nodes[!neighbour_nodes %in% adjusted_variables] %>%
        paste(collapse = ", ")
}
