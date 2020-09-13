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
#' @keywords internal
#' @seealso See the `vignette("description")` for a
#'   detailed description of the algorithm used to classify direct effects.
#'   [nc_model_estimates]
#'
classify_effects <- function(.tbl) {
    # TODO: Use an attribute as a "check"
    assert_is_data.frame(.tbl)
    check_tbl(.tbl)

    external_variable <- identify_external_variable(.tbl)

    xvar_model_estimates <- .tbl %>%
        round_p_values() %>%
        keep_xvar_model_estimates(external_variable)

    no_neighbour_models <- xvar_model_estimates %>%
        keep_no_neighbour_models(external_variable)

    neighbour_vs_no_neighbour_models <- xvar_model_estimates %>%
        dplyr::full_join(no_neighbour_models, by = c(external_variable, "index_node"))

    models_compared <- neighbour_vs_no_neighbour_models %>%
        add_comparison_columns()

    classify_direct_effects <- models_compared %>%
        # TODO: Add p-value threshold as argument?
        add_effects_column(external_variable) %>%
        keep_relevant_data()

    return(classify_direct_effects)
}


# Helpers -----------------------------------------------------------------

keep_xvar_model_estimates <- function(.tbl, main_x_var) {
    switch(main_x_var,
           outcome = .tbl %>%
               add_neighbours_by_model("index_node") %>%
               keep_main_x_var_estimates("index_node"),
           exposure = .tbl %>%
               add_neighbours_by_model(main_x_var) %>%
               keep_main_x_var_estimates(main_x_var)
           )
}


extract_neighbour_nodes <-
    function(term_var,
             adj_var,
             main_x_var) {

    adjusted_variables <- adj_var %>%
        unique() %>%
        strsplit(", ") %>%
        unlist() %>%
        paste(collapse = "|")

    main_x_var <- unique(main_x_var)
    term_var <- term_var[!term_var %in% main_x_var]

    neighbour_nodes <- term_var[term_var != "(Intercept)"]

    neighbour_nodes[!grepl(adjusted_variables, neighbour_nodes)] %>%
        paste(collapse = ", ")
}

identify_external_variable <- function(x) {
    grep("^(exposure|outcome)$",
         names(x),
         value = TRUE)
}

# Need to round since some p-values can be really small (basically zero),
# and others can be exactly zero. So need to assume both are same.
round_p_values <- function(.tbl) {
    mutate(.tbl, dplyr::across("p_value", ~ round(., 6)))
}

add_neighbours_by_model <- function(.tbl, main_x_var) {
    .tbl %>%
        dplyr::group_by(.data$model_id) %>%
        mutate(
            neighbour_vars = extract_neighbour_nodes(
                .data$term,
                .data$adjusted_vars,
                .data[[main_x_var]]
            )
        ) %>%
        dplyr::ungroup()
}

keep_main_x_var_estimates <- function(.tbl, main_x_var) {
    .tbl %>%
        dplyr::filter(.data[[main_x_var]] == .data$term) %>%
        mutate(adj_p_value = stats::p.adjust(.data$p_value, "fdr")) %>%
        select(-all_of(c("term", "model_id")))
}

keep_no_neighbour_models <- function(.tbl, main_x_var) {
    .tbl %>%
        dplyr::filter(dplyr::across(all_of("neighbour_vars"), ~ . == "")) %>%
        dplyr::rename(
            no_neighbours_adj_p_value = .data$adj_p_value,
            no_neighbours_p_value = .data$p_value,
            no_neighbours_estimate = .data$estimate
        ) %>%
        select(all_of(
            c(
                main_x_var,
                "index_node",
                "no_neighbours_adj_p_value",
                "no_neighbours_p_value",
                "no_neighbours_estimate"
            )
        ))
}

# TODO: Use another comparator here, like AIC or something?
add_comparison_columns <- function(.tbl) {
    .tbl %>%
        mutate(
            # nnm = no neighbour models
            # nm = neighbour models
            nnm_has_bigger_pval_than_nm = .data$adj_p_value <= .data$no_neighbours_adj_p_value,
            # For no neighbour model, have variable be NA
            nnm_has_bigger_pval_than_nm = dplyr::if_else(
                .data$neighbour_vars == "",
                NA,
                .data$nnm_has_bigger_pval_than_nm
            ),
            nnm_has_same_direction_as_nm = sign(.data$estimate) ==
                sign(.data$no_neighbours_estimate),
            # For no neighbour model, have variable be NA
            nnm_has_same_direction_as_nm = dplyr::if_else(
                .data$neighbour_vars == "",
                NA,
                .data$nnm_has_same_direction_as_nm
            )
        )
}

add_effects_column <- function(.tbl, ext_var, pvalue_threshold = 0.001) {
    .tbl %>%
        dplyr::group_by(.data[[ext_var]], .data$index_node) %>%
        mutate(effect = dplyr::case_when(
            direct_effect_logic(.data$nnm_has_bigger_pval_than_nm,
                                .data$nnm_has_same_direction_as_nm,
                                .data$estimate,
                                .data$std_error) ~ "direct",
            ambigious_effect_logic(.data$nnm_has_bigger_pval_than_nm,
                                   .data$nnm_has_same_direction_as_nm,
                                   .data$adj_p_value,
                                   pvalue_threshold) ~ "ambiguous",
            TRUE ~ "none"
        )) %>%
        dplyr::ungroup()
}

direct_effect_logic <- function(bigger_pval, same_dir, est, se) {
    all(bigger_pval, na.rm = TRUE) &
        all(same_dir, na.rm = TRUE) &
        # When standard error is smaller than a quarter of the estimate. (?)
        all(se < abs(est / 4), na.rm = TRUE)
    # TODO: Include pvalue threshold?
}

ambigious_effect_logic <- function(bigger_pval, same_dir, pvalue, pvalue_threshold = 0.001) {
    any(bigger_pval, na.rm = TRUE) &
        any(same_dir, na.rm = TRUE) &
        any(pvalue <= pvalue_threshold)
}

keep_relevant_data <- function(.tbl) {
    .tbl %>%
        dplyr::filter(.data$neighbour_vars == "") %>%
        select(-matches("no_neighbour|nnm_"),
               -all_of(c(
                   "statistic",
                   "adjusted_vars",
                   "neighbour_vars"
               )))
}

check_tbl <- function(.tbl) {
    if (!all(c("index_node", "estimate", "std_error") %in% names(.tbl))) {
        rlang::abort(c("The data frame provided does not contain the proper columns. Make sure to use the data frame generated from either of these:",
                       "nc_exposure_estimates()",
                       "nc_outcome_estimates()"))
    }
}

