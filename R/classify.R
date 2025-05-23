#' @title
#' Classify direct, ambiguous, or no effect between exposure and network nodes.
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' Using the output of the [nc_estimate_exposure_links()] or [nc_estimate_outcome_links()],
#' classify whether the exposure variable has a direct, ambiguous, or no link
#' with the index node metabolic variable, conditional on neighbouring metablic
#' variables and on potential confounders.
#'
#' @param data Multi-model estimates generated from [nc_estimate_exposure_links()] or
#'   [nc_estimate_outcome_links()] that contain the model summaries.
#' @inheritParams nc_estimate_links
#'
#' @return Outputs a [tibble][tibble::tibble-package] with model estimates
#'   between the exposure and the individual index network nodes, along with the
#'   classification of direct effect or not.
#' @keywords internal
#' @noRd
#' @seealso See the `vignette("description")` for a
#'   detailed description of the algorithm used to classify direct effects.
#'   See [nc_estimate_links] for examples on using NetCoupler.
#'
classify_effects <- function(data,
                             classify_option_list = classify_options()) {
    assert_data_frame(data)
    check_tbl(data)

    external_variable <- identify_external_variable(data)

    data <- data %>%
        dplyr::mutate(network_neighbours = purrr::map_chr(.data$network_neighbours,
                                                          paste0, collapse = ","))

    xvar_model_estimates <- data %>%
        round_p_values() %>%
        keep_xvar_model_estimates(external_variable)

    no_neighbour_models <- xvar_model_estimates %>%
        keep_no_neighbour_models(external_variable)

    neighbour_vs_no_neighbour_models <- xvar_model_estimates %>%
        dplyr::full_join(no_neighbour_models, by = c(external_variable, "index_node"))

    models_compared <- neighbour_vs_no_neighbour_models %>%
        add_comparison_columns()

    classify_direct_effects <- models_compared %>%
        add_effects_column(external_variable, classify_option_list) %>%
        keep_relevant_data()

    return(classify_direct_effects)
}

# Helpers -----------------------------------------------------------------

keep_xvar_model_estimates <- function(data, main_x_var) {
    switch(main_x_var,
           outcome = data %>%
               keep_main_x_var_estimates("index_node"),
           exposure = data %>%
               keep_main_x_var_estimates(main_x_var)
           )
}

identify_external_variable <- function(x) {
    grep("^(exposure|outcome)$",
         names(x),
         value = TRUE)
}

# Need to round since some p-values can be really small (basically zero),
# and others can be exactly zero. So need to assume both are same.
# TODO: Or do log of the p-value ie. like in GWAS studies.
round_p_values <- function(data) {
    mutate(data, dplyr::across("p_value", ~ round(., 6)))
}

keep_main_x_var_estimates <- function(data, main_x_var) {
    data %>%
        dplyr::filter(.data[[main_x_var]] == .data$term) %>%
        mutate(fdr_p_value = stats::p.adjust(.data$p_value, "fdr")) %>%
        select(-all_of(c("term", "model_id", "p_value")))
}

keep_no_neighbour_models <- function(data, main_x_var) {
    data %>%
        dplyr::filter(dplyr::if_all(all_of("network_neighbours"), ~ . == "")) %>%
        dplyr::rename(
            no_neighbours_fdr_p_value = "fdr_p_value",
            no_neighbours_estimate = "estimate"
        ) %>%
        select(all_of(
            c(
                main_x_var,
                "index_node",
                "no_neighbours_fdr_p_value",
                "no_neighbours_estimate"
            )
        ))
}

# TODO: Use another comparator here, like AIC or something?
add_comparison_columns <- function(data) {
    data %>%
        mutate(
            # nnm = no neighbour models
            # nm = neighbour models
            nnm_has_same_direction_as_nm = sign(.data$estimate) ==
                sign(.data$no_neighbours_estimate),
            # For no neighbour model, have variable be NA
            nnm_has_same_direction_as_nm = dplyr::if_else(
                .data$network_neighbours == "",
                NA,
                .data$nnm_has_same_direction_as_nm
            )
        )
}

add_effects_column <- function(data, ext_var, classify_option_list) {
    data %>%
        dplyr::group_by(.data[[ext_var]], .data$index_node) %>%
        mutate(effect = dplyr::case_when(
            direct_effect_logic(.data$no_neighbours_fdr_p_value,
                                .data$nnm_has_same_direction_as_nm,
                                .data$fdr_p_value,
                                classify_option_list = classify_option_list) ~ "direct",
            ambigious_effect_logic(.data$no_neighbours_fdr_p_value,
                                   .data$fdr_p_value,
                                   classify_option_list = classify_option_list) ~ "ambiguous",
            TRUE ~ "none"
        )) %>%
        dplyr::ungroup()
}

direct_effect_logic <- function(nnm_pvalue, same_dir, pvalue, classify_option_list) {
    # TODO: Add this later.
    # implementation <- rlang::arg_match(implementation)
    implementation <- "original"
    switch(
        implementation,
        # the criteria for direct effect identification were i) P < 0.05 in the
        # non-neighbor-adjusted model (rule1), ii) P < 0.1 in all neighbor-adjusted
        # models (rule 2), iii) same sign of the beta in all the models (rule 3).
        original =
            all(nnm_pvalue < classify_option_list$single_metabolite_threshold, na.rm = TRUE) &
            all(same_dir, na.rm = TRUE) &
            all(pvalue < classify_option_list$network_threshold, na.rm = TRUE)
        # TODO: Add other methods later.
    )
}

# TODO: Merge with direct above? depends on implementation
ambigious_effect_logic <- function(nnm_pvalue, pvalue, classify_option_list) {
    any(nnm_pvalue < classify_option_list$single_metabolite_threshold, na.rm = TRUE) |
        any(pvalue < classify_option_list$network_threshold, na.rm = TRUE)
}

keep_relevant_data <- function(data) {
    data %>%
        dplyr::filter(.data$network_neighbours == "") %>%
        select(-matches("no_neighbour|nnm_"),
               -all_of(c(
                   "statistic",
                   "adjusted_vars",
                   "network_neighbours"
               )))
}

check_tbl <- function(data) {
    if (!all(c("index_node", "estimate", "std_error") %in% names(data))) {
        rlang::abort(c("The data frame provided does not contain the proper columns. Make sure to use the data frame generated from either of these:",
                       "nc_estimate_exposure_links()",
                       "nc_estimate_outcome_links()"))
    }
}

