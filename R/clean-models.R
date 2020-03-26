.tidy_all_model_outputs <- function(.list, .names, .exponentiate) {
    .list %>%
        map_dfr(.tidy_model_output, .exponentiate = .exponentiate) %>%
        mutate(index_node = .names) %>%
        select_at(vars("index_node", everything()))
}

.tidy_model_output <- function(.object, .exponentiate) {
    model_id <- ids::random_id(1, bytes = 8)
    model_estimates <- .object %>%
        broom::tidy(exponentiate = .exponentiate, conf.int = TRUE) %>%
        # Give a unique id for the model.
        mutate(model_id = model_id)

    model_estimates <- model_estimates %>%
        .conditionally_add_model_summary(.object)

    model_estimates %>%
        select_at(vars("model_id", everything()))
}

.conditionally_add_model_summary <- function(.tbl, .object) {
    if (!any(class(.object) %in% c("lm", "glm"))) {
        return(.tbl)
    }

    if (requireNamespace("broom", quietly = TRUE)) {
        model_id <- unique(.tbl$model_id)

        model_summary <- .object %>%
            broom::glance() %>%
            mutate(model_id = model_id,
                   sample_size = stats::nobs(.object)) %>%
            select(any_of(c("model_id", "r_squared", "adj_r_squared",
                            "df", "logLik", "AIC", "BIC", "sample_size")))

        summary_with_estimates <- .tbl %>%
            dplyr::left_join(model_summary, by = "model_id")
        return(summary_with_estimates)
    } else {
        rlang::inform("If you'd like model summary statistics like AIC added to the results, please install the broom package.")
        return(.tbl)
    }
}
