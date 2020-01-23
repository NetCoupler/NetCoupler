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

    model_summary <- .object %>%
        broom::glance() %>%
        mutate(model_id = model_id,
               sample_size = stats::nobs(.object))

    dplyr::left_join(model_estimates, model_summary, by = "model_id") %>%
        select_at(vars("model_id", everything()))
}
