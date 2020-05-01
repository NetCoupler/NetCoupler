#' Compute model estimates between an external (exposure or outcome) variable
#' and a network.
#'
#' @name nc_model_estimates
#' @param .tbl The data.frame or tibble that contains the variables of interest,
#'   including the variables passed to the network.
#' @param .graph Output graph object from `nc_create_network()`.
#' @param .exposure,.outcome Character. The exposure or outcome variable of interest.
#' @param .adjustment_vars Optional. Variables to adjust for in the models.
#' @param .model_function A function for the model to use (e.g. [stats::lm()],
#'   [stats::glm()], [survival::coxph()]). Can be any model as long as the
#'   function has the arguments `formula` and `data`.
#' @param .model_arg_list Optional. A list containing the named arguments that
#'   will be passed to the model function. A simple example would be
#'   `list(family = binomial(link = "logit"))` to specify that the `glm` model
#'   is a logistic model and not a linear one. See the examples for more on the
#'   usage.
#' @param .exponentiate Logical. Whether to exponentiate the log estimates, as
#'   computed with e.g. logistic regression models.
#' @param .external_var Argument for internal function. The variable that links
#'   to the network variables ("external" to the network).
#' @param .external_side Argument for internal function. Character vector.
#'   Either "exposure" or "outcome", to indicate which side the external
#'   variable is on relative to the network.
#'
#'   - Exposure indicating the implied directionality is from the external
#'   variable to the network variable.
#'   - Outcome indicating the implied directionality is from the network variable
#'   to the external variable.
#'
#' @description
#' \lifecycle{experimental}
#'
#' TODO: Describe more here (fill out vignette first).
#'
#' @return Outputs a [tibble][tibble::tibble-package] that contains the model
#'   estimates from either the exposure or outcome side of the network.
#'
#' @examples
#'
#' library(dplyr)
#' metabolite_network <- simulated_data %>%
#'   select(matches("metabolite")) %>%
#'   nc_create_network()
#'
#' simulated_data %>%
#'   nc_exposure_estimates(
#'     .graph = metabolite_network,
#'     .exposure = "exposure",
#'     .adjustment_vars = "age",
#'     .model_function = lm
#'    )
#'
#' simulated_data %>%
#'   nc_outcome_estimates(
#'     .graph = metabolite_network,
#'     .outcome = "survival::Surv(survival_time, case_status)",
#'     .adjustment_vars = "age",
#'     .model_function = survival::coxph,
#'     .exponentiate = TRUE
#'    )
#'
NULL

#' @describeIn nc_model_estimates Computes the model estimates for the exposure side.
#' @export
nc_exposure_estimates <-
    function(.tbl,
             .graph,
             .exposure,
             .adjustment_vars = NA,
             .model_function,
             .model_arg_list = NULL,
             .exponentiate = FALSE) {
        .compute_model_estimates(
            .tbl = .tbl,
            .graph = .graph,
            .external_var = .exposure,
            .adjustment_vars = .adjustment_vars,
            .model_function = .model_function,
            .model_arg_list = .model_arg_list,
            .exponentiate = .exponentiate,
            .external_side = "exposure"
        )
    }

#' @describeIn nc_model_estimates Computes the model estimates for the exposure side.
#' @export
nc_outcome_estimates <-
    function(.tbl,
             .graph,
             .outcome,
             .adjustment_vars = NA,
             .model_function,
             .model_arg_list = NULL,
             .exponentiate = FALSE) {
        .compute_model_estimates(
            .tbl = .tbl,
            .graph = .graph,
            .external_var = .outcome,
            .adjustment_vars = .adjustment_vars,
            .model_function = .model_function,
            .model_arg_list = .model_arg_list,
            .exponentiate = .exponentiate,
            .external_side = "outcome"
        )
}

as_edge_tbl <- function(.edge_list) {
    nodes <- names(.edge_list)
    edge_table <- purrr::map_dfr(
        .edge_list,
        .single_edge_list_to_tbl,
        .id = "source_node",
        .nodes = nodes
    )
    return(edge_table)
}

#' @describeIn nc_model_estimates Internal function. Included to document algorithm.
#' @keywords internal
.compute_model_estimates <-
    function(.tbl,
             .graph,
             .external_var,
             .adjustment_vars = NA,
             .model_function,
             .model_arg_list = NULL,
             .exponentiate = FALSE,
             .external_side = c("exposure", "outcome")) {

    # TODO: Use tidy eval style input for variables.
    assert_is_data.frame(.tbl)
    assert_is_s4(.graph)
    assert_is_a_string(.external_var)
    # TODO: This check needs to be better constructed
    if (!any(is.na(.adjustment_vars)))
        assert_is_character(.adjustment_vars)
    if (!is.null(.model_arg_list))
        assert_is_list(.model_arg_list)
    assert_is_logical(.exponentiate)
    assert_is_function(.model_function)

    # TODO: Convert this .graph to an edge table from the beginning (as an input)
    edge_table <- as_edge_tbl(.graph@graph@edgeL)

    network_combinations <- .generate_all_network_combinations(edge_table)

    .external_side <- rlang::arg_match(.external_side)
    formula_list <- .generate_formula_list(
        .network_tbl = network_combinations,
        .ext_var = .external_var,
        .ext_side = .external_side,
        .adj_vars = .adjustment_vars
    )

    variables_to_keep <- formula_list %>%
        map(all.vars) %>%
        purrr::flatten_chr() %>%
        unique()

    model_data <- .tbl %>%
        select_at(variables_to_keep) %>%
        stats::na.omit()

    model_arg_list <- list(formula = formula_list)
    other_args <- list(data = model_data)
    if (!is.null(.model_arg_list))
        other_args <- c(other_args, .model_arg_list)

    model_tbl <- model_arg_list %>%
        purrr::pmap(purrr::lift_dl(.model_function), other_args) %>%
        purrr::map2_dfr(network_combinations$index_node,
                 .tidy_models, .exponentiate = .exponentiate)

    tidied_models <- model_tbl %>%
        mutate(
            external_var = .external_var,
            adjusted_vars = dplyr::if_else(
                !any(is.na(.adjustment_vars)),
                paste(.adjustment_vars, collapse = ", "),
                NA_character_
            )
        ) %>%
        dplyr::relocate(c("external_var", "index_node", "model_id")) %>%
        dplyr::rename_all(~ gsub("\\.", "_", .))

    return(tidied_models)
}

.single_edge_list_to_tbl <- function(.edges, .nodes) {
    tibble(target_node = .nodes[.edges$edges])
}

.all_neighbour_combinations <- function(.edge_table) {
    neighbours <- .edge_table$target_node
    all_combinations <- lapply(seq_along(neighbours),
                               combn,
                               x = neighbours,
                               simplify = FALSE)
    unlist(all_combinations, recursive = FALSE)
}

.tidy_models <- function(.object, .index_node, .exponentiate) {
    model_id <- ids::random_id(1, bytes = 8)
    model_estimates <- .object %>%
        broom::tidy(exponentiate = .exponentiate, conf.int = TRUE) %>%
        # Give a unique id for the model.
        mutate(model_id = model_id,
               index_node = .index_node)

    model_estimates %>%
        .conditionally_add_model_summary(.object)
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

.generate_all_network_combinations <- function(.edge_tbl) {
    split(.edge_tbl, .edge_tbl$source_node) %>%
        map(.all_neighbour_combinations) %>%
        imap_dfr(
            ~ tibble(index_node = .y, neighbours = .x) %>%
                dplyr::add_row(index_node = .y, neighbours = NULL)
        )
}

.generate_formula_list <-
    function(.network_tbl, .ext_var, .ext_side, .adj_vars) {
        xvars_prep <-
            list(.network_tbl$neighbours) %>%
            purrr::pmap(c, .adj_vars)

        external_input <- switch(
            .ext_side,
            exposure = list(y = .network_tbl$index_node,
                            x = .ext_var),
            outcome = list(x = .network_tbl$index_node,
                           y = .ext_var)
        )

        xvar_input <-
            list(xvars_prep, external_input$x) %>%
            purrr::pmap(c) %>%
            map(stats::na.omit)

        map2(xvar_input,
             external_input$y,
             stats::reformulate)
    }

