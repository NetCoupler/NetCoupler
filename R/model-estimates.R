#' Compute model estimates between an external (exposure or outcome) variable
#' and a network.
#'
#' @name nc_estimate_links
#' @param data The data.frame or tibble that contains the variables of interest,
#'   including the variables used to make the network.
#' @param edge_tbl Output graph object from `nc_estimate_network()`, converted
#'   to an edge table using `as_edge_tbl()`.
#' @param exposure,outcome Character. The exposure or outcome variable of interest.
#' @param adjustment_vars Optional. Variables to adjust for in the models.
#' @param model_function A function for the model to use (e.g. [stats::lm()],
#'   [stats::glm()], survival::coxph()). Can be any model as long as the
#'   function has the arguments `formula` and `data`.
#' @param model_arg_list Optional. A list containing the named arguments that
#'   will be passed to the model function. A simple example would be
#'   `list(family = binomial(link = "logit"))` to specify that the `glm` model
#'   is a logistic model and not a linear one. See the examples for more on the
#'   usage.
#' @param exponentiate Logical. Whether to exponentiate the log estimates, as
#'   computed with e.g. logistic regression models.
#' @param external_var Argument for internal function. The variable that links
#'   to the network variables ("external" to the network).
#' @param external_side Argument for internal function. Character vector.
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
#' - `nc_estimate_exposure_links()`: Computes the model estimates for the exposure side.
#' - `nc_estimate_outcome_links()`: Computes the model estimates for the exposure side.
#'
#' @return Outputs a [tibble][tibble::tibble-package] that contains the model
#'   estimates from either the exposure or outcome side of the network.
#'
#' @seealso The vignette ... has more details on how to use NetCoupler with
#'   different models.
#'
#' @examples
#'
#' \dontrun{
#' library(dplyr)
#' standardized_data <- simulated_data %>%
#'     nc_standardize(starts_with("metabolite"))
#'
#' metabolite_network <- simulated_data %>%
#'     nc_standardize(starts_with("metabolite"),
#'                    regressed_on = "age") %>%
#'     nc_estimate_network(starts_with("metabolite"))
#' edge_table <- as_edge_tbl(metabolite_network)
#'
#' standardized_data %>%
#'   nc_estimate_exposure_links(
#'     edge_tbl = edge_table,
#'     exposure = "exposure",
#'     model_function = lm
#'    )
#'
#' # Using a logistic regression
#' standardized_data %>%
#'   nc_estimate_outcome_links(
#'     edge_tbl = edge_table,
#'     outcome = "case_status",
#'     model_function = glm,
#'     adjustment_vars = "age",
#'     model_arg_list = list(family = binomial(link = "logit")),
#'     exponentiate = TRUE
#'   )
#'
#' # Adding random confounders to adjust for
#' standardized_data %>%
#'   mutate(Random = rnorm(nrow(.)),
#'          Sex = sample(rep(c("F", "M"), times = nrow(.) / 2))) %>%
#'   nc_estimate_exposure_links(
#'     edge_tbl = edge_table,
#'     exposure = "exposure",
#'     adjustment_vars = c("age", "Random", "Sex"),
#'     model_function = lm
#'    )
#' }
#'
NULL

#' @rdname nc_estimate_links
#' @export
nc_estimate_exposure_links <-
    function(data,
             edge_tbl,
             exposure,
             adjustment_vars = NA,
             model_function,
             model_arg_list = NULL,
             exponentiate = FALSE,
             implementation = c("updated", "original")) {
        multiple_models <- compute_model_estimates(
            data = data,
            edge_tbl = edge_tbl,
            external_var = exposure,
            adjustment_vars = adjustment_vars,
            model_function = model_function,
            model_arg_list = model_arg_list,
            exponentiate = exponentiate,
            external_side = "exposure"
        ) %>%
            dplyr::rename("exposure" = "external_var")
        list(
            full_model_list = multiple_models,
            classified_effects = classify_effects(multiple_models, implementation = implementation)
        )
    }

#' @rdname nc_estimate_links
#' @export
nc_estimate_outcome_links <-
    function(data,
             edge_tbl,
             outcome,
             adjustment_vars = NA,
             model_function,
             model_arg_list = NULL,
             exponentiate = FALSE,
             implementation = c("updated", "original")) {
        multiple_models <- compute_model_estimates(
            data = data,
            edge_tbl = edge_tbl,
            external_var = outcome,
            adjustment_vars = adjustment_vars,
            model_function = model_function,
            model_arg_list = model_arg_list,
            exponentiate = exponentiate,
            external_side = "outcome"
        ) %>%
            dplyr::rename("outcome" = "external_var")
        list(
            full_model_list = multiple_models,
            classified_effects = classify_effects(multiple_models, implementation = implementation)
        )
}

compute_model_estimates <-
    function(data,
             edge_tbl,
             external_var,
             adjustment_vars = NA,
             model_function,
             model_arg_list = NULL,
             exponentiate = FALSE,
             external_side = c("exposure", "outcome")) {

    # TODO: Use tidy eval style input for variables.
    assert_data_frame(data)
    assert_data_frame(edge_tbl)
    assert_character(external_var)
    # TODO: This check needs to be better constructed
    if (!any(is.na(adjustment_vars)))
        assert_character(adjustment_vars)
    if (!is.null(model_arg_list))
        assert_list(model_arg_list)
    assert_logical(exponentiate)
    assert_function(model_function)

    network_combinations <- generate_all_network_combinations(edge_tbl)

    external_side <- rlang::arg_match(external_side)
    # TODO: Change to not create formulas as environments (they take up a lot of space)
    formula_df <- generate_formula_df(
        network_object_tbl = network_combinations,
        ext_var = external_var,
        ext_side = external_side,
        adj_vars = adjustment_vars
    )

    variables_to_keep <- na.omit(unique(c(
        network_combinations$index_node,
        purrr::flatten_chr(network_combinations$neighbours),
        external_var,
        adjustment_vars
    )))

    model_data <- data %>%
        select(all_of(variables_to_keep)) %>%
        stats::na.omit()

    other_args <- list(data = model_data)
    if (!is.null(model_arg_list))
        other_args <- c(other_args, model_arg_list)

    model_tbl <- formula_df %>%
        purrr::pmap_dfr(
            run_model_and_tidy,
            model_function = model_function,
            model_args = other_args,
            exponentiate = exponentiate
        )

    tidied_models <- model_tbl %>%
        mutate(
            external_var = external_var,
            adjusted_vars = dplyr::if_else(
                !any(is.na(adjustment_vars)),
                paste(adjustment_vars, collapse = ", "),
                NA_character_
            )
        ) %>%
        dplyr::relocate(c("external_var", "index_node", "model_id")) %>%
        dplyr::rename_with(~ gsub("\\.", "_", .))

    return(tidied_models)
}

# Helpers -----------------------------------------------------------------

run_model_and_tidy <- function(x, y, index_node, model_function, model_args, exponentiate) {
    formula <- stats::reformulate(termlabels = x, response = y)
    function_with_args_as_list <- purrr::lift_dl(model_function)
    model_args <- c(list(formula = formula),
                    model_args)
    function_with_args_as_list(model_args) %>%
        tidy_models(index_node = index_node, exponentiate = exponentiate)
}

all_neighbour_combinations <- function(edge_tbl) {
    neighbours <- edge_tbl$target_node
    all_combinations <- lapply(seq_along(neighbours),
                               utils::combn,
                               x = neighbours,
                               simplify = FALSE)
    unlist(all_combinations, recursive = FALSE)
}

tidy_models <- function(model_object, index_node, exponentiate) {
    model_id <- ids::random_id(1, bytes = 8)
    model_estimates <- model_object %>%
        broom::tidy(exponentiate = exponentiate, conf.int = FALSE) %>%
        # Give a unique id for the model.
        mutate(model_id = model_id,
               index_node = index_node)

    model_estimates
    # For now, this code increases computing time for not much benefit (that I see)
    # model_estimates %>%
    #     conditionally_add_model_summary(model_object)
}

conditionally_add_model_summary <- function(data, model_object) {
    if (!any(class(model_object) %in% c("lm", "glm"))) {
        return(data)
    }

    if (requireNamespace("broom", quietly = TRUE)) {
        model_id <- unique(data$model_id)

        model_summary <- model_object %>%
            broom::glance() %>%
            mutate(model_id = model_id,
                   sample_size = stats::nobs(model_object)) %>%
            select(any_of(c("model_id", "r_squared", "adj_r_squared",
                            "df", "logLik", "AIC", "BIC", "sample_size")))

        summary_with_estimates <- data %>%
            dplyr::left_join(model_summary, by = "model_id")
        return(summary_with_estimates)
    } else {
        rlang::inform("If you'd like model summary statistics like AIC added to the results, please install the broom package.")
        return(data)
    }
}

generate_all_network_combinations <- function(edge_tbl) {
    split(edge_tbl, edge_tbl$source_node) %>%
        map(all_neighbour_combinations) %>%
        imap_dfr(
            ~ tibble(index_node = .y, neighbours = .x) %>%
                dplyr::add_row(index_node = .y, neighbours = NULL)
        )
}

generate_formula_df <-
    function(network_object_tbl, ext_var, ext_side, adj_vars) {
        xvars_prep <-
            list(network_object_tbl$neighbours) %>%
            purrr::pmap(c, adj_vars) %>%
            map(sort)

        external_input <- switch(
            ext_side,
            exposure = list(y = network_object_tbl$index_node,
                            x = ext_var),
            outcome = list(x = network_object_tbl$index_node,
                           y = ext_var)
        )

        xvar_input <- list(external_input$x, xvars_prep) %>%
            purrr::pmap(c) %>%
            map(unique) %>%
            map2(external_input$y, ~.x[!.x %in% .y]) %>%
            map(stats::na.omit)

        tibble(
            y = external_input$y,
            x = xvar_input,
            index_node = network_object_tbl$index_node
        ) %>%
            dplyr::distinct()
    }


# Putting this on hold, but keeping for now.
# de_vars <- NA
# exposure_estimates <- list()
# edge_table <- list()
# count <- 1
# while (length(de_vars) > 0) {
#     browser()
#     print(count)
#
#     edge_table[[count]] <- as_edge_tbl(metabolite_network) %>%
#         filter(!source_node %in% de_vars)
#
#     exposure_estimates[[count]] <- standardized_data %>%
#         nc_estimate_exposure_links(
#             edge_tbl = edge_table[[count]],
#             exposure = "exposure",
#             adjustment_vars = "age",
#             .direct_effect_vars = de_vars,
#             model_function = lm
#         ) %>%
#         nc_classify_effects()
#
#     de_vars_current <- exposure_estimates[[count]] %>%
#         filter(effect == "direct") %>%
#         pull(index_node) %>%
#         unique()
#
#     de_vars <- na.omit(c(de_vars, de_vars_current))
#
#     count <- count + 1
# }
