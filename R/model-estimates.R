#' Estimating pathways to an estimated DAG from an exposure.
#'
#' @description
#' \lifecycle{experimental}
#'
#' This algorithm estimates direct effects of a predefined exposure on each
#' network-variable for all causal models that agree with the input-network:
#' models are adjusted for all possible combinations of direct neighbors
#' (==variables in the adjacency set) -> Output is a multiset of possible
#' effects
#'
#' @param .exposure Character. The exposure variable of interest.
#'
#' @return Outputs a [tibble][tibble::tibble-package] with all the models computed and their
#'   estimates.
#'
#' @examples
#'
#' library(dplyr)
#' metabolite_network <- simulated_data %>%
#'   select(matches("metabolite")) %>%
#'   nc_create_network()
#' simulated_data %>%
#'   nc_model_estimates(
#'     .graph = metabolite_network,
#'     .exposure = "exposure",
#'     .adjustment_vars = "age",
#'     .model_function = lm
#'    )
#'
nc_model_estimates <-
    function(.tbl,
             .graph,
             .exposure,
             .adjustment_vars = NA,
             .model_function,
             .exponentiate = FALSE,
             ...) {

    assert_is_data.frame(.tbl)
    assert_is_s4(.graph)
    assert_is_a_string(.exposure)
    # TODO: This check needs to be better constructed
    if (!missing(.adjustment_vars) | !any(is.na(.adjustment_vars)))
        assert_is_character(.adjustment_vars)
    assert_is_function(.model_function)

    edge_list <- .graph@graph@edgeL
    edge_table <- as_edge_tbl(edge_list)

    # .external_var

    all_neighbour_combinations <- function(.edge_table) {
        neighbours <- .edge_table$target_node
        all_combinations <- lapply(
            seq_along(neighbours),
            combn,
            x = neighbours,
            simplify = FALSE
        )
        all_combinations <- unlist(all_combinations, recursive = FALSE)
        return(all_combinations)
    }

    split(edge_table, edge_table$source_node) %>%
        map(all_neighbour_combinations) %>%
        imap_dfr(
            ~ tibble(index_node = .y, neighbour_list = .x) %>%
                dplyr::add_row(index_node = .y, neighbour_list = NULL)
        ) %>%
        mutate(external = .external_var) %>%
        rename_at()

    generate_formula_list <- function(.tbl, .external_side = c("exposure", "outcome"),
                                      .adjustment_vars) {
        # external/index variable could be x or y
        .tbl %>%
            mutate(x_variables = ~ stats::na.omit(c(.)))
        x_variables <- stats::na.omit(c(neighbours, adjustments, .exposure))
        formula <- stats::reformulate(x_variables, response = y)
        return(formula_list)
    }
    all_possible_model_formulas <- network_edges %>%
        imap( ~ stats::reformulate(stats::na.omit(
            c(.graph@graph@nodes[.x$edges], .adjustment_vars, .exposure)),
            response = .y
        ))

    variables_to_keep <- all_possible_model_formulas %>%
        purrr::map(all.vars) %>%
        purrr::flatten_chr() %>%
        unique()

    # TODO: Right now only complete case in full dataset is allowed, need to consider at model-level
    .tbl <- .tbl %>%
        select_at(variables_to_keep) %>%
        stats::na.omit()

    all_possible_models <- all_possible_model_formulas %>%
        map( ~ {
            # TODO: Need to figure out how to pass other options to function
            .model_function(formula = .x,
                            data = .tbl,
                            na.action = "na.fail")
        }) %>%
        map(~ suppressMessages(MuMIn::dredge(.x, fixed = stats::na.omit(c(
                 .adjustment_vars, .exposure
             )))))

    # TODO: Have argument for threshold? For choosing number of models?
    # TODO: Don't know why but now message outputs about new name... might be a _dfr thing?
    tidied_models <- suppressMessages(
        imap_dfr(all_possible_models,
                 .extract_and_tidy_models,
                 .exponentiate = .exponentiate)
    )

    tidied_models <- tidied_models %>%
        mutate(
            exposure = .exposure,
            adjusted_vars = dplyr::if_else(
                !any(is.na(.adjustment_vars)),
                paste(.adjustment_vars, collapse = ", "),
                NA_character_
            )
        ) %>%
            # select_at(vars("index_node", everything())))
        dplyr::relocate(c("exposure", "index_node", "model_id")) %>%
        dplyr::rename_all(~ gsub("\\.", "_", .))

    return(tidied_models)
}

.single_edge_list_to_tbl <- function(.edges, .nodes) {
    tibble(target_node = .nodes[.edges$edges])
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
