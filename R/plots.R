#' Plot of the network of the metabolic variables.
#'
#' @description
#' \lifecycle{experimental}
#'
#' @param .tbl The data that was also used to generate the [nc_estimate_network()].
#' @param .graph The network graph object of the metabolic variable network.
#' @param .vars <[`tidy-select`][dplyr::dplyr_tidy_select]> Variables to include
#'   by using [dplyr::select()] style selection.
#' @param .edge_label_threshold Threshold set for edge weight value above which
#'   the edge label will be kept. This argument helps to reduce clutter in the
#'   graph.
#' @param .fn_node_rename Function to pass to rename the metabolic variables.
#'   Preferably use functions that search and replace patterns, like [gsub()] or
#'   [stringr::str_replace()].
#'
#' @return Outputs a `ggplot2` object of the metabolic network.
#' @export
#'
#' @seealso [nc_model_estimates]
#'
nc_plot_network <- function(.tbl,
                            .graph,
                            .vars = everything(),
                            .edge_label_threshold = 0.2,
                            .fn_node_rename = NULL) {

    if (is.null(.fn_node_rename))
        .fn_node_rename <- function(x) x
    assert_is_function(.fn_node_rename)
    .tbl <- select(.tbl, {{ .vars }})
    node_names <- unique(as_edge_tbl(.graph)[[1]])
    if (!all(names(.tbl) %in% node_names))
        rlang::abort("Both the data (given with `.tbl` and `.vars`) and the network graph must have the same variables. Right now they are different.")

    # TODO: Fix this to be tidier, there should be a better way to do it.
    if (!requireNamespace("ggplot2", quietly = TRUE))
        rlang::abort("Can't find ggplot2, please install it.")
    if (!requireNamespace("tidygraph", quietly = TRUE))
        rlang::abort("Can't find tidygraph, please install it.")

    graph_data_prep <- .tbl %>%
        compute_adjacency_graph(.graph = .graph) %>%
        tidygraph::as_tbl_graph() %>%
        tidygraph::activate("edges") %>%
        tidygraph::mutate(edge_label = dplyr::if_else(
            abs(.data$weight) > .edge_label_threshold,
            as.character(round(.data$weight, 2)),
            ""
        ))

    graph_data_prep %>%
        ggraph::ggraph("stress") +
        ggraph::geom_edge_diagonal(
            ggplot2::aes_string(
                label = "edge_label",
                colour = "weight",
                width = "abs(weight)"
            ),
            angle_calc = "along",
            label_dodge = grid::unit(0.2, "cm")
        ) +
        ggraph::geom_node_point(size = 2) +
        ggraph::scale_edge_colour_gradient2(mid = "gray80") +
        ggraph::scale_edge_width(guide = FALSE, range = c(0.75, 2)) +
        ggraph::geom_node_text(ggplot2::aes_string(label = ".fn_node_rename(name)"),
                               repel = TRUE) +
        ggraph::theme_graph(base_family = 'Helvetica')
}

.create_tbl_network_graph <- function(.tbl, .graph) {
    .tbl %>%
        select(all_of(.graph@graph@nodes)) %>%
        compute_adjacency_graph(.graph = .graph) %>%
        tidygraph::as_tbl_graph() %>%
        tidygraph::activate("edges")
}

.define_edge_label <- function(.tbl_graph, .edge_label_threshold = 0.2) {
    .tbl_graph %>%
        mutate(edge_label =
                   if_else(
                       abs(.data$weight) > .edge_label_threshold,
                       as.character(round(.data$weight, 2)),
                       ""
                   ))
}

plot_external_var <-
    function(.tbl,
             .graph,
             .tbl_model,
             .edge_label_threshold = 0.2,
             .external_var_side =c("outcome", "exposure")) {

    if (!requireNamespace("ggplot2", quietly = TRUE))
        rlang::abort("Can't find ggplot2, please install it.")
    if (!requireNamespace("tidygraph", quietly = TRUE))
        rlang::abort("Can't find tidygraph, please install it.")

    external_var <- match.arg(.external_var_side)
    # TODO: Extract the data processing from the plotting functionality
    tbl_graph <- .create_tbl_network_graph(.tbl, .graph)

    tbl_model_edges <- .tbl_model %>%
        mutate(
            to = as.numeric(as.factor(.data$index_node)),
            from = length(unique(.data[[external_var]])) + length(unique(.data$index_node)),
            effect = dplyr::na_if(.data$effect, "none"),
            estimate = if_else(is.na(.data$effect), NA_real_, .data$estimate)
        ) %>%
        select(all_of(c("from", "to", "estimate", "p_value", "effect")))

    tbl_graph_edges <- tbl_graph %>%
        tidygraph::activate("edges")

    tbl_edges <- dplyr::bind_rows(tbl_model_edges, as_tibble(tbl_graph_edges))

    tbl_graph_data <- tidygraph::tbl_graph(
        nodes = tibble(name = tbl_graph_edges %>%
                           tidygraph::activate("nodes") %>%
                           dplyr::pull(.data$name) %>%
                           dplyr::union(unique(
                               .tbl_model[[external_var]]
                           ))),
        edges = tbl_edges
    ) %>%
        tidygraph::activate("edges") %>%
        mutate(
            estimate = if_else(.data$effect == "none", NA_real_,
                               .data$estimate),
            weight = if_else(is.na(.data$weight), .data$estimate, .data$weight),
            effect = if_else(is.na(.data$effect), "direct", .data$effect)
        ) %>%
        .define_edge_label(.edge_label_threshold) %>%
        mutate(edge_label = if_else(is.na(.data$edge_label),
                                    "",
                                    .data$edge_label))

    nudge_to_side <- switch(.external_var_side,
                            exposure = min,
                            outcome = max)

    node_positions <- tbl_graph_data %>%
        ggraph::create_layout("gem") %>%
        mutate(y = if_else(.data$name == unique(.tbl_model[[external_var]]),
                           mean(.data$y), .data$y),
               x = if_else(.data$name == unique(.tbl_model[[external_var]]),
                           nudge_to_side(.data$x) * 2.25, .data$x))

    # TODO: Convert this into own geom object?
    tbl_graph_data %>%
        ggraph::ggraph("manual", x = node_positions$x, y = node_positions$y) +
        ggraph::geom_edge_diagonal(
            ggplot2::aes_string(
                label = "edge_label",
                colour = "weight",
                width = "abs(weight)",
                linetype = ".fct_rev(effect)",
                alpha = "effect"
            ),
            angle_calc = "along",
            label_dodge = grid::unit(0.25, "cm")
        ) +
        ggraph::geom_node_point(size = 2) +
        # ggraph::scale_edge_colour_distiller(palette = "RdBu") +
        ggraph::scale_edge_colour_gradient2(mid = "gray80") +
        ggraph::scale_edge_alpha_discrete(guide = FALSE, range = c(0.7, 1)) +
        ggraph::scale_edge_width(guide = FALSE, range = c(0.75, 2)) +
        ggraph::geom_node_text(ggplot2::aes_string(label = "name"),
                               repel = TRUE) +
        ggraph::scale_edge_linetype_discrete(name = "") +
        ggraph::theme_graph()
    }

#' Plot of the outcome or exposure side model estimation.
#'
#' @description
#' \lifecycle{experimental}
#'
#' Plots the results of the effect classification based on the model estimation,
#' linking the results into the network graph.
#'
#' @param .tbl The original data, with the metabolic variables that have been
#'   standardized.
#' @param .graph The graph object created from `nc_estimate_network()`.
#' @param .tbl_model The tibble of the model results obtained from
#'   `nc_classify_effects()`.
#' @param .edge_label_threshold Threshold to pass for the value to be added to
#'   the edge label.
#'
#' @return a [ggplot2][ggplot2::ggplot2-package] object showing the model
#'   estimation results linked with the network graph.
#' @export
#'
nc_plot_outcome_estimation <- function(.tbl,
                                       .graph,
                                       .tbl_model,
                                       .edge_label_threshold = 0.2) {
    plot_external_var(
        .tbl = .tbl,
        .graph = .graph,
        .tbl_model = .tbl_model,
        .edge_label_threshold = .edge_label_threshold,
        .external_var_side = "outcome"
    )

}

#' @describeIn nc_plot_outcome_estimation Plots the exposure side estimation.
#' @export
nc_plot_exposure_estimation <- function(.tbl,
                                        .graph,
                                        .tbl_model,
                                        .edge_label_threshold = 0.2) {
    plot_external_var(
        .tbl = .tbl,
        .graph = .graph,
        .tbl_model = .tbl_model,
        .edge_label_threshold = .edge_label_threshold,
        .external_var_side = "exposure"
    )
}
