#' Plot of the network of the metabolic variables.
#'
#' @description
#' \lifecycle{experimental}
#'
#' @param .tbl The data containing only the metabolic variables.
#' @param .graph The network graph object of the metabolic variable network.
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
#' @examples
#'
#' library(dplyr)
#' metabolite_data <- simulated_data %>%
#'   select(starts_with("metabolite"))
#' network <- metabolite_data %>%
#'   nc_create_network()
#' nc_plot_network(
#'   metabolite_data,
#'   network,
#'   .fn_node_rename = function(x) gsub("metabolite_", "M", x)
#' )
#'
nc_plot_network <- function(.tbl,
                            .graph,
                            .edge_label_threshold = 0.2,
                            .fn_node_rename = NULL) {

    if (is.null(.fn_node_rename))
        .fn_node_rename <- function(x) x
    assert_is_function(.fn_node_rename)

    # TODO: Fix this to be tidier, there should be a better way to do it.
    if (!requireNamespace("ggplot2", quietly = TRUE))
        rlang::abort("Can't find ggplot2, please install it.")
    if (!requireNamespace("tidygraph", quietly = TRUE))
        rlang::abort("Can't find tidygraph, please install it.")

    graph_data_prep <- .tbl %>%
        nc_adjacency_graph(.graph = .graph) %>%
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
        nc_adjacency_graph(.graph = .graph) %>%
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

nc_plot_external_var <- function(.tbl, .graph, .tbl_model, .edge_label_threshold = 0.2) {
    tbl_graph <- .create_tbl_network_graph(.tbl, .graph)

    tbl_model_edges <- .tbl_model %>%
        mutate(
            to = as.numeric(as.factor(.data$index_node)),
            from = length(unique(.data$outcome)) + length(unique(.data$index_node)),
            direct_effect = dplyr::na_if(.data$direct_effect, "none"),
            estimate = if_else(is.na(.data$direct_effect), NA_real_, .data$estimate)
        ) %>%
        select(all_of(c("from", "to", "estimate", "p_value", "direct_effect")))

    tbl_graph_edges <- tbl_graph %>%
        tidygraph::activate("edges") %>%
        as_tibble()

    tbl_edges <- dplyr::bind_rows(tbl_model_edges, tbl_graph_edges)

    tbl_graph_data <- tbl_graph(
        nodes = tibble(name = tbl_graph_edges %>%
                           tidygraph::activate("nodes") %>%
                           dplyr::pull(.data$name) %>%
                           dplyr::union(unique(
                               .tbl_model$outcome
                           ))),
        edges = tbl_edges
    ) %>%
        tidygraph::activate("edges") %>%
        mutate(
            estimate = if_else(.data$direct_effect == "none", NA_real_,
                               .data$estimate),
            weight = if_else(is.na(.data$weight), .data$estimate, .data$weight),
            direct_effect = if_else(is.na(.data$direct_effect),
                                    "direct", .data$direct_effect)
        ) %>%
        .define_edge_label(.edge_label_threshold) %>%
        mutate(edge_label = if_else(is.na(.data$edge_label),
                                    "",
                                    .data$edge_label))
    node_positions <- tbl_graph_data %>%
        ggraph::create_layout("gem") %>%
        mutate(y = if_else(.data$name == unique(.tbl_model$outcome),
                           mean(.data$y), .data$y),
               x = if_else(.data$name == unique(.tbl_model$outcome),
                           max(.data$x) * 2.25, .data$x))

    tbl_graph_data %>%
        ggraph::ggraph("manual", x = node_positions$x, y = node_positions$y) +
        ggraph::geom_edge_diagonal(
            ggplot2::aes_string(
                label = "edge_label",
                colour = "weight",
                width = "abs(weight)",
                linetype = "forcats::fct_rev(direct_effect)",
                alpha = "direct_effect"
            ),
            angle_calc = "along",
            label_dodge = grid::unit(0.25, "cm")
        ) +
        ggraph::geom_node_point(size = 2) +
        ggraph::scale_edge_colour_distiller(palette = "RdBu") +
        ggraph::scale_edge_alpha_discrete(guide = FALSE, range = c(0.7, 1)) +
        ggraph::scale_edge_width(guide = FALSE, range = c(0.75, 2)) +
        ggraph::geom_node_text(ggplot2::aes_string(label = "name"),
                               repel = TRUE) +
        ggraph::scale_edge_linetype_discrete(name = "") +
        ggraph::theme_graph()
}

