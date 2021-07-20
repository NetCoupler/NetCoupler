#' Plot of the network of the metabolic variables.
#'
#' @description
#' \lifecycle{experimental}
#'
#' @param data The data that was also used to generate the [nc_estimate_network()].
#' @param edge_tbl The network graph object of the metabolic variable network.
#' @param edge_label_threshold Threshold set for edge weight value above which
#'   the edge label will be kept. This argument helps to reduce clutter in the
#'   graph.
#' @param fn_node_rename Function to pass to rename the metabolic variables.
#'   Preferably use functions that search and replace patterns, like [gsub()] or
#'   [stringr::str_replace()].
#'
#' @return Outputs a `ggplot2` object of the metabolic network.
#' @export
#'
#' @seealso See [nc_estimate_links] for examples on using NetCoupler.
#'
nc_plot_network <- function(data,
                            edge_tbl,
                            edge_label_threshold = 0.2,
                            fn_node_rename = NULL) {

    if (is.null(fn_node_rename))
        fn_node_rename <- function(x) x
    assert_is_function(fn_node_rename)

    graph_data_prep <- nc_tbl_adjacency_graph(data, edge_tbl) %>%
        tidygraph::activate("edges") %>%
        tidygraph::mutate(edge_label = dplyr::if_else(
            abs(.data$weight) > edge_label_threshold,
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
        ggraph::scale_edge_colour_gradient2(mid = "gray80", limits = c(-1, 1)) +
        ggraph::scale_edge_width(guide = FALSE, range = c(0.75, 2)) +
        ggraph::geom_node_text(ggplot2::aes_string(label = "fn_node_rename(name)"),
                               repel = TRUE) +
        ggraph::theme_graph(base_family = 'Helvetica')
}

plot_external_var <-
    function(data,
             edge_tbl,
             data_model,
             edge_label_threshold = 0.2,
             external_var_side =c("outcome", "exposure")) {

    external_var <- rlang::arg_match(external_var_side)

    tbl_graph_edges <- nc_tbl_adjacency_graph(data, edge_tbl) %>%
        tidygraph::activate("edges")

    tbl_model_edges <- convert_model_data_to_model_edges(data_model, external_var)

    tbl_edges <- dplyr::bind_rows(tbl_model_edges, as_tibble(tbl_graph_edges))

    tbl_graph_data <- tidygraph::tbl_graph(
        nodes = node_names_as_tibble(tbl_graph_edges, data_model[[external_var]]),
        edges = tbl_edges
    ) %>%
        set_weights_and_tidy() %>%
        define_edge_label(edge_label_threshold) %>%
        mutate(edge_label = if_else(
            is.na(.data$edge_label) | .data$from == max(.data$from),
            "",
            .data$edge_label
        ))

    nudge_to_side <- switch(external_var_side,
                            exposure = min,
                            outcome = max)

    node_positions <- tbl_graph_data %>%
        ggraph::create_layout("stress") %>%
        mutate(y = if_else(.data$name == unique(data_model[[external_var]]),
                           mean(.data$y), .data$y),
               # Shift x axis over so nudge works
               x = scale(.data$x, scale = FALSE),
               x = if_else(.data$name == unique(data_model[[external_var]]),
                           nudge_to_side(.data$x) * 2, .data$x))

    # TODO: Convert this into own geom object?
    tbl_graph_data %>%
        ggraph::ggraph("manual", x = node_positions$x, y = node_positions$y) +
        ggraph::geom_edge_diagonal(
            ggplot2::aes_string(
                label = "edge_label",
                colour = "weight",
                width = "abs(weight)",
                linetype = "fct_rev(effect)",
                alpha = "effect"
            ),
            angle_calc = "along",
            label_dodge = grid::unit(0.25, "cm")
        ) +
        ggraph::geom_node_point(size = 2) +
        # ggraph::scale_edge_colour_distiller(palette = "RdBu") +
        ggraph::scale_edge_colour_gradient2(mid = "gray80", limit = c(-1, 1)) +
        ggraph::scale_edge_alpha_discrete(guide = FALSE, range = c(0.7, 0.9)) +
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
#' @param data The original data, with the metabolic variables that have been
#'   standardized.
#' @param edge_tbl The graph object created from `nc_estimate_network()`.
#' @param data_model The tibble of the model results obtained from
#'   `nc_classify_effects()`.
#' @param edge_label_threshold Threshold to pass for the value to be added to
#'   the edge label.
#'
#' @return a [ggplot2][ggplot2::ggplot2-package] object showing the model
#'   estimation results linked with the network graph.
#' @export
#'
nc_plot_outcome_estimation <- function(data,
                                       edge_tbl,
                                       data_model,
                                       edge_label_threshold = 0.2) {
    plot_external_var(
        data = data,
        edge_tbl = edge_tbl,
        data_model = data_model,
        edge_label_threshold = edge_label_threshold,
        external_var_side = "outcome"
    )

}

#' @describeIn nc_plot_outcome_estimation Plots the exposure side estimation.
#' @export
nc_plot_exposure_estimation <- function(data,
                                        edge_tbl,
                                        data_model,
                                        edge_label_threshold = 0.2) {
    plot_external_var(
        data = data,
        edge_tbl = edge_tbl,
        data_model = data_model,
        edge_label_threshold = edge_label_threshold,
        external_var_side = "exposure"
    )
}

# Helpers -----------------------------------------------------------------

convert_model_data_to_model_edges <- function(data, ext_var) {
    data %>%
        dplyr::arrange(.data$index_node) %>%
        mutate(
            to = as.numeric(as.factor(.data$index_node)),
            from = length(unique(.data[[ext_var]])) + length(unique(.data$index_node)),
            effect = dplyr::na_if(.data$effect, "none"),
            estimate = if_else(is.na(.data$effect), NA_real_, .data$estimate)
        ) %>%
        dplyr::filter(!is.na(.data$effect)) %>%
        select(all_of(c("from", "to", "estimate", "effect")))
}

set_weights_and_tidy <- function(data_graph) {
    data_graph %>%
        tidygraph::activate("edges") %>%
        mutate(
            estimate = dplyr::case_when(
                .data$effect == "ambiguous" & .data$estimate > 0 ~ 0.3,
                .data$effect == "ambiguous" & .data$estimate < 0 ~ -0.3,
                .data$effect == "direct" & .data$estimate > 0 ~ 0.7,
                .data$effect == "direct" & .data$estimate < 0 ~ -0.7
            ),
            weight = if_else(is.na(.data$weight), .data$estimate, .data$weight),
            effect = if_else(is.na(.data$effect), "direct", .data$effect)
        )
}

node_names_as_tibble <- function(data_edges, ext_vars) {
    tibble(
        name = data_edges %>%
            tidygraph::activate("nodes") %>%
            dplyr::pull(.data$name) %>%
            dplyr::union(unique(ext_vars))
    )
}
