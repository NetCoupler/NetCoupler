
#' Create an estimate of the metabolic network skeleton.
#'
#' Estimates the skeleton based on a family of DAGs without specifying the direction of edges.
#'
#' \lifecycle{experimental}
#'
#' Main NetCoupler network creator. Defaults to using the PC algorithm to calculate
#' possible edges.
#'
#' @param .data Data of the metabolic variables.
#' @param .alpha The alpha level to set.
#'
#' @return Outputs a DAG skeleton.
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' simulated_data %>%
#'   select(contains("metabolite")) %>%
#'   nc_create_network()
nc_create_network <- function(.data, .alpha = 0.05) {
    assert_is_data.frame(.data)
    assert_is_a_number(.alpha)

    # TODO: Determine if this is important.
    # pcor_matrix <- partial_corr_matrix(.data)

    # skel_est <- pc_skeleton_estimates(.data, .alpha)

    # TODO: Determine if this is important.
    # DAG_est <- pc_dag_estimates(.data, .alpha)

    # TODO: Determine if this is important.
    # adj_matrix <- adjacency_matrix(skel_est)

    # TODO: What does the w stand for?
    # TODO: Determine if this is important.
    # w_adj_matrix <- adj_matrix * pcor_matrix

    # TODO: Rename these to be more descriptive
    # network_list <- list(
        # TODO: Is this necessary to output?
        # pCor_mat = pcor_matrix,
        # skel_est = skel_est
        # TODO: Is this necessary to output?
        # DAG_est = DAG_est,
        # TODO: Is this necessary to output?
        # adj_matrix = adj_matrix,
        # TODO: Is this necessary to output?
        # w_adj_matrix = w_adj_matrix,
        # TODO: Is this necessary to output?
        # data = .data
    # )
    # return(network_list)
    pc_skeleton_estimates(.data, .alpha)
}


#' Compute the adjacency matrix of the graph with the data.
#'
#' @inheritParams nc_plot_network
#'
#' @return Outputs an `igraph` object from [igraph::graph.adjacency()].
#'
#' @examples
#'
#' library(dplyr)
#' metabolite_data <- simulated_data %>%
#'   select(starts_with("metabolite"))
#' network <- simulated_data %>%
#'   select(contains("metabolite")) %>%
#'   nc_create_network()
#' nc_adjacency_graph(metabolite_data, network) %>% class()
nc_adjacency_graph <- function(.data, .graph) {
    weighted_adjacency_matrix <- nc_adjacency_matrix(.graph) *
        nc_partial_corr_matrix(.data)

    igraph::graph.adjacency(weighted_adjacency_matrix,
                            weighted = TRUE, mode = "undirected")
}

#' Plots the network of the metabolic variables.
#'
#' @param .data
#' @param .graph
#'
#' @return Outputs a `ggplot2`
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(NetCoupler)
#' metabolite_data <- simulated_data %>%
#'   select(starts_with("metabolite"))
#' network <- simulated_data %>%
#'   select(contains("metabolite")) %>%
#'   nc_create_network()
#' nc_plot_network(metabolite_data, network)
#'
nc_plot_network <- function(.data, .graph) {
    .data %>%
        nc_adjacency_graph(.graph = .graph) %>%
        tidygraph::as_tbl_graph() %>%
        # tidygraph::activate(edges) %>%
        ggraph::ggraph("stress") +
        ggraph::geom_edge_bend(ggplot2::aes_string(
            colour = "weight",
            width = "abs(weight)"
        )) +
        ggraph::geom_node_point(size = 2) +
        ggraph::scale_edge_colour_gradient2(mid = "gray80") +
        ggraph::scale_edge_width(guide = FALSE, range = c(0.75, 2)) +
        ggraph::geom_node_label(ggplot2::aes_string(label = "name"), repel = TRUE) +
        ggraph::theme_graph()
}

#' Extract adjacency matrix from a DAG skeleton.
#'
#' Is generally a wrapper around calls to [igraph::get.adjacency()] and
#' [igraph::igraph.from.graphNEL()]. Transforms from a GraphNEL object in igraph.
#'
#' @param .dag_skeleton The PC DAG skeleton object.
#'
#' @return Outputs an adjacency matrix of the DAG skeleton.
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' skeleton_estimate <- simulated_data %>%
#'   select(contains("metabolite")) %>%
#'   nc_create_network()
#'
#' nc_adjacency_matrix(skeleton_estimate)
#'
nc_adjacency_matrix <- function(.dag_skeleton) {
    # TODO: Include a check here that it is a DAG skeleton..?
    igraph::get.adjacency(igraph::igraph.from.graphNEL(.dag_skeleton@graph))
}

#' Estimate Pearson's partial correlation coefficients.
#'
#' This function is a wrapper around [ppcor::pcor()] that extracts correlation
#' coefficient estimates, then adds the variable names to the column and row names.
#'
#' @param .data Input data of metabolic variables as matrix or data.frame.
#'
#' @return Outputs a matrix of partial correlation coefficients.
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' simulated_data %>%
#'   select(contains("metabolite")) %>%
#'   nc_partial_corr_matrix()
#'
nc_partial_corr_matrix <- function(.data) {
    # TODO: check if input data is gaussian
    pcor_matrix <- ppcor::pcor(.data)$estimate
    colnames(pcor_matrix) <- colnames(.data)
    rownames(pcor_matrix) <- colnames(.data)
    return(pcor_matrix)
}
