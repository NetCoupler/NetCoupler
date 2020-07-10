
#' @title
#' Create an estimate of the metabolic network skeleton.
#'
#' @description
#' \lifecycle{experimental}
#'
#' Main NetCoupler network creator.
#' Estimates the skeleton based on a family of DAGs without specifying the direction of edges.
#' Defaults to using the PC algorithm to calculate possible edges.
#'
#' @param .tbl Data of the metabolic variables.
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
nc_create_network <- function(.tbl, .alpha = 0.05) {
    assert_is_data.frame(.tbl)
    assert_is_a_number(.alpha)

    pc_skeleton_estimates(.tbl, .alpha)
}

#' Compute the adjacency matrix of the graph with the data.
#'
#' @description
#' \lifecycle{experimental}
#'
#' @inheritParams nc_plot_network
#'
#' @return Outputs an `igraph` object from [igraph::graph.adjacency()].
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' metabolite_data <- simulated_data %>%
#'   select(starts_with("metabolite"))
#' network <- metabolite_data %>%
#'   nc_create_network()
#' nc_adjacency_graph(metabolite_data, network)
nc_adjacency_graph <- function(.tbl, .graph) {
    weighted_adjacency_matrix <- nc_adjacency_matrix(.graph) *
        nc_partial_corr_matrix(.tbl)

    igraph::graph.adjacency(weighted_adjacency_matrix,
                            weighted = TRUE, mode = "undirected")
}

#' Extract adjacency matrix from a DAG skeleton.
#'
#' @description
#' \lifecycle{experimental}
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
#' @description
#' \lifecycle{experimental}
#'
#' This function is a wrapper around [ppcor::pcor()] that extracts correlation
#' coefficient estimates, then adds the variable names to the column and row names.
#'
#' @param .tbl Input data of metabolic variables as matrix or data.frame.
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
nc_partial_corr_matrix <- function(.tbl) {
    # TODO: check if input data is gaussian
    pcor_matrix <- ppcor::pcor(.tbl)$estimate
    colnames(pcor_matrix) <- colnames(.tbl)
    rownames(pcor_matrix) <- colnames(.tbl)
    return(pcor_matrix)
}
