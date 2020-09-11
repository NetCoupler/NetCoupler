
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
        round(nc_partial_corr_matrix(.tbl), digits = 3)

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

#' Estimate equivalence class of DAG from the PC algorithm.
#'
#' Is mostly a wrapper around [pcalg::pc()]. Estimates an order-independent
#' skeleton.
#'
#' @param .tbl Input data, samples by metabolite matrix or as data.frame.
#' @param .alpha Significance level threshold applied to each test.
#'
#' @return Outputs a `pcAlgo` object.
#'
#' @examples
#'
#' \dontrun{
#' library(dplyr)
#' simulated_data %>%
#' select(contains("metabolite")) %>%
#' NetCoupler:::pc_dag_estimates()
#' }
#'
pc_dag_estimates <- function(.tbl, .alpha = 0.01) {
    number_samples <- nrow(.tbl)
    metabolite_names <- colnames(.tbl)

    pcalg::pc(
        suffStat = list(C = stats::cor(.tbl), n = number_samples),
        indepTest = pcalg::gaussCItest,
        labels = metabolite_names,
        skel.method = "stable",
        alpha = .alpha,
        fixedGaps = NULL,
        fixedEdges = NULL,
        verbose = FALSE,
        maj.rule = FALSE,
        solve.confl = FALSE
    )
}

#' Estimate order-independent PC-stable skeleton of a DAG.
#'
#' Uses the PC-algorithm and is mostly a wrapper around [pcalg::skeleton()].
#'
#' @param .tbl Input metabolic data.
#' @param .alpha Significance level threshold applied to each test.
#'
#' @return A DAG skeleton object.
#'
#' @examples
#'
#' \dontrun{
#' library(dplyr)
#' simulated_data %>%
#' select(contains("metabolite")) %>%
#' NetCoupler:::pc_skeleton_estimates()
#' }
#'
pc_skeleton_estimates <- function(.tbl, .alpha = 0.01) {
    number_samples <- nrow(.tbl)
    metabolite_names <- colnames(.tbl)

    # TODO: Confirm that this does this.
    pcalg::skeleton(
        suffStat = list(C = stats::cor(.tbl, use = "complete.obs"), n = number_samples),
        # Test conditional independence of Gaussians via Fisher's Z
        indepTest = pcalg::gaussCItest,
        labels = metabolite_names,
        method = "stable",
        alpha = .alpha,
        fixedGaps = NULL,
        fixedEdges = NULL,
        verbose = FALSE
    )
}
