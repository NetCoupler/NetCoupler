#' @title
#' Create an estimate of the metabolic network as an undirected graph.
#'
#' @description
#' \lifecycle{experimental}
#'
#' The main NetCoupler network creator.
#' Uses the input data to estimate the underlying undirected graph.
#' The default uses the PC algorithm, implemented within NetCoupler
#' with [pc_estimate_undirected_graph()]
#' Defaults to using the PC algorithm to calculate possible edges.
#' Any missing values in the input data are removed by this function,
#' since some computations can't handle missingness.
#'
#' @param data Data that would form the underlying network.
#' @param cols <[`tidy-select`][dplyr::dplyr_tidy_select]> Variables to include
#'   by using [dplyr::select()] style selection.
#' @param alpha The alpha level to use to test whether an edge exists or not.
#'   Default is 0.01.
#'
#' @return Outputs a [tidygraph::tbl_graph()] with the start and end nodes, as
#'   well as the edge weights.
#' @export
#'
#' @seealso See [nc_estimate_links] for examples on using NetCoupler and
#'   [pc_estimate_undirected_graph] for more details on the PC-algorithm network
#'   estimation method.
#'
nc_estimate_network <- function(data, cols = everything(), alpha = 0.01) {
    assert_data_frame(data)
    assert_number(alpha)

    subset_data <- data %>%
        select({{cols}}) %>%
        na.omit()

    # TODO: Not sure why but this only has 9 of the 12 variables
    tbl_network <- subset_data %>%
        pc_estimate_undirected_graph(alpha) %>%
        as_tbl_graph.pcAlgo()

    compute_weighted_adjacency_graph(
        subset_data,
        tbl_network
    )
}

#' Convert network graphs to edge tables as tibbles/data.frames.
#'
#' @description
#' \lifecycle{experimental}
#'
#' @param network_object Network graph from e.g. [nc_estimate_network()].
#'
#' @return A [tibble][tibble::tibble-package], with at least two columns:
#'
#' - `source_node`: The starting node (variable).
#' - `target_node`: The ending node (variable) that links with the source node.
#' - `adjacency_weight`: (Optional) The "weight" given to the edge, which
#' represents the strength of the link between two nodes.
#'
#' @export
#'
#' @seealso See [nc_estimate_links] for examples on using NetCoupler.
#'
as_edge_tbl <- function(network_object) {
    UseMethod("as_edge_tbl", network_object)
}

#' @export
as_edge_tbl.tbl_graph <- function(network_object) {
    dplyr::bind_rows(
        network_object %>%
            igraph::as_data_frame() %>%
            dplyr::rename(source_node = from, target_node = to,
                   adjacency_weight = weight),
        network_object %>%
            igraph::as_data_frame() %>%
            dplyr::rename(source_node = to, target_node = from,
                   adjacency_weight = weight)
    ) %>%
        unique()
}

#' @export
as_edge_tbl.default <- function(network_object) {
    if (checkmate::test_data_frame(network_object))
        return(network_object)
    else
        rlang::abort("We don't know how to handle the object given as `network_object`. This function only can accept `tbl_graph()` objects for now.")
}

#' Estimate the undirected graph of the metabolic data.
#'
#' @description
#' Uses the PC-algorithm and is mostly a wrapper around [pcalg::skeleton()].
#'
#' @details
#' This function estimates the "skeleton of a DAG", meaning a graph without
#' arrowheads, aka an undirected graph.
#' The default estimation method used is the "PC-stable" method, which estimates
#' the *order-independent* skeleton of the DAG, meaning the order of the
#' variables given does not impact the results (older versions of the algorithm
#' were order-dependent). The method also assumes no latent variables.
#'
#' An edge is determined by testing for conditional dependence between two
#' nodes based on the [pcalg::gaussCItest()]. Conditional *independence* exists
#' when the nodes have zero partial correlation determined from a p-value based
#' hypothesis test against the correlation matrix of the data from the nodes.
#' The estimated edges exists between the *start* and *end* nodes when the
#' *start* and *end* variables are conditionally dependent given the subset of
#' remaining variables.
#'
#' @param data Input numeric data that forms the basis of the underlying graph.
#' @param alpha Significance level threshold applied to each test to determine
#'   conditional dependence for if an edge exists.
#'
#' @return A `pcAlgo` object that contains the DAG skeleton, aka undirected graph.
#' @keywords internal
#' @seealso The help documentation of [pcalg::skeleton()] has more details.
#'
pc_estimate_undirected_graph <- function(data, alpha = 0.01) {
    number_samples <- nrow(data)
    metabolite_names <- colnames(data)

    pcalg::skeleton(
        suffStat = list(C = stats::cor(data, use = "complete.obs"), n = number_samples),
        indepTest = pcalg::gaussCItest,
        labels = metabolite_names,
        method = "stable",
        alpha = alpha,
        fixedGaps = NULL,
        fixedEdges = NULL,
        verbose = FALSE
    )
}

# Helpers -----------------------------------------------------------------

# Convert the pcAlgo object to a tidygraph tbl_graph object.
as_tbl_graph.pcAlgo <- function(pc_graph) {
    pc_graph %>%
        pcalg::getGraph() %>%
        igraph::graph_from_graphnel() %>%
        tidygraph::as_tbl_graph()
}

#' Compute the weighted adjacency matrix and use to create the graph with weighting.
#'
#' @param data Input data.
#' @param network_graph The output object from [nc_estimate_network()].
#'
#' @return Outputs a [tidygraph::tbl_graph()] object.
#' @keywords internal
#' @noRd
#'
compute_weighted_adjacency_graph <- function(data, network_graph) {
    # Calculate the weighted adjacency matrix
    # Rounding fixes a problem with very small numbers,
    # this forces them to be zero.
    weighted_adj_matrix <- igraph::as_adjacency_matrix(network_graph) *
        round(compute_partial_corr_matrix(data), digits = 3)

    weighted_adj_matrix %>%
        as.matrix() %>%
        igraph::graph_from_adjacency_matrix(weighted = TRUE,
                                            mode = "undirected") %>%
        tidygraph::as_tbl_graph()
}

#' Estimate Pearson's partial correlation coefficients.
#'
#' This function is a wrapper around [ppcor::pcor()] that extracts correlation
#' coefficient estimates, then adds the variable names to the column and row names.
#'
#' @param data Input data of metabolic variables as matrix or data.frame.
#'
#' @return Outputs a matrix of partial correlation coefficients.
#' @keywords internal
#' @noRd
#'
compute_partial_corr_matrix <- function(data) {
    pcor_matrix <- ppcor::pcor(data)$estimate
    colnames(pcor_matrix) <- colnames(data)
    rownames(pcor_matrix) <- colnames(data)
    pcor_matrix
}
