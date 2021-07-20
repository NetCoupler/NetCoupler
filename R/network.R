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
#' @param data Data of the metabolic variables.
#' @param cols <[`tidy-select`][dplyr::dplyr_tidy_select]> Variables to include
#'   by using [dplyr::select()] style selection.
#' @param alpha The alpha level to set. Default is 0.05.
#'
#' @return Outputs a DAG skeleton.
#' @export
#'
#' @seealso See [nc_estimate_links] for examples on using NetCoupler.
#'
nc_estimate_network <- function(data, cols = everything(), alpha = 0.05) {
    assert_is_data.frame(data)
    assert_is_a_number(alpha)

    data %>%
        select({{cols}}) %>%
        pc_skeleton_estimates(alpha)
}

#' Convert network graph to edge table.
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
#'
#' @export
#'
#' @seealso See [nc_estimate_links] for examples on using NetCoupler.
#'
as_edge_tbl <- function(network_object) {
    UseMethod("as_edge_tbl", network_object)
}

#' @export
as_edge_tbl.igraph <- function(network_object) {
    network_object %>%
        tidygraph::as_tbl_graph()
}

#' @export
as_edge_tbl.tbl_graph <- function(network_object) {
    nodes <- network_object %>%
        tidygraph::activate("nodes") %>%
        mutate(id = dplyr::row_number(.data$name)) %>%
        tidygraph::as_tibble()

    edges <- network_object %>%
        tidygraph::activate("edges") %>%
        tidygraph::as_tibble()

    tibble(
        source_node = edges %>%
            dplyr::left_join(nodes, by = c("from" = "id")) %>%
            dplyr::pull(.data$name),
        target_node = edges %>%
            dplyr::left_join(nodes, by = c("to" = "id")) %>%
            dplyr::pull(.data$name),
        adjacency_weight = edges$weight
    )
}


#' @export
as_edge_tbl.data.frame <- function(network_object) {
    # if (names(network_object) %in% c("source_node", "target_node"))
}

#' @export
as_edge_tbl.default <- function(network_object) {
    rlang::abort("The `network_object` object is not from the pcalg package. We currently do not have support for other network packages.")
}

#' @export
as_edge_tbl.pcAlgo <- function(network_object) {
    network_object <- network_object@graph@edgeL
    nodes <- names(network_object)
    edge_table <- purrr::map_dfr(
        network_object,
        single_network_to_tbl,
        .id = "source_node",
        nodes = nodes
    )
    return(edge_table)
}

#' Compute the adjacency matrix of the graph with the data.
#'
#' @description
#' \lifecycle{experimental}
#'
#' @inheritParams nc_plot_network
#'
#' @return Outputs an `igraph` object from [igraph::graph_from_adjacency_matrix()].
#' @keywords internal
#'
compute_adjacency_graph <- function(data, edge_tbl) {
    # TODO: This may change underlying graph connections, check into this.
    weighted_adjacency_matrix <- compute_adjacency_matrix(edge_tbl) *
        round(compute_partial_corr_matrix(data), digits = 3)

    igraph::graph_from_adjacency_matrix(
        as.matrix(weighted_adjacency_matrix),
        weighted = TRUE,
        mode = "undirected"
    )
}

#' Extract adjacency matrix from a DAG skeleton.
#'
#' @description
#' \lifecycle{experimental}
#'
#' Is generally a wrapper around calls to [igraph::as_adjacency_matrix()] and
#' [igraph::graph_from_graphnel()]. Transforms from a GraphNEL object in igraph.
#'
#' @param dag_skeleton The PC DAG skeleton object.
#'
#' @return Outputs an adjacency matrix of the DAG skeleton.
#' @keywords internal
#'
compute_adjacency_matrix <- function(dag_skeleton) {
    # TODO: Include a check here that it is a DAG skeleton..?
    from_skeleton <- igraph::graph_from_graphnel(dag_skeleton@graph)
    igraph::as_adjacency_matrix(from_skeleton)
}

#' Estimate Pearson's partial correlation coefficients.
#'
#' @description
#' \lifecycle{experimental}
#'
#' This function is a wrapper around [ppcor::pcor()] that extracts correlation
#' coefficient estimates, then adds the variable names to the column and row names.
#'
#' @param data Input data of metabolic variables as matrix or data.frame.
#'
#' @return Outputs a matrix of partial correlation coefficients.
#' @keywords internal
#'
compute_partial_corr_matrix <- function(data) {
    pcor_matrix <- ppcor::pcor(data)$estimate
    colnames(pcor_matrix) <- colnames(data)
    rownames(pcor_matrix) <- colnames(data)
    return(pcor_matrix)
}

# TODO: Don't know what this does or is for.
#' Estimate equivalence class of DAG from the PC algorithm.
#'
#' Is mostly a wrapper around [pcalg::pc()]. Estimates an order-independent
#' skeleton.
#'
#' @param data Input data, samples by metabolite matrix or as data.frame.
#' @param alpha Significance level threshold applied to each test.
#'
#' @return Outputs a `pcAlgo` object.
#' @keywords internal
#'
pc_dag_estimates <- function(data, alpha = 0.01) {
    number_samples <- nrow(data)
    metabolite_names <- colnames(data)

    pcalg::pc(
        suffStat = list(C = stats::cor(data), n = number_samples),
        indepTest = pcalg::gaussCItest,
        labels = metabolite_names,
        skel.method = "stable",
        alpha = alpha,
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
#' @param data Input metabolic data.
#' @param alpha Significance level threshold applied to each test.
#'
#' @return A DAG skeleton object.
#' @keywords internal
#'
pc_skeleton_estimates <- function(data, alpha = 0.01) {
    number_samples <- nrow(data)
    metabolite_names <- colnames(data)

    # TODO: Confirm that this does this.
    pcalg::skeleton(
        suffStat = list(C = stats::cor(data, use = "complete.obs"), n = number_samples),
        # Test conditional independence of Gaussians via Fisher's Z
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

single_network_to_tbl <- function(edges, nodes) {
    tibble(target_node = nodes[edges$edges])
}
