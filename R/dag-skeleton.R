#' Estimate equivalence class of DAG from the PC algorithm.
#'
#' Is mostly a wrapper around [pcalg::pc()]. Estimates an order-independent
#' skeleton.
#'
#' @param .data Input data, samples by metabolite matrix or as data.frame.
#' @param .alpha Significance level threshold applied to each test.
#'
#' @return Outputs a `pcAlgo` object.
#'
#' @examples
#'
#' library(dplyr)
#' simulated_data %>%
#' select(contains("metabolite")) %>%
#' pc_dag_estimates()
#'
pc_dag_estimates <- function(.data, .alpha = 0.01) {
    number_samples <- nrow(.data)
    metabolite_names <- colnames(.data)

    pcalg::pc(
        suffStat = list(C = stats::cor(.data), n = number_samples),
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

#' Extract adjacency matrix from a DAG skeleton.
#'
#' Is generally a wrapper around calls to [igraph::get.adjacency()] and
#' [igraph::igraph.from.graphNEL()]. Transforms from a GraphNEL object in igraph.
#'
#' @param .dag_skeleton The PC DAG skeleton object.
#'
#' @return Outputs an adjacency matrix of the DAG skeleton.
#'
#' @examples
#' library(dplyr)
#' skeleton_estimate <- simulated_data %>%
#' select(contains("metabolite")) %>%
#' pc_skeleton_estimates()
#'
#' adjacency_matrix(skeleton_estimate)
#'
adjacency_matrix <- function(.dag_skeleton) {
    # TODO: Include a check here that it is a DAG skeleton..?
    igraph::get.adjacency(igraph::igraph.from.graphNEL(.dag_skeleton@graph))
}

#' Estimate order-independent PC-stable skeleton of a DAG.
#'
#' Uses the PC-algorithm and is mostly a wrapper around [pcalg::skeleton()].
#'
#' @param .data Input metabolic data.
#' @param .alpha Significance level threshold applied to each test.
#'
#' @return DAG skeleton object.
#'
#' @examples
#'
#' library(dplyr)
#' simulated_data %>%
#' select(contains("metabolite")) %>%
#' pc_skeleton_estimates()
#'
pc_skeleton_estimates <- function(.data, .alpha = 0.01) {
    number_samples <- nrow(.data)
    metabolite_names <- colnames(.data)

    # TODO: Confirm that this does this.
    pcalg::skeleton(
        suffStat = list(C = stats::cor(.data), n = number_samples),
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

#' Estimate Pearson's partial correlation coefficients.
#'
#' This function is a wrapper around [ppcor::pcor()] that extracts correlation
#' coefficient estimates, then adds the variable names to the column and row names.
#'
#' @param .data Input data of metabolic variables as matrix or data.frame.
#'
#' @return Outputs a matrix of partial correlation coefficients.
#'
#' @examples
#'
#' library(dplyr)
#' simulated_data %>%
#' select(contains("metabolite")) %>%
#' partial_corr_matrix()
#'
partial_corr_matrix <- function(.data) {
    # TODO: check if input data is gaussian
    pcor_matrix <- ppcor::pcor(.data)$estimate
    colnames(pcor_matrix) <- colnames(.data)
    rownames(pcor_matrix) <- colnames(.data)
    return(pcor_matrix)
}
