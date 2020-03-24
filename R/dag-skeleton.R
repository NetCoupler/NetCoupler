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
#' \dontrun{
#' library(dplyr)
#' simulated_data %>%
#' select(contains("metabolite")) %>%
#' NetCoupler:::pc_dag_estimates()
#' }
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

#' Estimate order-independent PC-stable skeleton of a DAG.
#'
#' Uses the PC-algorithm and is mostly a wrapper around [pcalg::skeleton()].
#'
#' @param .data Input metabolic data.
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

