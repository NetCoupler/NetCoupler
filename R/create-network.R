
#' Create an estimate of the metabolic network skeleton.
#'
#' Main NetCoupler network creator. Estimates the skeleton based on a family of
#' DAGs without specifying the direction of edges.
#'
#' @param .data Data of the metabolic variables.
#' @param .alpha The alpha level to set.
#'
#' @return Outputs a DAG skeleton.
#' @export
#'
#' @examples
#' library(dplyr)
#' simulated_data %>%
#' select(contains("metabolite")) %>%
#' nc_create_network()
nc_create_network <- function(.data, .alpha) {

    pcor_matrix <- partial_corr_matrix(.data)
    skel_est <- pc_skeleton_estimates(.data, .alpha)
    DAG_est <- pc_dag_estimates(.data, .alpha)
    adj_matrix <- adjacency_matrix(skel_est)
    # TODO: What does the w stand for?
    w_adj_matrix <- adj_matrix * pcor_matrix

    # TODO: Rename these to be more descriptive
    network_list <- list(
        pCor_mat = pcor_matrix,
        skel_est = skel_est,
        DAG_est = DAG_est,
        adj_matrix = adj_matrix,
        w_adj_matrix = w_adj_matrix,
        data = .data
    )
    return(network_list)
}
