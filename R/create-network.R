
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
    assert_is_number(.alpha)

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
