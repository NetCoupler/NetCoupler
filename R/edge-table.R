#' Convert network graph to edge table.
#'
#' @description
#' \lifecycle{experimental}
#'
#' @param .edge_list Network graph from e.g. [nc_create_network()]
#'
#' @return A [tibble][tibble::tibble-package], with at least two columns:
#'
#' - `source_node`: The starting node (variable).
#' - `target_node`: The ending node (variable) that links with the source node.
#'
#' @export
#'
#' @seealso For examples of usage, see [nc_model_estimates()].
#'
as_edge_tbl <- function(.edge_list) {
    # TODO: Convert this to S3 method in case other network methods are different
    if (class(.edge_list) != "pcAlgo") {
        rlang::abort("The .edge_list is not from the pcalgo package. We currently do not have support for other network packages.")
    }
    .edge_list <- .edge_list@graph@edgeL
    nodes <- names(.edge_list)
    edge_table <- purrr::map_dfr(
        .edge_list,
        .single_edge_list_to_tbl,
        .id = "source_node",
        .nodes = nodes
    )
    return(edge_table)
}

.single_edge_list_to_tbl <- function(.edges, .nodes) {
    tibble(target_node = .nodes[.edges$edges])
}
