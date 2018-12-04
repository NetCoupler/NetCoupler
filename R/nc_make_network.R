
#' Title
#'
#' @param .data
#' @param .alpha
#'
#' @return
#' @export
#'
#' @examples
#'  test_net<-nc_make_network1(simulated_data,.001)
#'
nc_make_network1 <- function(.data,.alpha,.network_variables,.filter_variable,.filter_value,.adjustment_variable) {

    #Standardize metabolites on exposure, remove all variables that are not metbolites
    filtered_data<- dplyr::filter(.data, UQ(rlang::sym(.filter_variable))== .filter_value)%>%dplyr::filter_at(vars(noquote(.adjustment_variable)), all_vars(!is.na(.)))

    exposureadjusted_simulated_metabolite_data<-filtered_data%>%
        tibble::as_tibble() %>%
        purrr::modify_at(noquote(.network_variables), ~ residuals(lm(.x ~ as.matrix(dplyr::select(filtered_data,noquote(.adjustment_variable)))))) %>%
        dplyr::select(noquote(.network_variables))


    #Estimate the skeleton (family of DAGs without specification of the direction of edges)
    est.pcor.skel.DAG.adj<-function(.dat,.alpha_val, .network_variables){
        #estimate Pearson's partial correlation coefficients
        pCor_mat<-ppcor::pcor(.dat)$estimate
        colnames(pCor_mat)<-colnames(data.frame(dplyr::select(.dat,noquote(.network_variables))))
        rownames(pCor_mat)<-colnames(data.frame(dplyr::select(.dat,noquote(.network_variables))))

        number_observations<-nrow(.dat)
        metabolite_names<-colnames(data.frame(dplyr::select(.dat,noquote(.network_variables))))

        #estimate order-independent "PC-stable" skeleton of DAG using PC-algorithm
        skel_est<-pcalg::skeleton(suffStat = list(C=cor(.dat),n = number_observations),
                                  indepTest = gaussCItest,
                                  labels = metabolite_names,
                                  method = "stable",
                                  alpha = .alpha_val,
                                  fixedGaps = NULL,
                                  fixedEdges = NULL,
                                  verbose = FALSE)
        #estimate order-independent DAG using PC-algorithm
        DAG_est<-pcalg::pc(suffStat = list(C=cor(.dat), n = number_observations),
                           indepTest = gaussCItest,
                           labels = metabolite_names,
                           skel.method = "stable",
                           alpha = .alpha_val,
                           fixedGaps = NULL,
                           fixedEdges = NULL,
                           verbose = FALSE,
                           maj.rule = FALSE,
                           solve.confl = FALSE)

        #return adjacency matrix of DAG skeleton
        adj_matrix<-igraph::get.adjacency(igraph.from.graphNEL(skel_est@graph))

        w_adj_matrix<-adj_matrix*pCor_mat
        return(list(pCor_mat = pCor_mat, skel_est = skel_est, DAG_est = DAG_est, adj_matrix = adj_matrix, w_adj_matrix=w_adj_matrix, data=exposureadjusted_simulated_metabolite_data))

    }

    DAG<-est.pcor.skel.DAG.adj(.dat = exposureadjusted_simulated_metabolite_data, .alpha_val = .alpha, .network_variables = .network_variables)

}


