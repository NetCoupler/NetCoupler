
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
nc_make_network1 <- function(.data,.alpha) {
     #Standardize metabolites on exposure, remove all variables that are not metbolites
     exposureadjusted_simulated_metabolite_data<-.data %>%tibble::as_tibble() %>%purrr::modify_at(c(4:ncol(.data)), ~ residuals(lm(.x ~ as.matrix(dplyr::select(.data,exposure))))) %>%dplyr::select(-survival_time, -case_status, -exposure)

     #Estimate the skeleton (family of DAGs without specification of the direction of edges)
     est.pcor.skel.DAG.adj<-function(dat,alpha_val){
         #estimate Pearson's partial correlation coefficients
         pCor_mat<-ppcor::pcor(dat)$estimate
         colnames(pCor_mat)<-colnames(dat)
         rownames(pCor_mat)<-colnames(dat)

         number_observations<-nrow(dat)
         metabolite_names<-colnames(dat)
         #estimate order-independent "PC-stable" skeleton of DAG using PC-algorithm
         skel_est<-pcalg::skeleton(suffStat = list(C=cor(dat),n = number_observations),
                            indepTest = gaussCItest,
                            labels = metabolite_names, method = "stable", alpha = alpha_val, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE)
         #estimate order-independent DAG using PC-algorithm
         DAG_est<-pcalg::pc(suffStat = list(C=cor(dat),n = number_observations),
                     indepTest = gaussCItest, labels = metabolite_names, skel.method = "stable",
                     alpha = alpha_val, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, maj.rule = FALSE, solve.confl = FALSE)

         adj_matrix<-igraph::get.adjacency(igraph.from.graphNEL(skel_est@graph))          #return adjacency matrix of DAG skeleton

         w_adj_matrix<-adj_matrix*pCor_mat
         return(list(pCor_mat = pCor_mat, skel_est = skel_est, DAG_est = DAG_est, adj_matrix = adj_matrix, w_adj_matrix=w_adj_matrix))

     }

     DAG<-est.pcor.skel.DAG.adj(dat=exposureadjusted_simulated_metabolite_data,alpha_val=.alpha)


}


