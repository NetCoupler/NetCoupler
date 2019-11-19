#test git workflow
#' Title
#'
#' @param .data
#' @param .alpha
#'
#' @return
#' @export
#'
#' @examples
nc_make_network1 <- function(.data,.alpha,.network_variables) {
  

  est.pcor.skel.DAG.adj<-function(dat,alpha_val,.network_variables){
    
    #dat: samples x metabolites data matrix
    
    library(ppcor)
    library(pcalg)
    library(igraph)
    pCor_mat<-ppcor::pcor(dat)$estimate                             #estimate Pearson's partial correlation coefficients
    colnames(pCor_mat)<-colnames(data.frame(dplyr::select(dat,noquote(.network_variables))))
    rownames(pCor_mat)<-colnames(data.frame(dplyr::select(dat,noquote(.network_variables))))
    
    number_observations<-nrow(dat)
    metabolite_names<-colnames(data.frame(dplyr::select(dat,noquote(.network_variables))))
    
    #estimate order-independent "PC-stable" skeleton of DAG using PC-algorithm
    skel_est<-pcalg::skeleton(suffStat = list(C=cor(dat),n = number_observations),
                              indepTest = gaussCItest,
                              labels = metabolite_names,
                              method = "stable",
                              alpha = alpha_val,
                              fixedGaps = NULL,
                              fixedEdges = NULL,
                              verbose = FALSE)
    #estimate order-independent DAG using PC-algorithm
    DAG_est<-pcalg::pc(suffStat = list(C=cor(dat), n = number_observations),
                       indepTest = gaussCItest,
                       labels = metabolite_names,
                       skel.method = "stable",
                       alpha = alpha_val,
                       fixedGaps = NULL,
                       fixedEdges = NULL,
                       verbose = FALSE,
                       maj.rule = FALSE,
                       solve.confl = FALSE)
    adj_matrix<-get.adjacency(igraph.from.graphNEL(skel_est@graph))          #return adjacency matrix of DAG skeleton
    
    w_adj_matrix<-adj_matrix*pCor_mat
    return(list(pCor_mat = pCor_mat, skel_est = skel_est, DAG_est = DAG_est, adj_matrix = adj_matrix, w_adj_matrix=w_adj_matrix))
    
  }
  metabolite_data<-dplyr::select(.data, noquote(.network_variables))
  ##################
  ##################
  # pcalg
  DAG<-est.pcor.skel.DAG.adj(dat =metabolite_data, alpha_val = .alpha, .network_variables = .network_variables)
  plot(DAG$skel_est)

  ##################
  # analyse igraph : let's deal with the clustering later
  #iDAG<- graph.adjacency(DAG$w_adj_matrix, weighted=TRUE)
 

  #cluster_DAG<-cluster_optimal(iDAG, weights =E(iDAG)$weight)
  
  #membership( cluster_DAG)
  #sizes( cluster_DAG)
  #plot( cluster_DAG,iDAG)
}  


