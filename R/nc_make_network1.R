
#' Title
#'
#' @param .data
#' @param .alpha
#'
#' @return
#' @export
#'
#' @examples
nc_make_network1 <- function(.data,.alpha) {
  
  
  est.pcor.skel.DAG.adj<-function(dat,alpha_val){
    
    #dat: samples x metabolites data matrix
    
    library(ppcor)
    library(pcalg)
    library(igraph)
    pCor_mat<-ppcor::pcor(dat)$estimate                             #estimate Pearson's partial correlation coefficients
    colnames(pCor_mat)<-colnames(dat)
    rownames(pCor_mat)<-colnames(dat)
    
    n_samples<-nrow(dat)                                            #number of samples
    V_met<-colnames(dat)                                            #labels equal to node names i.e. metabolites
    
    skel_est<-skeleton(suffStat = list(C=cor(dat),n = n_samples),     #estimate order-independent "PC-stable" skeleton of DAG using PC-algorithm
                       indepTest = gaussCItest,                        #test conditional independence of Gaussians via Fisher's Z
                       labels = V_met, method = "stable", alpha = alpha_val, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE)
    
    DAG_est<-pc(suffStat = list(C=cor(dat),n = n_samples),           #estimate equivalence class of DAG using PC-algorithm
                indepTest = gaussCItest, labels = V_met, skel.method = "stable",          #order-independent skeleton
                alpha = alpha_val, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, maj.rule = FALSE, solve.confl = FALSE)
    
    adj_matrix<-get.adjacency(igraph.from.graphNEL(skel_est@graph))          #return adjacency matrix of DAG skeleton
    
    w_adj_matrix<-adj_matrix*pCor_mat
    return(list(pCor_mat = pCor_mat, skel_est = skel_est, DAG_est = DAG_est, adj_matrix = adj_matrix, w_adj_matrix=w_adj_matrix))
    
  }
  ##################
  ##################
  # pcalg
  DAG<-est.pcor.skel.DAG.adj(dat=.data,alpha_val=.alpha)
  plot(DAG$skel_est)

  ##################
  # analyse igraph
  iDAG<- graph.adjacency(DAG$w_adj_matrix, weighted=TRUE)
 

  cluster_DAG<-cluster_optimal(iDAG, weights =E(iDAG)$weight)
  
  membership( cluster_DAG)
  sizes( cluster_DAG)
  plot( cluster_DAG,iDAG)
}  


