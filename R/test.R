X_DAG<-est.pcor.skel.DAG.adj(dat=x,alpha_val=0.05)
plot(X_DAG$skel_est)


X_GLas<-glasso(as.matrix(cov(x)), rho=.2)
adjM_GLas<-X_GLas$wi
colnames(adjM_GLas)<-colnames(x)
rownames(adjM_GLas)<-colnames(x)
GLas_net<-graph_from_adjacency_matrix(adjM_GLas, mode = c("undirected"), weighted = TRUE, diag = FALSE,
                                      add.colnames = NULL, add.rownames = NA)
tkplot(GLas_net, vertex.size=20, vertex.color="white")
