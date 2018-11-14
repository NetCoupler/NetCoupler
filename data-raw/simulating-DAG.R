# A DAG
library(dagitty)
library(glasso)
library(magrittr)

graph <-
    dagitty("dag{
        a -> b[beta=.3]
        b -> c[beta=.2]
        c -> d[beta=.5]
        d -> e[beta=.3]
        e -> f[beta=.6]
        b [beta=.2]
        b -> k [beta=.3]
        k -> l [beta=.3]
        l -> m [beta=.5]
        m -> o [beta=.3]
        o -> p [beta=.5]
        a -> v [beta=.2]
        c -> v [beta=.5]
        d -> v [beta=.3]
        e -> m [beta=.5]
        k -> o [beta=.3]
        o -> p [beta=.65]
        }")

simulated_data <- simulateSEM(graph, N = 2500) %>%
    setNames(paste0("metabolite_", 1:length(.)))

devtools::use_data(simulated_data, overwrite = TRUE)


X_DAG<-est.pcor.skel.DAG.adj(dat=x,alpha_val=0.05)
plot(X_DAG$skel_est)


X_GLas<-glasso(as.matrix(cov(x)), rho=.2)
adjM_GLas<-X_GLas$wi
colnames(adjM_GLas)<-colnames(x)
rownames(adjM_GLas)<-colnames(x)
GLas_net<-graph_from_adjacency_matrix(adjM_GLas, mode = c("undirected"), weighted = TRUE, diag = FALSE,
                                     add.colnames = NULL, add.rownames = NA)
tkplot(GLas_net, vertex.size=20, vertex.color="white")


