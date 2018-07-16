# PC_skel <- skeleton(
#   suffStat = list(C = cor(PC_data_SC), n = n),
#   indepTest = gaussCItest, labels = V, method = "stable", # indep.test: partial correlations
#   alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE
# )
#
# PC_DAG <- pc(
#   suffStat = list(C = cor(PC_data_SC), n = n),
#   indepTest = gaussCItest, labels = V, skel.method = "stable", # indep.test: partial correlations
#   alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, maj.rule = F, solve.confl = F
# )
#
# # get adjacency matrix (transform GraphNEL in Igraph)
# adjM_PC <- get.adjacency(graphNEL2igraph(PC_skel@graph))
#
# NW_PC_EPIC <- PC_DAG@graph
#
# # GET CONNECTED COMPONENTS PER EXPOSURE
# # Connected component(CC):
# # CC Cluster(>=2) of metabolites where
# # Met~Exp in non-metabolite adjusted model
# NC.PC_res <- list()
# nam <- c()
# ALL_CC <- list()
# # For each exposure do
# for (k in Exp)
# {
#   print(k)
#   EXP_PC[[k]]$Metabolite <- as.character(EXP_PC[[k]]$Metabolite)
#   # delete edges from adjacency matrix whenever non-exposure associated node is part of the node-pair
#   adjM2CoCo <- adjM_PC
#   nonAssociated_Nodes <- c(unlist(EXP_PC[[k]] %>% dplyr::filter(Assoc == 0) %>% dplyr::select(Metabolite)))
#   Associated_Nodes <- c(unlist(EXP_PC[[k]] %>% dplyr::filter(Assoc != 0) %>% dplyr::select(Metabolite)))
#   adjM2CoCo[nonAssociated_Nodes, ] <- 0
#   adjM2CoCo[, nonAssociated_Nodes] <- 0
#
#   # extract connected components
#   c <- clusters(graph_from_adjacency_matrix(adjM2CoCo))
#   c$membership
#   n_cc <- max(c$membership)
#   CC_PC <- list()
#   nam <- c()
#   j <- c(0)
#   # For each connected component (CC) assign values to Nodes
#   # "True" if in CC, "False" otherwise
#   for (i in 1:n_cc)
#   {
#     # who's in there?
#     idx <- c$membership == i
#
#     if (sum(idx) > 1) {
#       print(i)
#       j <- j + 1
#       # print(j)
#       # no singleton
#       # do something with this one
#       # checking algorithm
#       # C = B[idx,idx]
#       CC_PC[j] <- list(assign(paste0(k, "CC", sep = "_", j), idx)) # adment Output to  list
#       nam <- c(nam, paste0(k, "CC", sep = "_", j)) # define name for list object
#       # cat(idx)
#       # cat('\n')
#       # cat('=====\n')
#     }
#     names(CC_PC) <- nam # assign names to listed objects
#     # transform into dataframe
#     CC <- data.frame(CC_PC[1]) %>% draw_rownames() %>% dplyr::select(Metabolite)
#     for (i in 1:length(CC_PC))
#     {
#       CC_temp <- data.frame()
#       CC_temp <- data.frame(CC_PC[i]) %>% draw_rownames()
#       CC_temp$Metabolite <- as.character(CC_temp$Metabolite)
#       CC <- dplyr::left_join(CC, CC_temp, by = "Metabolite")
#     }
#     # Combine with multi_model output
#     assign(paste0("CC", sep = ".", k, sep = ".", "PC"), dplyr::left_join(EXP_PC[[k]], CC, by = "Metabolite"))
#     NC.PC_res[k] <- list(assign(paste0("CC", sep = ".", k, sep = ".", "PC"), dplyr::left_join(EXP_PC[[k]], CC, by = "Metabolite")))
#
#     rm(CC)
#   }
#   # Keep a list of CC as well
#   assign(paste0("CC", sep = "_", k, sep = "_", "PC"), CC_PC)
#   ALL_CC <- append(ALL_CC, assign(paste0("CC", sep = "_", k, sep = "_", "PC"), CC_PC))
#   rm(CC_PC)
# }
# #
# # GET CONNECTED COMPONENTS PER Outcome
# # Connected component(CC):
# # CC Cluster(>=2) of metabolites where
# # Met~Exp in non-metabolite adjusted model
#
# # For each Outcome do
# for (k in Out)
# {
#   # delete edges from adjacency matrix whenever non-exposure associated node is part of the node-pair
#   adjM2CoCo <- adjM_PC
#   nonAssociated_Nodes <- c(unlist(OUT_PC[[k]] %>% filter(Assoc == 0) %>% dplyr::select(Metabolite)))
#   adjM2CoCo[nonAssociated_Nodes, ] <- 0
#   adjM2CoCo[, nonAssociated_Nodes] <- 0
#
#   # extract connected components
#   c <- clusters(graph_from_adjacency_matrix(adjM2CoCo))
#   c$membership
#   n_cc <- max(c$membership)
#
#   CC_PC <- list()
#   nam <- c()
#   j <- c(0)
#   # For each connected component (CC) assign values to Nodes
#   # "True" if in CC, "False" otherwise
#   for (i in 1:n_cc)
#   {
#     # who's in there?
#     idx <- c$membership == i
#     # print(idx)
#
#     if (sum(idx) > 1) {
#       j <- j + 1
#       # print(j)
#       # no singleton
#       # do something with this one
#       # checking algorithm
#       # C = B[idx,idx]
#       CC_PC[j] <- list(assign(paste0(k, "CC", sep = "_", j), idx)) # adment Output to  list
#       nam <- c(nam, paste0(k, "CC", sep = "_", j)) # define name for list object
#       # cat(idx)
#       # cat('\n')
#       # cat('=====\n')
#     }
#     names(CC_PC) <- nam # assign names to listed objects
#     # print(nam)
#     # transform into dataframe
#     CC <- data.frame(CC_PC[1]) %>% draw_rownames_out() %>% dplyr::select(Metabolite)
#     CC$Metabolite <- as.character(CC$Metabolite)
#
#     for (i in 1:length(CC_PC))
#     {
#       CC_temp <- data.frame()
#       CC_temp <- data.frame(CC_PC[i]) %>% draw_rownames_out()
#       CC_temp$Metabolite <- as.character(CC_temp$Metabolite)
#       CC <- left_join(CC, CC_temp, by = "Metabolite")
#     }
#     # Combine with multi_model output
#     assign(paste0("CC", sep = ".", k, sep = ".", "PC"), left_join(OUT_PC[[k]], CC, by = "Metabolite"))
#     NC.PC_res[k] <- list(assign(paste0("CC", sep = ".", k, sep = ".", "PC"), left_join(OUT_PC[[k]], CC, by = "Metabolite")))
#     rm(CC)
#   }
#   # Keep a list of CC as well
#   assign(paste0("CC", sep = "_", k, sep = "_", "PC"), CC_PC)
#   ALL_CC <- append(ALL_CC, assign(paste0("CC", sep = "_", k, sep = "_", "PC"), CC_PC))
#   rm(CC_PC)
# }
