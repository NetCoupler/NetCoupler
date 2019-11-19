
# Get connected components per exposure#
#blablabla
get.con.comp <- function(exposure_names, exposure_list, adjM_norename, met_group) {

  # exposure_names: list of investigated exposure names
  # exposure_list: list of summary statistics for all metabolites including round-number information
  # adjM_norename: not renamed adjacency matrix of metabolite DAG
  # met_group: name of investigated metabolite group

  NC_res <- list()
  all_CC <- list()

  # for each exposure in exposure_list do:
  for (k in exposure_names) {
    exposure_list[[k]]$Metabolite <- as.character(exposure_list[[k]]$Metabolite) # rename true metabolite names as characters

    # delete edges from adjacency matrix whenever node, which is not associated with exposure, is part of node-pair:
    adjM2CoCo <- adjM_norename
    nonAssoc_nodes <- c(unlist(exposure_list[[k]] %>% dplyr::filter(Assoc == 0) %>% dplyr::select(Metabolite))) # nodes not associated with exposure
    Assoc_nodes <- c(unlist(exposure_list[[k]] %>% dplyr::filter(Assoc != 0) %>% dplyr::select(Metabolite))) # nodes associated with exposure
    adjM2CoCo[nonAssoc_nodes, ] <- 0 # erase edges containing not associated nodes
    adjM2CoCo[, nonAssoc_nodes] <- 0 # erase edges containing not associated nodes

    # extract connected components:
    c <- clusters(graph_from_adjacency_matrix((adjM2CoCo))) # information about which associated nodes are connected to which nodes
    n_cc <- max(c$membership) # number of different clusters (why not c$no?)

    connect_comp <- list()
    nam <- c()
    j <- c(0) # ??????????

    # for each connected component (CC) assign logical indicator to nodes: TRUE if in CC, FALSE otherwise:
    for (i in 1:n_cc) { # loop over all different clusters
      idx <- c$membership == i # logical indicator if metabolite is in cluster i

      if (sum(idx) > 1) { # if cluster contains more than one metabolite then these metabolites are connected components (?)
        j <- j + 1 # ?????????????????????

        connect_comp[j] <- list(assign(paste0(k, "CC", sep = "_", j), idx)) # add output to list
        nam <- c(nam, paste0(k, "CC", sep = "_", j)) # define name for list object: exposureCC_j
      }

      names(connect_comp) <- nam # assign names to listed objects
      CC <- data.frame(connect_comp[1]) %>% draw.rownames() %>% dplyr::select(Metabolite) # get metabolite names

      for (i in 1:length(connect_comp)) { # loop over list objects
        CC_temp <- data.frame(connect_comp[i]) %>% draw.rownames()
        CC_temp$Metabolite <- as.character(CC_temp$Metabolite)
        CC <- dplyr::left_join(CC, CC_temp, by = "Metabolite")
      }

      # combine with multi-model output (stored in exposure_list):
      # assign(paste0("CC",sep=".",k,sep=".",met_group),dplyr::left_join(exposure_list[[k]],CC,by="Metabolite"))           #?????????????????????
      NC_res[k] <- list(assign(paste0("CC", sep = ".", k, sep = ".", met_group), dplyr::left_join(exposure_list[[k]], CC, by = "Metabolite")))
    }

    # store all CCs in list:
    # assign(paste0("CC",sep="_",k,sep="_",met_group),connect_comp)                          #???????????????????
    all_CC <- append(all_CC, assign(paste0("CC", sep = "_", k, sep = "_", met_group), connect_comp)) # list of connected components for all outcomes
  }

  return(list(NC_res = NC_res, all_CC = all_CC))
}

