# NETcoupler general R-commands#

# longer analysis version#

rm(list = ls())

# Rename feature names in order to avoid clash with glm#

rename.met <- function(dat) {

  # dat: samples x metabolites data matrix

  Ll <- paste("NM", c(1:dim(dat)[2]), sep = "") # generate shorter metabolite names

  names_mapping <- cbind(colnames(dat), Ll) # mapping of old and new metabolite names
  colnames(names_mapping) <- c("Metabolite", "Outcome")

  data_renamed <- dat
  colnames(data_renamed) <- Ll # is character!

  return(list(data_renamed = data_renamed, names_mapping = names_mapping))
}

# Obtain partial correlation matrix, DAG skeleton, DAG and adjacency matrix for DAG skeleton#

est.pcor.skel.DAG.adj <- function(dat, alpha_val) {

  # dat: samples x metabolites data matrix

  # check if input data is gaussian

  pCor_mat <- ppcor::pcor(dat)$estimate # estimate Pearson's partial correlation coefficients
  colnames(pCor_mat) <- colnames(dat)
  rownames(pCor_mat) <- colnames(dat)

  n_samples <- nrow(dat) # number of samples
  V_met <- colnames(dat) # labels equal to node names i.e. metabolites

  skel_est <- skeleton(
    suffStat = list(C = cor(dat), n = n_samples), # estimate order-independent "PC-stable" skeleton of DAG using PC-algorithm
    indepTest = gaussCItest, # test conditional independence of Gaussians via Fisher's Z
    labels = V_met, method = "stable", alpha = alpha_val, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE
  )

  DAG_est <- pc(
    suffStat = list(C = cor(dat), n = n_samples), # estimate equivalence class of DAG using PC-algorithm
    indepTest = gaussCItest, labels = V_met, skel.method = "stable", # order-independent skeleton
    alpha = alpha_val, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, maj.rule = FALSE, solve.confl = FALSE
  )

  adj_matrix <- get.adjacency(graphNEL2igraph(skel_est@graph)) # return adjacency matrix of DAG skeleton

  return(list(pCor_mat = pCor_mat, skel_est = skel_est, DAG_est = DAG_est, adj_matrix = adj_matrix))
}

# Net Coupler IN function: exposure -> metabolome#

# return: for each metabolite all possible models with different adjacency subsets and complete list of exposures and covariates

net.coupler.in <- function(graph_skel, dat, dat_compl, exp_dat, DE, glmulti_method) {

  # graph_skel: estimated DAG skeleton of samples x metabolites data matrix
  # dat: renamed samples x metabolites data matrix
  # exp_dat: exposure/phenotype data
  # DE: indicator if direct effects were already identified
  # glmulti_method: choose method for glmulti-function

  cat("*************************************************************************************************** \n")
  cat("This algorithm estimates direct effects of a predefined exposure on each network-variable for all   \n")
  cat("causal models that agree with the input-network: models are adjusted for all possible combinations  \n")
  cat("of direct neighbors (==variables in the adjacency set) -> Output is a multiset of possible effects \n")
  cat("*************************************************************************************************** \n")

  exp_names <- colnames(exp_dat) # names of exposures/covariates

  model_details_all <- list(NULL) # prepare empty list to store output
  node_names <- colnames(dat) # create vector "nodes" with node names as strings

  for (i in seq(along = node_names)) { # create empty list with slots for each network-variable, i.e. metabolite

    model_details <- list(NULL)
    model_details_all[[i]] <- model_details
  }
  names(model_details_all) <- node_names

  for (i in node_names) { # net.coupler.in loop over all metabolite nodes

    # prepare separate datasets for each metabolite:#

    met_outcome <- i

    # check if met_outcome is character string:
    if (is.character(met_outcome) == FALSE) stop("'met_outcome' is not a character string as required")

    # select data on met_outcome within samples x metabolites data matrix, store in "met_out_data":
    met_out_data <- dat[, met_outcome]

    # check if met_out_data is numeric:
    if (is.numeric(met_out_data) == FALSE) stop("'met_out_data' is not a numeric vector as required")

    # create vector with integers indicating adjacent variables, i.e. metabolites in skeleton:
    edge_list <- slot(graph_skel@graph, "edgeL") # extract edge-list from skeleton
    adjset <- c(edge_list[[met_outcome]]) # extract adjacency set for selected node/metabolite
    adjset <- c(adjset[[1]]) # extract indices of adjacency set

    if ((length(adjset) == 0) == FALSE) { # check if selected metabolite has adjacency set

      if (length(adjset) >= 15) stop(paste("adjacency set for metabolite", i, " larger than 15, please lower adjacency set", sep = ""))

      print(i)
      print("neighbors")

      # select data on adjacency set, store in "adjset_data":
      adjset_data <- subset(dat_compl, select = c(adjset))

      # check if adjset_data is numeric:
      if (is.numeric(adjset_data[1, 1]) == FALSE) stop("'adjset_data' is not a matrix of numeric variables as required")

      # create vector with names of adjacency set as strings:
      adjset_names <- colnames(adjset_data)

      # for loop with already identified direct effects:
      if (is.null(DE) == FALSE) {
        adjset_names <- setdiff(adjset_names, as.vector(DE))
        adjset_data <- adjset_data[, c(adjset_names)]
      }

      # check if adjset_names is vector of character strings:
      if (is.character(adjset_names) == FALSE) stop("'adjset_names' is not a vector of character strings as required")

      # combine data on exposure (including already identified direct effect metabolites), metabolite outcome, and adjacency set, store as dataframe ("modeldata"), input data for glm:
      modeldata <- data.frame(cbind(met_out_data, exp_dat, adjset_data))

      # estimate multimodel coefficients#

      # modify fitting function of glm so that it always includes one predefined set of variables (here: exposure + covariates), while subsetting adjacent metabolite set:
      glm.redefined <- function(formula, data, always="", ...) {
        glm(as.formula(paste(deparse(formula), always)), data = data, ...)
      }

      # fit all possible causal models using glmulti: iterate over all possible combinations of adjacent metabolites, fixed exposure + covariates set

      start_time_glmulti <- Sys.time()
      glmulti_obj <- glmulti(
        y = "met_out_data", # response variable to be predicted: selected metabolite
        xr = c(adjset_names[1:length(adjset_names)]), # predictor variables to be tested in iteration
        data = modeldata, # dataframe containing data
        level = 1, # only main effects (terms of order 1) are used to build the candidate set
        confsetsize = 2^(length(adjset_names)), # number of models to be looked for, dependent on size of adjacency set
        fitfunction = glm.redefined, # fitting function to be used: here: glm-function with fixed predictor set "always"
        always = paste0("+", paste0(exp_names, collapse = "+ ")), # defines "always"-set: complete set of exposures and covariates
        includeobjects = TRUE, # fitted models are included as objects
        intercept = TRUE, # no intercept included in candidate models
        plotty = FALSE, # no progress of the IC profile is plotted
        report = FALSE, # report about progress at run time is given
        method = glmulti_method # choose glmulti method
      )
      end_time_glmulti <- Sys.time()
      print(paste("run time for glmulti:", as.numeric(end_time_glmulti - start_time_glmulti, units = "secs")))

      # output summaries:

      glmulti_obj_objects <- glmulti_obj@objects # collect glmulti objects, i.e. all fitted models (?)

      nbmds <- c(1:glmulti_obj@nbmods) # number of fitted models of glmulti-function

      model_details <- list(NULL)

      for (j in seq(along = nbmds)) { # generate list for model details
        model_details[[j]] <- list(NULL)
      }
      names(model_details) <- paste("Model", nbmds, sep = "_")

      for (j in seq(along = nbmds)) { # collect all model details for different adjacency subsets for specific metabolite
        model_summary <- glmulti_obj_objects[[j]]
        model_details[[j]] <- list(
          Model = paste("Model", j, "of", length(nbmds)),
          Model_summary = model_summary
        )
      }

      model_details_all[[met_outcome]] <- list(
        Model_summaries = model_details, # collect all model details for different adjacency subsets for specific metabolite in complete model list for all metabolites
        Number_of_Models = length(nbmds),
        Adj_set = paste(adjset_names, collapse = ", "),
        Outcome_metabolite = met_outcome
      )
    } else {
      print(i)
      print("no neighbors")

      # combine data on exposure (including already identified direct effect metabolites) and metabolite outcome, store as dataframe ("modeldata"), input data for glm:
      modeldata <- data.frame(cbind(met_out_data, exp_dat))

      # estimate multimodel coefficients#

      # modify fitting function of glm so that it always includes one predefined set of variables (here: exposure + covariates), while subsetting adjacent metabolite set:
      glm.redefined <- function(formula, data, always="", ...) {
        glm(as.formula(paste(deparse(formula), always)), data = data, ...)
      }

      # fit all possible causal models using glmulti: iterate over all possible combinations of adjacent metabolites, fixed exposure + covariates set

      start_time_glmulti <- Sys.time()
      glmulti_obj <- glmulti(
        y = "met_out_data", # response variable to be predicted: selected metabolite
        xr = exp_names[1], # predictor variables to be tested in iteration: here use first covariate as dummy variable
        data = modeldata, # dataframe containing data
        level = 1, # only main effects (terms of order 1) are used to build the candidate set
        # confsetsize = 1,                                      #number of models to be looked for, dependent on size of adjacency set
        fitfunction = glm.redefined, # fitting function to be used: here: glm-function with fixed predictor set "always"
        always = paste0("+", paste0(exp_names[-1], collapse = "+ ")), # defines "always"-set: complete set of exposures and covariates
        includeobjects = TRUE, # fitted models are included as objects
        intercept = TRUE, # no intercept included in candidate models
        plotty = FALSE, # no progress of the IC profile is plotted
        report = FALSE, # report about progress at run time is given
        method = glmulti_method # choose glmulti method
      )
      end_time_glmulti <- Sys.time()
      print(paste("run time for glmulti:", as.numeric(end_time_glmulti - start_time_glmulti, units = "secs")))

      # output summaries:

      glmulti_obj_objects <- glmulti_obj@objects # collect glmulti objects, i.e. all fitted models (?)

      nbmds <- c(1:glmulti_obj@nbmods) # number of fitted models of glmulti-function

      model_details <- list(NULL)

      for (j in seq(along = nbmds)) { # generate list for model details
        model_details[[j]] <- list(NULL)
      }
      names(model_details) <- paste("Model", nbmds, sep = "_")

      for (j in seq(along = nbmds)) { # collect all model details for different adjacency subsets for specific metabolite
        model_summary <- glmulti_obj_objects[[j]]
        model_details[[j]] <- list(
          Model = paste("Model", j, "of", length(nbmds)),
          Model_summary = model_summary
        )
      }

      model_details[[1]] <- NULL

      model_details_all[[met_outcome]] <- list(
        Model_summaries = model_details, # collect all model details for different adjacency subsets for specific metabolite in complete model list for all metabolites
        Number_of_Models = 1,
        Adj_set = "no adjacent metabolites",
        Outcome_metabolite = met_outcome
      )
    }
  }

  return(model_details_all)
}

# Extract exposure coefficients per outcome, i.e. metabolite#

getExp.coef.peroutcome <- function(object, outcome, exposure) {

  # object: output of net.coupler
  # outcome: specific outcome, i.e. metabolite, to evaluate
  # exposure: name of considered exposure variable

  # create vector containing integers from 1 to number of outcome-specific models:
  nbm <- c(1:length(object[[outcome]]$Model_summaries))
  Nbmds <- max(nbm)

  # define empty dataframes for output:
  mm_coef_temp1 <- structure(list(
    Model = as.character(),
    Nbmds = as.numeric(),
    Outcome = as.character(),
    Exposure = as.character(),
    Covariables = as.character(),
    Estimate = as.numeric(),
    SE = as.numeric(),
    tval = as.numeric(),
    P = as.numeric()
  ),
  class = "data.frame"
  )

  mm_coef_temp2 <- structure(list(
    Model = as.character(),
    Nbmds = as.numeric(),
    Outcome = as.character(),
    Exposure = as.character(),
    Covariables = as.character(),
    Estimate = as.numeric(),
    SE = as.numeric(),
    tval = as.numeric(),
    P = as.numeric()
  ),
  class = "data.frame"
  )

  # loop along number of outcome(metabolite)-specific models:
  for (i in seq(along = nbm)) {
    # get exposure-effect estimates (beta coefficient, SE, t-value, p-value) from single model and write into dataframe
    SUM <- summary(object[[outcome]]$Model_summaries[[i]]$Model_summary)

    Cov <- rownames(SUM$coefficients)

    # Cov<-dplyr::filter(data.frame(row.names(SUM$coefficients)),row.names(SUM$coefficients)!=exposure)                     #record specific covariates of this model

    mm_coef_temp2 <- data.frame(
      Model = as.character(object[[outcome]]$Model_summaries[[i]]$Model),
      Nbmds = Nbmds,
      Outcome = as.character(object[[outcome]]$Outcome_metabolite),
      Exposure = exposure,
      Covariables = as.character(paste(Cov, collapse = ", ")),
      Estimate = as.numeric(SUM$coefficients[paste(exposure), 1]),
      SE = as.numeric(SUM$coefficients[paste(exposure), 2]),
      tval = as.numeric(SUM$coefficients[paste(exposure), 3]),
      P = as.numeric(SUM$coefficients[paste(exposure), 4])
    )

    # bind information to an outcome-specific dataframe:
    mm_coef_temp1 <- bind_rows(mm_coef_temp1, mm_coef_temp2)
  }

  return(mm_coef_temp1)
}

# Extract exposure coefficients for specific exposure and outcomes, i.e. metabolites (?)#

getExp.coef <- function(object, outcome, exposure) {

  # object: output of net.coupler
  # outcome: list of metabolites to evaluate, e.g. all
  # exposure: name of considered exposure variable

  cat("*************************************************************************************************** \n")
  cat("This function produces a table of effect estimates of the exposure (betas, SE, t-value, p-value) on \n")
  cat("a selected (set of) outcome(s) [specified as character-vector of variable-names] for all possible   \n")
  cat("causal models based on conditional independence criteria encoded in the input-network => MULTISET OF\n")
  cat("POSSIBLE EFFECTS PER EXPOSURE-OUTCOME PAIRS: NUMBER OF POSSIBLE EFFECTS ~ SIZE OF THE ADJECENCY SET!\n")
  cat("***************************************************************************************************")

  # define empty dataframes for output:

  mm_coef <- structure(list(
    Model = as.character(),
    Nbmds = as.numeric(),
    Outcome = as.character(),
    Exposure = as.character(),
    Covariables = as.character(),
    Estimate = as.numeric(),
    SE = as.numeric(),
    tval = as.numeric(),
    P = as.numeric()
  ),
  class = "data.frame"
  )

  mm_coef_temp <- structure(list(
    Model = as.character(),
    Nbmds = as.numeric(),
    Outcome = as.character(),
    Exposure = as.character(),
    Covariables = as.character(),
    Estimate = as.numeric(),
    SE = as.numeric(),
    tval = as.numeric(),
    P = as.numeric()
  ),
  class = "data.frame"
  )

  # for each exposure-outcome pair: get exposure coefficients from all possible models and write to table:

  for (j in outcome) {
    mm_coef_temp <- data.frame(getExp.coef.peroutcome(object = object, outcome = j, exposure = exposure))
    mm_coef <- bind_rows(mm_coef, mm_coef_temp)
  }

  return(mm_coef)
}

# Summary statistics and multiple testing adjustment#

mult.stat <- function(sum_netin, MinMod, adjust_method, rule_1 = 11, rule_1_cut = 0.1, rule_2 = 5, rule_2_cut = 0.05, rule_3 = 4, rule_4 = 3, ass_rule1=9, round_number) {

  # sum_netin: output of getExp.coef after adding original metabolite names
  # MinMod: covariates of minimal model without any adjacency set
  # adjust_method: select adjustment method(s), e.g. "fdr" or c("fdr","bonferroni"), see p.adjust.methods for available options
  # rule= "set DE=1 or =2 if ..."
  # rule_1: determines first argument of rule (e.g. adjusted p-value of marginal model < rule_1_cut)
  # rule_2: determines second argument of rule (e.g. upperP < rule_2_cut)
  # rule_3: determines first variable of third argument of rule (e.g. highEst)
  # rule_4: determines second variable of third argument of rule (e.g. lowEst)
  # possible arguments for rule_...: adjusted p-values of marginal model=11=rule_1
  #                                 largest p-value across all adjacency models=upperP=5=rule_2
  #                                 highest estimated beta across all adjacency models=highEst=4=rule_3
  #                                 lowest estimated beta across all adjacency models=lowEst=3=rule_4
  # default rule: if(adjusted p< 0.1 & upperP<0.05 & (sign(highEst)==sign(lowEst)) & lowEst>0){DE=1}
  #              if(adjusted p< 0.1 & upperP<0.05 & (sign(highEst)==sign(lowEst)) & lowEst<0){DE=2}else{DE=0}
  # ass_rule1: determines first variable of second argument of association rule (e.g. bestGuess=9)
  # default association rule: if(adjusted p<0.1 & bestGuess>0){Assoc=1}
  #                          if(adjusted p<0.1 & bestGuess<0){Assoc=2}else{Assoc=0}
  # round_number: actual round number in ambiguous metabolites loop

  # return: for each metabolite: average beta, min beta, max beta, highest p-value, lowest p-value across all models and betas and p-values for minimal model (marg):
  sum_netin_sum <- dplyr::left_join(sum_netin %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgEst = mean(Estimate), lowEst = min(Estimate), highEst = max(Estimate), upperP = max(P), lowerP = min(P), Nbmds = mean(Nbmds)),
    sum_netin %>% dplyr::filter(Covariables == MinMod) %>% dplyr::select(Outcome = Outcome, Metabolite = Metabolite, bestGuess = Estimate, marg_P = P),
    by = "Outcome"
  )

  # multiple testing correction:
  p_adjust_M <- p.adjust.methods[p.adjust.methods %in% adjust_method] # select multiple testing correction method(s)
  p_adj <- sapply(p_adjust_M, function(meth) {
    p.adjust(sum_netin_sum$marg_P, meth)
  }) # calculate adjusted p-values for marg_P

  sum_netin_sum_adj <- bind_cols(sum_netin_sum, data.frame(p_adj))
  colnames(sum_netin_sum_adj)[11] <- adjust_method

  # is there an effect of specified exposure on specific metabolite?: determine effect-indicator (DE!=0: direct effect of exposure on outcome, DE=0: ambiguous if adjusted marg_p<0.1):
  sum_netin_sum_adj <- mutate(sum_netin_sum_adj,
    DE = derivedFactor(
      "1" = (sum_netin_sum_adj[, rule_1] < rule_1_cut & sum_netin_sum_adj[, rule_2] < rule_2_cut & (sign(sum_netin_sum_adj[, rule_3]) == sign(sum_netin_sum_adj[, rule_4])) & sum_netin_sum_adj[, rule_4] > 0),
      "2" = (sum_netin_sum_adj[, rule_1] < rule_1_cut & sum_netin_sum_adj[, rule_2] < rule_2_cut & (sign(sum_netin_sum_adj[, rule_3]) == sign(sum_netin_sum_adj[, rule_4])) & sum_netin_sum_adj[, rule_4] < 0),
      .method = "first", .default = 0
    )
  )

  # is there any effect of specified exposure on specific metabolite: determine association-indicator (Assoc!=0: there is an effect (no differenciation between direct or indirect), Assoc=0: there is no effect)
  sum_netin_sum_adj <- mutate(sum_netin_sum_adj,
    Assoc = derivedFactor(
      "1" = (sum_netin_sum_adj[, rule_1] < rule_1_cut & sum_netin_sum_adj[, ass_rule1] > 0),
      "2" = (sum_netin_sum_adj[, rule_1] < rule_1_cut & sum_netin_sum_adj[, ass_rule1] < 0),
      .method = "first", .default = 0
    )
  )

  # collect ambiguous effects:
  amb <- sum_netin_sum_adj %>% filter(Assoc != 0 & DE == 0)

  # collect direct effects:
  direct <- sum_netin_sum_adj %>% filter(Assoc != 0 & DE != 0)

  # collect summary statistics including round information
  sum_netin_sum_adj_FV <- sum_netin_sum_adj %>% dplyr::mutate(round = as.numeric(round_number), Assoc_FV = as.character(Assoc), DE_FV = as.character(DE))

  return(list(sum_netin_sum_adj = sum_netin_sum_adj, amb = amb, direct = direct, sum_netin_sum_adj_FV = sum_netin_sum_adj_FV))
}

# Get rownames as variable in dplyr#

draw.rownames <- function(dat) {
  dat %>% do(mutate(., Metabolite = rownames(.)))
}

draw.rownames.out <- function(dat) {
  dat %>% do(mutate(., Metabolite = as.factor(rownames(.))))
}

# Get connected components per exposure#

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

# Ambiguous metabolites loop#

amb.met.loop <- function(exp_dat, graph_skel, dat, dat_compl, DE, glmulti_method, exposure, met_map, adjust_method, round_number, adjM_norename, met_group, CC_ind) {

  # exp_dat: exposure/phenotype data
  # graph_skel: estimated DAG skeleton of samples x metabolites data matrix
  # dat: renamed samples x metabolites data matrix
  # dat_compl: renamed samples x metabolites data matrix (typically same as dat)
  # DE: indicator if direct effects were already identified, here set to NA
  # glmulti_method: choose method for glmulti-function
  # exposure: name of considered exposure variable
  # met_map: mapping information between old and new metabolite names
  # adjust_method: select adjustment method(s), e.g. "fdr" or c("fdr","bonferroni"), see p.adjust.methods for available options
  # round_number: actual round number in ambiguous metabolites loop, initiate with round_number=1
  # adjM_norename: not renamed adjacency matrix of metabolite DAG
  # met_group: name of investigated metabolite group

  netin_amb <- list()
  netin_direct <- list()
  netin_sum <- list()
  con_comp <- list()
  MinMod_list <- list()

  MinMod_exp <- c(paste(colnames(exp_dat), collapse = ", ")) # record covariates for minimal model (without any adjacency set)

  repeat{
    print(paste(round_number, ". iteration", sep = ""))
    if (round_number > 1) print(paste0("Connected component ", CC_ind, sep = ""))

    # estimate direct effects of predefined exposure on each network-variable, causal models that agree with the input-network: models are adjusted for all possible combinations of direct neighbors (==variables in adjacency set) -> Output is multiset of possible effects:
    net_coupler_in_AC <- net.coupler.in(graph_skel = graph_skel, dat = dat, dat_compl = dat_compl, exp_dat = exp_dat, DE = DE, glmulti_method = glmulti_method)

    # return results (e.g. p-values, etc.) for specific exposure:
    sum_netin <- getExp.coef(object = net_coupler_in_AC, outcome = colnames(dat), exposure = exposure)

    # get original metabolite names back:
    sum_netin <- merge(sum_netin, as.data.frame(met_map), by = "Outcome")

    # save current minimal model:
    MinMod_list[[length(MinMod_list) + 1]] <- MinMod_exp
    names(MinMod_list)[length(MinMod_list)] <- paste0(round_number, ". iteration")

    sum_stat_netin <- mult.stat(sum_netin = sum_netin, MinMod = MinMod_exp, adjust_method = adjust_method, round_number = round_number) # calculate summary statistics and determine direct and ambiguous effects

    netin_amb[[length(netin_amb) + 1]] <- sum_stat_netin$amb # summary statistics for metabolites on which predefined exposure has ambiguous effect
    netin_direct[[length(netin_direct) + 1]] <- sum_stat_netin$direct # summary statistics for metabolites on which predefined exposure has direct effect
    netin_sum[[length(netin_sum) + 1]] <- sum_stat_netin$sum_netin_sum_adj_FV # summary statistics for meatbolites including round-number information
    names(netin_amb)[length(netin_amb)] <- paste0(round_number, ". iteration")
    names(netin_direct)[length(netin_direct)] <- paste0(round_number, ". iteration")
    names(netin_sum)[length(netin_sum)] <- paste0(round_number, ". iteration")

    numb_DE <- dim(netin_direct[[length(netin_direct)]])[1]
    numb_AMB <- dim(netin_amb[[length(netin_amb)]])[1]

    # check if direct effects!=0 and/or ambiguous effects!=0, otherwise quit:
    if (numb_DE == 0 | numb_AMB == 0) {
      cat(paste("\n ***no direct and/or ambiguous effects identified -> stop algorithm after ", round_number, ". iteration*** \n", sep = ""))
      break
    }

    # if there are direct and ambiguous effects of predefined exposure on metabolites -> extract connected components:
    exposure_list <- list(netin_sum[[length(netin_sum)]])
    names(exposure_list) <- exposure

    # extract connected components for direct and ambiguous effects for predefined exposure:
    con_comp[[length(con_comp) + 1]] <- get.con.comp(exposure_names = c(exposure), exposure_list = exposure_list, adjM_norename = adjM_norename, met_group = met_group)
    names(con_comp)[length(con_comp)] <- paste0(round_number, ". iteration")

    CC_NC_res <- con_comp[[length(con_comp)]]$NC_res[[1]] # original summary statistics for all metabolites including membership in specific connected components for redmeat

    EXPOSURECC <- CC_NC_res %>% dplyr::select(Metabolite, contains("CC")) # Metabolites and membership in specific connected component for redmeat

    # check if current connected component index number CC_ind>dim(EXPOSURECC)[2]-1, if TRUE quit:
    if (CC_ind > dim(EXPOSURECC)[2] - 1) {
      cat(paste("\n ***maximum number of connected components reached -> stop algorithm after evaluating ", CC_ind - 1, ". connected component*** \n", sep = ""))
      break
    }

    names_EXPOSURECC <- EXPOSURECC[which(EXPOSURECC[, CC_ind + 1] == TRUE), 1] # Metabolites which are in specific connected component for redmeat !!!attention, only for one connected component!!!!, loop if more than one?

    # intersection of metabolites with direct effect from predefined exposure and connected components: are there direct effect metabolites which are also in specific connected components, i.e. are connected to at least one other metabolite?
    DE_norename <- intersect(as.character(unlist(names_EXPOSURECC)), as.character(netin_direct[[length(netin_direct)]]$Metabolite))

    DE <- c(DE, met_map[which(met_map[, 1] %in% DE_norename), 2]) # get back short metabolite names

    # intersection of metabolites with ambiguous effect from predefined exposure and connected components: are there ambiguous effect metabolites which are also in specific connected components, i.e. are connected to at least one other metabolite?
    AMB_norename <- intersect(as.character(unlist(names_EXPOSURECC)), as.character(netin_amb[[length(netin_amb)]]$Metabolite))

    AMB <- met_map[which(met_map[, 1] %in% AMB_norename), 2] # get back short metabolite names

    numb_DE_intersect <- length(DE)
    numb_AMB_intersect <- length(AMB)

    # check if intersect(direct effects,connected components)!=0 and/or intersect(ambiguous effects,connected components)!=0, otherwise quit:
    if (numb_DE_intersect == 0 | numb_AMB_intersect == 0) {
      cat(paste("\n ***no direct and/or ambiguous effects in specific connected component -> stop algorithm after ", round_number, ". iteration*** \n", sep = ""))
      break
    }

    # define new minimal model, now including connected component direct effect metabolite as fixed variable:
    MinMod_exp <- paste0(paste0(DE, collapse = ", "), ", ", MinMod_exp)

    # define new metabolite data matrix only including ambiguous effect metabolites from predefined exposure which are also in specific connected component:
    dat <- as.matrix(dat[, c(which(colnames(dat) %in% AMB)), drop = FALSE])

    # add connected component direct effect metabolite data to exp_dat:
    exp_dat <- cbind(dat_compl[, c(which(colnames(dat_compl) %in% DE))], exp_dat)
    colnames(exp_dat)[1:length(DE)] <- DE

    round_number <- round_number + 1
  }

  return(list(netin_direct = netin_direct, netin_amb = netin_amb, netin_sum = netin_sum, con_comp = con_comp, MinMod_list = MinMod_list))
}

# Ambiguous metabolites loop for several connected components#

amb.met.loop.CC <- function(exp_dat, graph_skel, dat, dat_compl, DE, glmulti_method, exposure, met_map, adjust_method, round_number, adjM_norename, met_group) {

  # exp_dat: exposure/phenotype data
  # graph_skel: estimated DAG skeleton of samples x metabolites data matrix
  # dat: renamed samples x metabolites data matrix
  # dat_compl: renamed samples x metabolites data matrix (typically same as dat)
  # DE: indicator if direct effects were already identified, here set to NA
  # glmulti_method: choose method for glmulti-function
  # exposure: name of considered exposure variable
  # met_map: mapping information between old and new metabolite names
  # adjust_method: select adjustment method(s), e.g. "fdr" or c("fdr","bonferroni"), see p.adjust.methods for available options
  # round_number: actual round number in ambiguous metabolites loop, initiate with round_number=1
  # adjM_norename: not renamed adjacency matrix of metabolite DAG
  # met_group: name of investigated metabolite group

  amb_met_loop <- list()

  CC_numb_max <- dim(dat)[2]

  for (i in 1:CC_numb_max) {
    amb_met_loop[[i]] <- amb.met.loop(exp_dat = exp_dat, graph_skel = graph_skel, dat = dat, dat_compl = dat_compl, DE = DE, glmulti_method = glmulti_method, exposure = exposure, met_map = met_map, adjust_method = adjust_method, round_number = round_number, adjM_norename = adjM_norename, met_group = met_group, CC_ind = i)
    names(amb_met_loop)[i] <- paste0("conn_comp_", i)

    if (length(amb_met_loop[[i]]$con_comp) == 0) {
      break
    }
    if (length(amb_met_loop[[i]]$con_comp[[1]]$all_CC) < i) {
      break
    }
  }

  return(amb_met_loop)
}
