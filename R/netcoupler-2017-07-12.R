# HZ Clemens Wittenbecher netcoupler.in for phosphatidylcholines#

# longer analysis version#

# Obtain partial correlation matrix, DAG skeleton, DAG and adjacency matrix for DAG skeleton#

est.pcor.skel.DAG.adj <- function(dat) {

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
    labels = V_met, method = "stable", alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE
  )

  DAG_est <- pc(
    suffStat = list(C = cor(dat), n = n_samples), # estimate equivalence class of DAG using PC-algorithm
    indepTest = gaussCItest, labels = V_met, skel.method = "stable", # order-independent skeleton
    alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, maj.rule = FALSE, solve.confl = FALSE
  )

  adj_matrix <- get.adjacency(graphNEL2igraph(skel_est@graph)) # return adjacency matrix of DAG skeleton

  return(list(pCor_mat = pCor_mat, skel_est = skel_est, DAG_est = DAG_est, adj_matrix = adj_matrix))
}

# Net Coupler IN function: exposure -> metabolome#

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

# Import data + data analysis#

# metabolite data (e.g. phosphatidylcholines) data for subcohort:

met_data_SC <- readRDS("C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/R_data/STR_GPL_SC.rds")

dim(met_data_SC) # metabolite data with rows=patients and columns=metabolites


is.numeric(as.matrix(met_data_SC))


# filter out metabolites containing "_aa_":
met_data_SC <- dplyr::select(met_data_SC, contains("_aa_"))

dim(met_data_SC)


Exp_data_SC <- readRDS("C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/R_data/EXP_SC.rds") # phenotype data for 2092 patients, including diet information, lifestyle, etc.

is.numeric(as.matrix(Exp_data_SC))


dim(Exp_data_SC)


met_data_SC_rename <- rename.met(dat = met_data_SC)$data_renamed # rename metabolites with short names
met_mapping <- rename.met(dat = met_data_SC)$names_mapping # mapping information between old and new metabolite names

met_estimates_norename <- est.pcor.skel.DAG.adj(dat = met_data_SC)
met_skel_norename <- met_estimates_norename$skel_est # estimate DAG skeleton for non-renamed matrix
met_DAG_norename <- met_estimates_norename$DAG_est # estimate DAG for non-renamed matrix
met_adj_mat_norename <- met_estimates_norename$adj_matrix # estimate adjacency matrix for non-renamed matrix

met_estimates_rename <- est.pcor.skel.DAG.adj(dat = met_data_SC_rename)
met_skel_rename <- met_estimates_rename$skel_est # estimate DAG skeleton for renamed matrix
met_DAG_rename <- met_estimates_rename$DAG_est # estimate DAG for renamed matrix
met_adj_mat_rename <- met_estimates_rename$adj_matrix # estimate adjacency matrix for renamed matrix

# visualisation of skeletons and DAGs?

# longer analysis version#

MinMod_exp <- c(paste(colnames(Exp_data_SC), collapse = ", ")) # record covariates for minimal model without any adjacency set

# estimate direct effects of predefined exposure on each network-variable, causal models that agree with the input-network: models are adjusted for all possible combinations of direct neighbors (==variables in adjacency set) -> Output is multiset of possible effects:
net_coupler_in_PC <- net.coupler.in(graph_skel = met_skel_rename, dat = met_data_SC_rename, dat_compl = met_data_SC_rename, exp_dat = Exp_data_SC, DE = NA)

# 1.) return results (e.g. p-values, etc.) for whole-grain bread:
sum_netin_WGB <- getExp.coef(object = net_coupler_in_PC, outcome = colnames(met_data_SC_rename), exposure = "WGBperMJ")

# get original metabolite names back:
sum_netin_WGB <- merge(sum_netin_WGB, as.data.frame(met_mapping), by = "Outcome")

sum_stat_netin_WGB <- mult.stat(sum_netin = sum_netin_WGB, MinMod = MinMod_exp, adjust_method = "fdr", round_number = 1) # calculate summary statistics and determine direct and ambiguous effects

netin_WGB_amb <- sum_stat_netin_WGB$amb # summary statistics for metabolites on which WGB has ambiguous effect
netin_WGB_direct <- sum_stat_netin_WGB$direct # summary statistics for metabolites on which WGB has direct effect
netin_WGB_sum <- sum_stat_netin_WGB$sum_netin_sum_adj_FV # summary statistics for all metabolites including round-number information

netin_WGB_direct
# A tibble: 0 x 13
# ... with 13 variables: Outcome <chr>, avgEst <dbl>, lowEst <dbl>, highEst <dbl>, upperP <dbl>, lowerP <dbl>, Nbmds <dbl>, Metabolite <fctr>, bestGuess <dbl>, marg_P <dbl>, fdr <dbl>, DE <fctr>, Assoc <fctr>

# no direct effect of WGB on metabolites -> ambiguous effects still classified as ambiguous

# 2.) return results (e.g. p-values, etc.) for redmeat:
sum_netin_redmeat <- getExp.coef(object = net_coupler_in_PC, outcome = colnames(met_data_SC_rename), exposure = "TMperMJ")

# get original metabolite names back:
sum_netin_redmeat <- merge(sum_netin_redmeat, as.data.frame(met_mapping), by = "Outcome")

sum_stat_netin_redmeat <- mult.stat(sum_netin = sum_netin_redmeat, MinMod = MinMod_exp, adjust_method = "fdr", round_number = 1) # calculate summary statistics and determine direct and ambiguous effects

netin_redmeat_amb <- sum_stat_netin_redmeat$amb # summary statistics for metabolites on which redmeat has ambiguous effect
netin_redmeat_direct <- sum_stat_netin_redmeat$direct # summary statistics for metabolites on which redmeat has direct effect
netin_redmeat_sum <- sum_stat_netin_redmeat$sum_netin_sum_adj_FV # summary statistics for meatbolites including round-number information

# there are direct effects of redmeat on metabolites -> extract connected components:

# extract connected components for direct effects for redmeat:
con_comp_redmeat <- get.con.comp(exposure_names = c("Redmeat"), exposure_list = list(Redmeat = netin_redmeat_sum), adjM_norename = met_adj_mat_norename, met_group = "PC")

print.data.frame(con_comp_redmeat$NC_res$Redmeat)
# Outcome       avgEst      lowEst      highEst      upperP       lowerP Nbmds   Metabolite    bestGuess       marg_P          fdr DE Assoc round Assoc_FV DE_FV RedmeatCC_1 RedmeatCC_2 RedmeatCC_3
# 1      NM1 -0.089139617 -0.14810258 -0.040740688 0.330629325 3.022918e-03     4 rPC_aa_C28_1 -0.148102579 3.022918e-03 1.165574e-02  0     2     1        2     0        TRUE       FALSE       FALSE
# 2     NM10  0.053565307 -0.10843123  0.172759089 0.981629101 5.618086e-11    32 rPC_aa_C34_4  0.003632206 9.444157e-01 9.664208e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 3     NM11  0.131546327 -0.02042764  0.283578115 0.933678908 4.543472e-13    32 rPC_aa_C36_0  0.277431952 3.494881e-08 5.941297e-07  0     1     1        1     0       FALSE        TRUE       FALSE
# 4     NM12  0.034095881 -0.02965084  0.129162623 0.910157298 8.392635e-04    32 rPC_aa_C36_1  0.005611105 9.101573e-01 9.664208e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 5     NM13  0.051615957  0.04823385  0.054400871 0.307391736 2.707360e-02     4 rPC_aa_C36_2  0.052245461 3.073917e-01 5.481595e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 6     NM14  0.003093814 -0.03599697  0.065554697 0.652479804 2.702216e-02     8 rPC_aa_C36_3 -0.024448494 6.317093e-01 7.954858e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 7     NM15  0.099624690 -0.02349022  0.186826208 0.711746506 5.893559e-09    32 rPC_aa_C36_4  0.149205752 3.085342e-03 1.165574e-02  0     1     1        1     0       FALSE       FALSE        TRUE
# 8     NM16 -0.016142578 -0.06140433  0.022586982 0.957515945 3.262945e-02    16 rPC_aa_C36_5 -0.007397815 8.821536e-01 9.664208e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 9     NM17 -0.045001726 -0.09111519 -0.010748364 0.766712338 5.651188e-04    16 rPC_aa_C36_6 -0.036381765 4.748149e-01 7.019002e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 10    NM18  0.204422580  0.08603102  0.325098232 0.001768967 9.751138e-22    64 rPC_aa_C38_0  0.323648267 1.389026e-10 4.722688e-09  1     1     1        1     1       FALSE        TRUE       FALSE
# 11    NM19  0.091540870  0.05243775  0.163780050 0.280703917 9.842993e-04     8 rPC_aa_C38_1  0.161980274 1.542653e-03 7.492886e-03  0     1     1        1     0       FALSE        TRUE       FALSE
# 12     NM2 -0.110187146 -0.19980810 -0.027820347 0.319109283 4.964666e-11    64 rPC_aa_C30_0 -0.196002905 9.080851e-05 7.718724e-04  0     2     1        2     0        TRUE       FALSE       FALSE
# 13    NM20 -0.038500566 -0.13825477  0.033515690 0.768350798 6.038372e-05    16 rPC_aa_C38_3  0.014888891 7.683508e-01 9.329974e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 14    NM21  0.137678149  0.08350960  0.208089762 0.002952730 9.750752e-09     8 rPC_aa_C38_4  0.208089762 3.951107e-05 4.477922e-04  1     1     1        1     1       FALSE       FALSE        TRUE
# 15    NM22  0.061649027  0.05035766  0.071317197 0.322446741 4.247927e-05     4 rPC_aa_C38_5  0.050357662 3.224467e-01 5.481595e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 16    NM23 -0.011529947 -0.14237235  0.097044167 0.946400629 8.540117e-05    32 rPC_aa_C38_6  0.070792212 1.500978e-01 3.402217e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 17    NM24  0.057474067  0.02227120  0.105051816 0.405970116 4.547447e-03     8 rPC_aa_C40_2  0.097821825 5.541365e-02 1.712786e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 18    NM25  0.020194385 -0.09935420  0.113543181 0.945579478 4.280859e-03    32 rPC_aa_C40_3  0.079834650 1.178603e-01 3.082499e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 19    NM26  0.035014694 -0.08952135  0.096030815 0.552996469 8.013584e-05    16 rPC_aa_C40_4  0.074674522 1.353727e-01 3.287622e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 20    NM27 -0.068231451 -0.09092276 -0.025165360 0.621212742 4.871264e-07    16 rPC_aa_C40_5 -0.025165360 6.212127e-01 7.954858e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 21    NM28  0.009291667 -0.03248565  0.048707412 0.608438064 1.669249e-01     4 rPC_aa_C40_6  0.030844005 5.411393e-01 7.666140e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 22    NM29 -0.124373804 -0.25178849 -0.002155635 0.966420792 1.923137e-11    16 rPC_aa_C42_0 -0.002155635 9.664208e-01 9.664208e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 23     NM3 -0.029058813 -0.16234554  0.098857173 0.983628997 1.744621e-06    32 rPC_aa_C32_0 -0.051519327 3.064552e-01 5.481595e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 24    NM30 -0.011567728 -0.17312117  0.058108114 0.888889966 3.651265e-06    16 rPC_aa_C42_1  0.051502590 3.101994e-01 5.481595e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 25    NM31 -0.084539704 -0.17687962  0.010686539 0.879260896 3.506331e-06    64 rPC_aa_C42_2 -0.010358419 8.371960e-01 9.664208e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 26    NM32  0.080479578  0.02979746  0.124468301 0.454165564 9.588317e-04    16 rPC_aa_C42_4  0.090966031 7.858114e-02 2.226466e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 27    NM33  0.004505882 -0.04349376  0.055802183 0.916374107 9.228398e-02     8 rPC_aa_C42_5  0.027967245 5.838146e-01 7.939879e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 28    NM34 -0.050630720 -0.09434677 -0.010952225 0.763332708 1.072353e-02    16 rPC_aa_C42_6 -0.036621113 4.736716e-01 7.019002e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 29     NM4 -0.102752046 -0.16697629 -0.028085103 0.407283364 1.321437e-07     8 rPC_aa_C32_1 -0.166853792 6.696243e-04 3.794538e-03  0     2     1        2     0        TRUE       FALSE       FALSE
# 30     NM5 -0.112017569 -0.18411158 -0.028800902 0.419714800 3.364941e-11     8 rPC_aa_C32_2 -0.180987884 5.948714e-04 3.794538e-03  0     2     1        2     0        TRUE       FALSE       FALSE
# 31     NM6 -0.017299875 -0.10809096  0.047699544 0.718467199 1.516390e-02     8 rPC_aa_C32_3 -0.067040478 1.940454e-01 4.123465e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 32     NM7 -0.006797895 -0.14164445  0.102204170 0.988300703 1.913732e-05    16 rPC_aa_C34_1 -0.041618577 3.983917e-01 6.450151e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 33     NM8 -0.026832642 -0.10083981  0.032224169 0.958833025 5.614083e-04    16 rPC_aa_C34_2 -0.002660506 9.588330e-01 9.664208e-01  0     0     1        0     0       FALSE       FALSE       FALSE
# 34     NM9 -0.040967592 -0.11081642  0.031174952 0.540416046 2.235802e-03     8 rPC_aa_C34_3 -0.110816417 3.037999e-02 1.032920e-01  0     0     1        0     0       FALSE       FALSE       FALSE

CC.redmeat.PC <- con_comp_redmeat$NC_res$Redmeat # original summary statistics for all metabolites including membership in specific connected components for redmeat
REDMEATCC <- CC.redmeat.PC %>% dplyr::select(Metabolite, contains("CC")) # Metabolites and membership in specific connected component for redmeat

names_REDMEATCC_1 <- REDMEATCC[which(REDMEATCC[, 2] == TRUE), 1] # metabolites which are member of first connected component
names_REDMEATCC_2 <- REDMEATCC[which(REDMEATCC[, 3] == TRUE), 1] # metabolites which are member of second connected component
names_REDMEATCC_3 <- REDMEATCC[which(REDMEATCC[, 4] == TRUE), 1] # metabolites which are member of third connected component

# intersection of metabolites with direct effect from redmeat and connected components: are there direct effect metabolites which are also in specific connected components, i.e. are connected to at least one other metabolite?

DE1_redmeat_CC1 <- intersect(as.character(unlist(names_REDMEATCC_1)), as.character(netin_redmeat_direct$Metabolite))

DE1_redmeat_CC1
# character(0)

# no direct effect metabolite in CC1 -> stop algorithm for connected component 1

DE1_redmeat_CC2 <- intersect(as.character(unlist(names_REDMEATCC_2)), as.character(netin_redmeat_direct$Metabolite))

DE1_redmeat_CC2_rename <- met_mapping[which(met_mapping[, 1] %in% DE1_redmeat_CC2), 2] # get back short metabolite names

# intersection of metabolites with ambiguous effect from redmeat and connected components: are there ambiguous effect metabolites which are also in specific connected components, i.e. are connected to at least one other metabolite?
AMB1_redmeat_CC2 <- intersect(as.character(unlist(names_REDMEATCC_2)), as.character(netin_redmeat_amb$Metabolite))

AMB1_redmeat_CC2_rename <- met_mapping[which(met_mapping[, 1] %in% AMB1_redmeat_CC2), 2]

# define new minimal model, now including connected component direct effect metabolite as fixed variable:
MinMod_exp_DE1_redmeat_CC2 <- paste0(paste0(DE1_redmeat_CC2_rename, collapse = ", "), ", ", MinMod_exp)

# define new metabolite data matrix only including ambiguous effect metabolites from redmeat which are also in specific connected component:
met_data_SC_rename_AMB1_redmeat_CC2 <- met_data_SC_rename[, c(which(colnames(met_data_SC_rename) %in% AMB1_redmeat_CC2_rename))]

# add connected component direct effect metabolite data to Exp_data_SC:
Exp_data_SC_DE1_redmeat_CC2 <- cbind(met_data_SC_rename[, c(which(colnames(met_data_SC_rename) == DE1_redmeat_CC2_rename))], Exp_data_SC)
colnames(Exp_data_SC_DE1_redmeat_CC2) <- c(DE1_redmeat_CC2_rename, colnames(Exp_data_SC))

# repeat net.coupler.in with new set of fixed variables ("always"-set consists now of all exposures, all covariates, and all connected components direct effects):
net_coupler_in_PC_redmeat_CC2 <- net.coupler.in(graph_skel = met_skel_rename, dat = met_data_SC_rename_AMB1_redmeat_CC2, dat_compl = met_data_SC_rename, exp_dat = Exp_data_SC_DE1_redmeat_CC2, DE = DE1_redmeat_CC2_rename)

# return results (e.g. p-values, etc.) for redmeat for second net.coupler.in round for CC2:
sum_netin_redmeat_CC2 <- getExp.coef(object = net_coupler_in_PC_redmeat_CC2, outcome = colnames(met_data_SC_rename_AMB1_redmeat_CC2), exposure = "TMperMJ")

# get original metabolite names back:
sum_netin_redmeat_CC2 <- merge(sum_netin_redmeat_CC2, as.data.frame(met_mapping), by = "Outcome")

sum_stat_netin_redmeat_CC2 <- mult.stat(sum_netin = sum_netin_redmeat_CC2, MinMod = MinMod_exp_DE1_redmeat_CC2, adjust_method = "fdr", round_number = 2) # calculate summary statistics and determine direct and ambiguous effects

netin_redmeat_CC2_amb <- sum_stat_netin_redmeat_CC2$amb # summary statistics for metabolites on which redmeat has ambiguous effect
netin_redmeat_CC2_direct <- sum_stat_netin_redmeat_CC2$direct # summary statistics for metabolites on which redmeat has direct effect
netin_redmeat_CC2_sum <- sum_stat_netin_redmeat_CC2$sum_netin_sum_adj_FV # summary statistics for meatbolites including round-number information

netin_redmeat_CC2_amb
# A tibble: 0 x 13
# ... with 13 variables: Outcome <chr>, avgEst <dbl>, lowEst <dbl>, highEst <dbl>, upperP <dbl>, lowerP <dbl>, Nbmds <dbl>, Metabolite <fctr>, bestGuess <dbl>, marg_P <dbl>, fdr <dbl>, DE <fctr>, Assoc <fctr>
netin_redmeat_CC2_direct
# A tibble: 0 x 13
# ... with 13 variables: Outcome <chr>, avgEst <dbl>, lowEst <dbl>, highEst <dbl>, upperP <dbl>, lowerP <dbl>, Nbmds <dbl>, Metabolite <fctr>, bestGuess <dbl>, marg_P <dbl>, fdr <dbl>, DE <fctr>, Assoc <fctr>#

# no additional direct or ambiguous effect metabolites identified, stop algorithm for this connected component

DE1_redmeat_CC3 <- intersect(as.character(unlist(names_REDMEATCC_3)), as.character(netin_redmeat_direct$Metabolite))

DE1_redmeat_CC3_rename <- met_mapping[which(met_mapping[, 1] %in% DE1_redmeat_CC3), 2] # get back short metabolite names

# intersection of metabolites with ambiguous effect from redmeat and connected components: are there ambiguous effect metabolites which are also in specific connected components, i.e. are connected to at least one other metabolite?
AMB1_redmeat_CC3 <- intersect(as.character(unlist(names_REDMEATCC_3)), as.character(netin_redmeat_amb$Metabolite))

AMB1_redmeat_CC3_rename <- met_mapping[which(met_mapping[, 1] %in% AMB1_redmeat_CC3), 2]

# define new minimal model, now including connected component direct effect metabolite as fixed variable:
MinMod_exp_DE1_redmeat_CC3 <- paste0(paste0(DE1_redmeat_CC3_rename, collapse = ", "), ", ", MinMod_exp)

# define new metabolite data matrix only including ambiguous effect metabolites from redmeat which are also in specific connected component:
met_data_SC_rename_AMB1_redmeat_CC3 <- as.matrix(met_data_SC_rename[, c(which(colnames(met_data_SC_rename) %in% AMB1_redmeat_CC3_rename)), drop = FALSE])

# add connected component direct effect metabolite data to Exp_data_SC:
Exp_data_SC_DE1_redmeat_CC3 <- cbind(met_data_SC_rename[, c(which(colnames(met_data_SC_rename) == DE1_redmeat_CC3_rename))], Exp_data_SC)
colnames(Exp_data_SC_DE1_redmeat_CC3) <- c(DE1_redmeat_CC3_rename, colnames(Exp_data_SC))

# repeat net.coupler.in with new set of fixed variables ("always"-set consists now of all exposures, all covariates, and all connected components direct effects):
net_coupler_in_PC_redmeat_CC3 <- net.coupler.in(graph_skel = met_skel_rename, dat = met_data_SC_rename_AMB1_redmeat_CC3, dat_compl = met_data_SC_rename, exp_dat = Exp_data_SC_DE1_redmeat_CC3, DE = DE1_redmeat_CC3_rename)

# return results (e.g. p-values, etc.) for redmeat for second net.coupler.in round for CC2:
sum_netin_redmeat_CC3 <- getExp.coef(object = net_coupler_in_PC_redmeat_CC3, outcome = colnames(met_data_SC_rename_AMB1_redmeat_CC3), exposure = "TMperMJ")

# get original metabolite names back:
sum_netin_redmeat_CC3 <- merge(sum_netin_redmeat_CC3, as.data.frame(met_mapping), by = "Outcome")

sum_stat_netin_redmeat_CC3 <- mult.stat(sum_netin = sum_netin_redmeat_CC3, MinMod = MinMod_exp_DE1_redmeat_CC3, adjust_method = "fdr", round_number = 2) # calculate summary statistics and determine direct and ambiguous effects

netin_redmeat_CC3_amb <- sum_stat_netin_redmeat_CC3$amb # summary statistics for metabolites on which redmeat has ambiguous effect
netin_redmeat_CC3_direct <- sum_stat_netin_redmeat_CC3$direct # summary statistics for metabolites on which redmeat has direct effect
netin_redmeat_CC3_sum <- sum_stat_netin_redmeat_CC3$sum_netin_sum_adj_FV # summary statistics for meatbolites including round-number information

netin_redmeat_CC3_amb
# A tibble: 0 x 13
# ... with 13 variables: Outcome <chr>, avgEst <dbl>, lowEst <dbl>, highEst <dbl>, upperP <dbl>, lowerP <dbl>, Nbmds <dbl>, Metabolite <fctr>, bestGuess <dbl>, marg_P <dbl>, fdr <dbl>, DE <fctr>, Assoc <fctr>
netin_redmeat_CC3_direct
# A tibble: 0 x 13
# ... with 13 variables: Outcome <chr>, avgEst <dbl>, lowEst <dbl>, highEst <dbl>, upperP <dbl>, lowerP <dbl>, Nbmds <dbl>, Metabolite <fctr>, bestGuess <dbl>, marg_P <dbl>, fdr <dbl>, DE <fctr>, Assoc <fctr>

# no additional direct or ambiguous effect metabolites identified, stop algorithm for this connected component

# final results:

netin_redmeat_sum_CC <- rbind(netin_redmeat_CC2_sum, netin_redmeat_CC3_sum)

netin_redmeat_sum_1_2 <- merge(netin_redmeat_sum, netin_redmeat_sum_CC, by = "Metabolite", suffixes = c("_1", "_2"))

netin_redmeat_sum <- netin_redmeat_sum[, c(8, 1:7, 9:16)]
colnames(netin_redmeat_sum)[2:16] <- paste0(colnames(netin_redmeat_sum)[2:16], sep = "_", "1")

final_netin_redmeat <- bind_rows(netin_redmeat_sum[-which(netin_redmeat_sum$Outcome_1 %in% netin_redmeat_sum_CC$Outcome), ], netin_redmeat_sum_1_2)

# add connected component membership:
final_netin_redmeat <- dplyr::left_join(final_netin_redmeat, REDMEATCC, by = "Metabolite")

write.xlsx(final_netin_redmeat, file = "C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/aaPC_redmeat_final_test.xls")

# 3.) return results (e.g. p-values, etc.) for coffee:
sum_netin_coffee <- getExp.coef(object = net_coupler_in_PC, outcome = colnames(met_data_SC_rename), exposure = "CofCup")

# get original metabolite names back:
sum_netin_coffee <- merge(sum_netin_coffee, as.data.frame(met_mapping), by = "Outcome")

sum_stat_netin_coffee <- mult.stat(sum_netin = sum_netin_coffee, MinMod = MinMod_exp, adjust_method = "fdr", round_number = 1) # calculate summary statistics and determine direct and ambiguous effects

netin_coffee_amb <- sum_stat_netin_coffee$amb # summary statistics for metabolites on which coffee has ambiguous effect
netin_coffee_direct <- sum_stat_netin_coffee$direct # summary statistics for metabolites on which coffee has direct effect
netin_coffee_sum <- sum_stat_netin_coffee$sum_netin_sum_adj_FV # summary statistics for meatbolites including round-number information

# there are direct effects of coffee on metabolites -> extract connected components:

# extract connected components for direct effects for coffee:
con_comp_coffee <- get.con.comp(exposure_names = c("Coffee"), exposure_list = list(Coffee = netin_coffee_sum), adjM_norename = met_adj_mat_norename, met_group = "PC")

print.data.frame(con_comp_coffee$NC_res$Coffee)
# Outcome        avgEst        lowEst      highEst      upperP       lowerP Nbmds   Metabolite    bestGuess       marg_P        fdr DE Assoc round Assoc_FV DE_FV CoffeeCC_1
# 1      NM1  0.0084839568 -0.0149494262  0.039377021 0.836651745 2.962971e-01     4 rPC_aa_C28_1 -0.009286606 0.8366517452 0.88894248  0     0     1        0     0      FALSE
# 2     NM10 -0.0423546804 -0.1341632261 -0.012158701 0.556338877 4.186687e-03    32 rPC_aa_C34_4 -0.134163226 0.0043765478 0.05936127  0     2     1        2     0       TRUE
# 3     NM11  0.0350312196 -0.0082402234  0.085290982 0.986188566 1.162773e-02    32 rPC_aa_C36_0  0.031393834 0.4878415521 0.72115708  0     0     1        0     0      FALSE
# 4     NM12  0.0076096473 -0.0566179210  0.039533571 0.973958037 7.980634e-02    32 rPC_aa_C36_1 -0.031186351 0.4872869222 0.72115708  0     0     1        0     0      FALSE
# 5     NM13  0.0451177736  0.0345975886  0.059976821 0.414824195 2.213230e-02     4 rPC_aa_C36_2  0.037680399 0.4148241953 0.70569001  0     0     1        0     0      FALSE
# 6     NM14 -0.0100910849 -0.0534263343  0.040802192 0.917903398 5.986016e-02     8 rPC_aa_C36_3 -0.050460383 0.2732301796 0.57312914  0     0     1        0     0      FALSE
# 7     NM15 -0.0042755362 -0.0578633645  0.038720356 0.995465259 9.680467e-02    32 rPC_aa_C36_4 -0.055472207 0.2225961213 0.56889069  0     0     1        0     0      FALSE
# 8     NM16 -0.0069150368 -0.0683978636  0.033051323 0.641335579 8.055660e-02    16 rPC_aa_C36_5 -0.056448129 0.2103357754 0.56889069  0     0     1        0     0      FALSE
# 9     NM17 -0.0372107817 -0.1085932974  0.007760875 0.926078543 6.438721e-03    16 rPC_aa_C36_6 -0.108593297 0.0182107122 0.12383284  0     0     1        0     0      FALSE
# 10    NM18 -0.0046502815 -0.0416027729  0.066596052 0.996147505 4.375051e-02    64 rPC_aa_C38_0  0.011121993 0.8060742208 0.88894248  0     0     1        0     0      FALSE
# 11    NM19 -0.0189724882 -0.0285986028 -0.004753614 0.915558761 5.136120e-01     8 rPC_aa_C38_1 -0.018486582 0.6885835201 0.83613713  0     0     1        0     0      FALSE
# 12     NM2 -0.0070233383 -0.0888416935  0.033432047 0.991703569 2.656991e-02    64 rPC_aa_C30_0 -0.088841693 0.0490937285 0.27819779  0     0     1        0     0      FALSE
# 13    NM20 -0.0007314511 -0.0371947129  0.011060290 0.982355189 4.151118e-01    16 rPC_aa_C38_3 -0.037194713 0.4151117708 0.70569001  0     0     1        0     0      FALSE
# 14    NM21  0.0033401094 -0.0345908573  0.017839180 0.883789268 3.679643e-01     8 rPC_aa_C38_4 -0.034590857 0.4482566727 0.72115708  0     0     1        0     0      FALSE
# 15    NM22 -0.0183878938 -0.0546585272  0.002823182 0.851350105 2.342491e-01     4 rPC_aa_C38_5 -0.054658527 0.2342491088 0.56889069  0     0     1        0     0      FALSE
# 16    NM23 -0.0321111325 -0.1016160987  0.001251200 0.986142581 2.709050e-03    32 rPC_aa_C38_6 -0.081671739 0.0659458167 0.27843283  0     0     1        0     0      FALSE
# 17    NM24  0.0005820456 -0.0086654868  0.016300332 0.983853721 7.081491e-01     8 rPC_aa_C40_2  0.016300332 0.7235499410 0.84829993  0     0     1        0     0      FALSE
# 18    NM25  0.0318922931  0.0118102469  0.082340461 0.734804240 1.301641e-02    32 rPC_aa_C40_3  0.027821405 0.5459851907 0.74253986  0     0     1        0     0      FALSE
# 19    NM26 -0.0047690596 -0.0487211497  0.030564342 0.962349959 2.804648e-01    16 rPC_aa_C40_4 -0.048721150 0.2804648359 0.57312914  0     0     1        0     0      FALSE
# 20    NM27 -0.0005126620 -0.0489290704  0.014919495 0.922186392 2.873096e-01    16 rPC_aa_C40_5 -0.048929070 0.2873095834 0.57312914  0     0     1        0     0      FALSE
# 21    NM28 -0.0323349399 -0.0779796744 -0.003194923 0.854653055 8.713226e-02     4 rPC_aa_C40_6 -0.077979674 0.0871322579 0.27843283  0     0     1        0     0      FALSE
# 22    NM29  0.0584777464  0.0380329962  0.078384864 0.332257595 6.215013e-03    16 rPC_aa_C42_0  0.078384864 0.0900812086 0.27843283  0     0     1        0     0      FALSE
# 23     NM3  0.0016009036 -0.0502705931  0.032317176 0.990743022 1.041128e-01    32 rPC_aa_C32_0 -0.046805431 0.3034213069 0.57312914  0     0     1        0     0      FALSE
# 24    NM30 -0.0116481230 -0.0434046164  0.022450337 0.844822628 7.442513e-02    16 rPC_aa_C42_1  0.022450337 0.6241172173 0.78592538  0     0     1        0     0      FALSE
# 25    NM31  0.0167052411 -0.0006175696  0.052620968 0.999833517 1.767983e-01    64 rPC_aa_C42_2  0.024321764 0.5930700115 0.77555309  0     0     1        0     0      FALSE
# 26    NM32  0.0119297811 -0.0097574262  0.032531339 0.926142716 3.384697e-01    16 rPC_aa_C42_4  0.004326435 0.9261427157 0.93482610  0     0     1        0     0      FALSE
# 27    NM33 -0.0597927977 -0.0977533231 -0.023326212 0.506123143 3.195701e-03     8 rPC_aa_C42_5 -0.078394345 0.0890675227 0.27843283  0     0     1        0     0      FALSE
# 28    NM34  0.0129112795 -0.0433320631  0.052874923 0.750051437 4.920608e-02    16 rPC_aa_C42_6 -0.028966300 0.5301608676 0.74253986  0     0     1        0     0      FALSE
# 29     NM4 -0.1029311164 -0.1583227827 -0.072300673 0.004100293 6.800801e-05     8 rPC_aa_C32_1 -0.158322783 0.0003506137 0.01192087  2     2     1        2     2      FALSE
# 30     NM5 -0.0414341039 -0.1327900048 -0.013125900 0.555200710 5.237759e-03     8 rPC_aa_C32_2 -0.132790005 0.0052377592 0.05936127  0     2     1        2     0       TRUE
# 31     NM6  0.0422135892 -0.0035997681  0.091874345 0.928533025 4.850752e-03     8 rPC_aa_C32_3  0.010705827 0.8182761699 0.88894248  0     0     1        0     0      FALSE
# 32     NM7 -0.0067320424 -0.0789324766  0.057536808 0.267322399 1.251074e-02    16 rPC_aa_C34_1 -0.078932477 0.0761425279 0.27843283  0     0     1        0     0      FALSE
# 33     NM8  0.0194471325 -0.0271564655  0.048865073 0.949482752 7.312927e-02    16 rPC_aa_C34_2  0.003805232 0.9348260977 0.93482610  0     0     1        0     0      FALSE
# 34     NM9 -0.0642663199 -0.1199641283 -0.008188597 0.764397736 1.785879e-04     8 rPC_aa_C34_3 -0.112367032 0.0150428440 0.12383284  0     0     1        0     0      FALSE

CC.coffee.PC <- con_comp_coffee$NC_res$Coffee # original summary statistics for all metabolites including membership in specific connected components for coffee
COFFEECC <- CC.coffee.PC %>% dplyr::select(Metabolite, contains("CC")) # Metabolites and membership in specific connected component for coffee

names_COFFEECC_1 <- COFFEECC[which(COFFEECC[, 2] == TRUE), 1] # metabolites which are member of first connected component

# intersection of metabolites with direct effect from coffee and connected components: are there direct effect metabolites which are also in specific connected components, i.e. are connected to at least one other metabolite?

DE1_coffee_CC1 <- intersect(as.character(unlist(names_COFFEECC_1)), as.character(netin_coffee_direct$Metabolite))

DE1_coffee_CC1
# character(0)

# no direct effect metabolite in CC1 -> stop algorithm for connected component 1

#
# shorter analysis version#

# Obtain partial correlation matrix, DAG skeleton, DAG and adjacency matrix for DAG skeleton#

est.pcor.skel.DAG.adj <- function(dat) {

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
    labels = V_met, method = "stable", alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE
  )

  DAG_est <- pc(
    suffStat = list(C = cor(dat), n = n_samples), # estimate equivalence class of DAG using PC-algorithm
    indepTest = gaussCItest, labels = V_met, skel.method = "stable", # order-independent skeleton
    alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, maj.rule = FALSE, solve.confl = FALSE
  )

  adj_matrix <- get.adjacency(graphNEL2igraph(skel_est@graph)) # return adjacency matrix of DAG skeleton

  return(list(pCor_mat = pCor_mat, skel_est = skel_est, DAG_est = DAG_est, adj_matrix = adj_matrix))
}

# Net Coupler IN function: exposure -> metabolome#

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

amb.met.loop <- function(exp_dat, graph_skel, dat, dat_compl, DE, exposure, met_map, adjust_method, round_number, adjM_norename, met_group, CC_ind) {

  # exp_dat: exposure/phenotype data
  # graph_skel: estimated DAG skeleton of samples x metabolites data matrix
  # dat: renamed samples x metabolites data matrix
  # dat_compl: renamed samples x metabolites data matrix (typically same as dat)
  # DE: indicator if direct effects were already identified, here set to NA
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
    net_coupler_in_AC <- net.coupler.in(graph_skel = graph_skel, dat = dat, dat_compl = dat_compl, exp_dat = exp_dat, DE = DE)

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

amb.met.loop.CC <- function(exp_dat, graph_skel, dat, dat_compl, DE, exposure, met_map, adjust_method, round_number, adjM_norename, met_group) {

  # exp_dat: exposure/phenotype data
  # graph_skel: estimated DAG skeleton of samples x metabolites data matrix
  # dat: renamed samples x metabolites data matrix
  # dat_compl: renamed samples x metabolites data matrix (typically same as dat)
  # DE: indicator if direct effects were already identified, here set to NA
  # exposure: name of considered exposure variable
  # met_map: mapping information between old and new metabolite names
  # adjust_method: select adjustment method(s), e.g. "fdr" or c("fdr","bonferroni"), see p.adjust.methods for available options
  # round_number: actual round number in ambiguous metabolites loop, initiate with round_number=1
  # adjM_norename: not renamed adjacency matrix of metabolite DAG
  # met_group: name of investigated metabolite group

  amb_met_loop <- list()

  CC_numb_max <- dim(dat)[2]

  for (i in 1:CC_numb_max) {
    amb_met_loop[[i]] <- amb.met.loop(exp_dat = exp_dat, graph_skel = graph_skel, dat = dat, dat_compl = dat_compl, DE = DE, exposure = exposure, met_map = met_map, adjust_method = adjust_method, round_number = round_number, adjM_norename = adjM_norename, met_group = met_group, CC_ind = i)
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

#

# Import data + data analysis#

# metabolite data (e.g. phosphatidylcholines) data for subcohort:

met_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/GPL/STR_GPL_SC.rds")
dim(met_data_SC) # metabolite data with rows=patients and columns=metabolites


is.numeric(as.matrix(met_data_SC))


# filter out metabolites containing "_aa_":
met_data_SC <- dplyr::select(met_data_SC, contains("_aa_"))

dim(met_data_SC)


Exp_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/GPL/EXP_SC.rds") # phenotype data for 2092 patients, including diet information, lifestyle, etc.

is.numeric(as.matrix(Exp_data_SC))


dim(Exp_data_SC)


met_data_SC_rename <- rename.met(dat = met_data_SC)$data_renamed # rename metabolites with short names
met_mapping <- rename.met(dat = met_data_SC)$names_mapping # mapping information between old and new metabolite names

met_estimates_norename <- est.pcor.skel.DAG.adj(dat = met_data_SC)
met_skel_norename <- met_estimates_norename$skel_est # estimate DAG skeleton for non-renamed matrix
met_DAG_norename <- met_estimates_norename$DAG_est # estimate DAG for non-renamed matrix
met_adj_mat_norename <- met_estimates_norename$adj_matrix # estimate adjacency matrix for non-renamed matrix

met_estimates_rename <- est.pcor.skel.DAG.adj(dat = met_data_SC_rename)
met_skel_rename <- met_estimates_rename$skel_est # estimate DAG skeleton for renamed matrix
met_DAG_rename <- met_estimates_rename$DAG_est # estimate DAG for renamed matrix
met_adj_mat_rename <- met_estimates_rename$adj_matrix # estimate adjacency matrix for renamed matrix

# visualisation of skeletons and DAGs?

# 1.) return results (e.g. p-values, etc.) for whole-grain bread:

net_coupler_in_WGB <- amb.met.loop.CC(exp_dat = Exp_data_SC, graph_skel = met_skel_rename, dat = met_data_SC_rename, dat_compl = met_data_SC_rename, DE = NULL, exposure = "WGBperMJ", met_map = met_mapping, adjust_method = "fdr", round_number = 1, adjM_norename = met_adj_mat_norename, met_group = "PC")

final_netin_wgb <- net_coupler_in_WGB$conn_comp_1$netin_sum$`1. iteration`

print.data.frame(net_coupler_in_WGB$conn_comp_1$netin_sum$`1. iteration`)
# Outcome       avgEst      lowEst      highEst     upperP       lowerP Nbmds   Metabolite    bestGuess       marg_P          fdr DE Assoc round Assoc_FV DE_FV
# 1      NM1 -0.009823579 -0.05104039  0.013918336 0.99424760 3.290368e-01     4 rPC_aa_C28_1 -0.051040394 3.290368e-01 4.766413e-01  0     0     1        0     0
# 2     NM10 -0.014991381 -0.08251914  0.017052430 0.93995846 7.915478e-02    32 rPC_aa_C34_4 -0.082519139 1.308019e-01 2.340666e-01  0     0     1        0     0
# 3     NM11  0.049914839 -0.03632126  0.091128743 0.98230198 9.697192e-04    32 rPC_aa_C36_0 -0.019674100 7.080062e-01 7.522566e-01  0     0     1        0     0
# 4     NM12 -0.137331187 -0.25673591 -0.062726075 0.05228627 1.126354e-09    32 rPC_aa_C36_1 -0.256735907 9.025377e-07 3.068628e-05  0     2     1        2     0
# 5     NM13  0.008954285 -0.08107810  0.102473100 0.18231356 8.357738e-03     4 rPC_aa_C36_2 -0.081078100 1.307447e-01 2.340666e-01  0     0     1        0     0
# 6     NM14 -0.058333168 -0.13721641 -0.013970775 0.55102145 2.680405e-03     8 rPC_aa_C36_3 -0.137216405 1.032136e-02 5.013231e-02  0     2     1        2     0
# 7     NM15 -0.005634925 -0.10457207  0.066430886 0.98475153 4.770072e-02    32 rPC_aa_C36_4 -0.104572070 4.770072e-02 1.474386e-01  0     0     1        0     0
# 8     NM16 -0.077081239 -0.15428204 -0.035542397 0.22638423 1.491172e-03    16 rPC_aa_C36_5 -0.154282041 3.211742e-03 2.729981e-02  0     2     1        2     0
# 9     NM17  0.010378102 -0.07782015  0.054480728 0.80002653 1.229944e-02    16 rPC_aa_C36_6 -0.077820154 1.447564e-01 2.423261e-01  0     0     1        0     0
# 10    NM18 -0.099046509 -0.14342981 -0.046301419 0.22513714 4.229929e-06    64 rPC_aa_C38_0 -0.110528639 3.568428e-02 1.213266e-01  0     0     1        0     0
# 11    NM19 -0.031891285 -0.05989936 -0.009530607 0.84989410 2.633548e-01     8 rPC_aa_C38_1 -0.059899356 2.633548e-01 4.070028e-01  0     0     1        0     0
# 12     NM2  0.021443186 -0.08863843  0.082374675 0.99859806 1.742806e-02    64 rPC_aa_C30_0 -0.088638429 9.075150e-02 2.340666e-01  0     0     1        0     0
# 13    NM20 -0.009682427 -0.11144477  0.076310661 0.89921278 2.021218e-02    16 rPC_aa_C38_3 -0.111444768 3.551515e-02 1.213266e-01  0     0     1        0     0
# 14    NM21 -0.008765847 -0.07630158  0.029815510 0.89554111 1.496720e-01     8 rPC_aa_C38_4 -0.076301576 1.496720e-01 2.423261e-01  0     0     1        0     0
# 15    NM22 -0.062249062 -0.14288768 -0.014841228 0.63279878 7.433345e-03     4 rPC_aa_C38_5 -0.142887679 7.433345e-03 4.212229e-02  0     2     1        2     0
# 16    NM23 -0.023479721 -0.08414678  0.001946311 0.96697679 5.438213e-02    32 rPC_aa_C38_6 -0.084146784 1.026480e-01 2.340666e-01  0     0     1        0     0
# 17    NM24  0.042856065  0.01572987  0.075194514 0.74613879 4.601429e-02     8 rPC_aa_C40_2  0.017318240 7.461388e-01 7.687491e-01  0     0     1        0     0
# 18    NM25  0.005588523 -0.03829563  0.056591009 0.87222708 1.341270e-01    32 rPC_aa_C40_3 -0.023822823 6.560455e-01 7.195338e-01  0     0     1        0     0
# 19    NM26  0.025113345 -0.04877068  0.043385742 0.99990474 7.196271e-02    16 rPC_aa_C40_4 -0.048770683 3.520051e-01 4.766413e-01  0     0     1        0     0
# 20    NM27  0.020988713 -0.08068126  0.104104551 0.87163838 1.060589e-03    16 rPC_aa_C40_5 -0.080681260 1.307333e-01 2.340666e-01  0     0     1        0     0
# 21    NM28 -0.024965627 -0.08006170  0.007775272 0.85370103 1.302429e-01     4 rPC_aa_C40_6 -0.080061705 1.302429e-01 2.340666e-01  0     0     1        0     0
# 22    NM29  0.035754705 -0.02043624  0.124676207 0.84041474 1.306283e-03    16 rPC_aa_C42_0  0.032923331 5.395689e-01 6.391986e-01  0     0     1        0     0
# 23     NM3 -0.043906418 -0.14724607  0.046093611 0.91068471 8.153746e-04    32 rPC_aa_C32_0 -0.145915292 5.753555e-03 3.912417e-02  0     2     1        2     0
# 24    NM30  0.068457727  0.02081664  0.129382694 0.46067872 3.655382e-04    16 rPC_aa_C42_1  0.048477135 3.620692e-01 4.766413e-01  0     0     1        0     0
# 25    NM31 -0.073369476 -0.11342670 -0.005573029 0.90351829 1.182350e-03    64 rPC_aa_C42_2 -0.082919419 1.166428e-01 2.340666e-01  0     0     1        0     0
# 26    NM32  0.023162193 -0.02245694  0.079899265 0.99078855 8.594509e-02    16 rPC_aa_C42_4  0.002647664 9.610239e-01 9.610239e-01  0     0     1        0     0
# 27    NM33  0.007104666 -0.02536650  0.031370457 0.98308216 4.410631e-01     8 rPC_aa_C42_5 -0.025366504 6.354236e-01 7.195338e-01  0     0     1        0     0
# 28    NM34  0.001872159 -0.03408488  0.036153832 0.92457212 3.075569e-01    16 rPC_aa_C42_6 -0.032404709 5.451988e-01 6.391986e-01  0     0     1        0     0
# 29     NM4 -0.063223744 -0.17955699  0.032429462 0.72001948 4.782640e-04     8 rPC_aa_C32_1 -0.179556990 4.782640e-04 5.420326e-03  0     2     1        2     0
# 30     NM5  0.040925513 -0.03819695  0.072447310 0.48864574 1.020594e-02     8 rPC_aa_C32_2 -0.038196949 4.886457e-01 6.153317e-01  0     0     1        0     0
# 31     NM6 -0.042765574 -0.09709138 -0.003622449 0.91609042 6.871028e-02     8 rPC_aa_C32_3 -0.097091378 7.276552e-02 2.061690e-01  0     0     1        0     0
# 32     NM7 -0.086569558 -0.24415275 -0.035485191 0.17150587 2.421440e-06    16 rPC_aa_C34_1 -0.244152750 2.421440e-06 4.116448e-05  0     2     1        2     0
# 33     NM8  0.044985630 -0.04899088  0.074157872 0.63230833 2.004062e-02    16 rPC_aa_C34_2 -0.048990880 3.644904e-01 4.766413e-01  0     0     1        0     0
# 34     NM9 -0.057770393 -0.12836024 -0.009895959 0.70958472 1.907803e-03     8 rPC_aa_C34_3 -0.128360244 1.672905e-02 7.109847e-02  0     2     1        2     0

# no direct effect of WGB on metabolites -> ambiguous effects still classified as ambiguous

# 2.) return results (e.g. p-values, etc.) for redmeat:

net_coupler_in_redmeat <- amb.met.loop.CC(exp_dat = Exp_data_SC, graph_skel = met_skel_rename, dat = met_data_SC_rename, dat_compl = met_data_SC_rename, DE = NULL, exposure = "TMperMJ", met_map = met_mapping, adjust_method = "fdr", round_number = 1, adjM_norename = met_adj_mat_norename, met_group = "PC")

# evaluation:

# number of connected components and membership for redmeat:
CC_memb <- as.data.frame(net_coupler_in_redmeat$conn_comp_1$con_comp$`1. iteration`$all_CC)
CC_memb <- cbind(rownames(CC_memb), CC_memb)
colnames(CC_memb)[1] <- "Metabolite"
rownames(CC_memb) <- NULL

# first iteration, i.e. first net.coupler.in for first connected component:
net_coupler_in_redmeat_initial_CC1 <- net_coupler_in_redmeat$conn_comp_1$netin_sum$`1. iteration`
#*** no direct and/or ambiguous effects in this specific connected component -> stop algorithm after 1. iteration***

# there are direct and ambiguous metabolites in CC2 -> continue with 2. iteration:
net_coupler_in_redmeat_iteration2_CC2 <- net_coupler_in_redmeat$conn_comp_2$netin_sum$`2. iteration`
#*** no direct and/or ambiguous effects identified -> stop algorithm after 2. iteration***

# there are direct and ambiguous metabolites in CC3 -> continue with 2. iteration:
net_coupler_in_redmeat_iteration2_CC3 <- net_coupler_in_redmeat$conn_comp_3$netin_sum$`2. iteration`
#*** no direct and/or ambiguous effects identified -> stop algorithm after 2. iteration***

net_redmeat_sum_CC_short <- rbind(net_coupler_in_redmeat_iteration2_CC2, net_coupler_in_redmeat_iteration2_CC3)

net_redmeat_sum_1_2_short <- merge(net_coupler_in_redmeat_initial_CC1, net_redmeat_sum_CC_short, by = "Metabolite", suffixes = c("_1", "_2"))

net_coupler_in_redmeat_initial_CC1 <- net_coupler_in_redmeat_initial_CC1[, c(8, 1:7, 9:16)]
colnames(net_coupler_in_redmeat_initial_CC1)[2:16] <- paste0(colnames(net_coupler_in_redmeat_initial_CC1)[2:16], sep = "_", "1")

final_netin_redmeat_short <- bind_rows(net_coupler_in_redmeat_initial_CC1[-which(net_coupler_in_redmeat_initial_CC1$Outcome_1 %in% net_redmeat_sum_CC_short$Outcome), ], net_redmeat_sum_1_2_short)

# add connected component membership:
final_netin_redmeat_short <- dplyr::left_join(final_netin_redmeat_short, CC_memb, by = "Metabolite")

write.xlsx(final_netin_redmeat_short, file = "C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/aaPC_redmeat_final_test_ambloop_14_7_2017.xls")

# 3.) return results (e.g. p-values, etc.) for coffee:

net_coupler_in_coffee <- amb.met.loop.CC(exp_dat = Exp_data_SC, graph_skel = met_skel_rename, dat = met_data_SC_rename, dat_compl = met_data_SC_rename, DE = NULL, exposure = "CofCup", met_map = met_mapping, adjust_method = "fdr", round_number = 1, adjM_norename = met_adj_mat_norename, met_group = "PC")

# direct an ambiguous metabolites not in same connected component, stop algorithm

# number of connected components and membership for coffee:
CC_memb_cof <- as.data.frame(net_coupler_in_coffee$conn_comp_1$con_comp$`1. iteration`$all_CC)
CC_memb_cof <- cbind(rownames(CC_memb_cof), CC_memb_cof)
colnames(CC_memb_cof)[1] <- "Metabolite"
rownames(CC_memb_cof) <- NULL

net_coupler_in_coffee_initial <- net_coupler_in_coffee$conn_comp_1$netin_sum$`1. iteration`

# add connected component membership:
final_netin_coffee <- dplyr::left_join(net_coupler_in_coffee_initial, CC_memb_cof, by = "Metabolite")

final_netin_coffee
# Outcome        avgEst        lowEst      highEst      upperP       lowerP Nbmds   Metabolite    bestGuess       marg_P        fdr DE Assoc round Assoc_FV DE_FV CofCupCC_1
# 1      NM1  0.0084839568 -0.0149494262  0.039377021 0.836651745 2.962971e-01     4 rPC_aa_C28_1 -0.009286606 0.8366517452 0.88894248  0     0     1        0     0      FALSE
# 2     NM10 -0.0423546804 -0.1341632261 -0.012158701 0.556338877 4.186687e-03    32 rPC_aa_C34_4 -0.134163226 0.0043765478 0.05936127  0     2     1        2     0       TRUE
# 3     NM11  0.0350312196 -0.0082402234  0.085290982 0.986188566 1.162773e-02    32 rPC_aa_C36_0  0.031393834 0.4878415521 0.72115708  0     0     1        0     0      FALSE
# 4     NM12  0.0076096473 -0.0566179210  0.039533571 0.973958037 7.980634e-02    32 rPC_aa_C36_1 -0.031186351 0.4872869222 0.72115708  0     0     1        0     0      FALSE
# 5     NM13  0.0451177736  0.0345975886  0.059976821 0.414824195 2.213230e-02     4 rPC_aa_C36_2  0.037680399 0.4148241953 0.70569001  0     0     1        0     0      FALSE
# 6     NM14 -0.0100910849 -0.0534263343  0.040802192 0.917903398 5.986016e-02     8 rPC_aa_C36_3 -0.050460383 0.2732301796 0.57312914  0     0     1        0     0      FALSE
# 7     NM15 -0.0042755362 -0.0578633645  0.038720356 0.995465259 9.680467e-02    32 rPC_aa_C36_4 -0.055472207 0.2225961213 0.56889069  0     0     1        0     0      FALSE
# 8     NM16 -0.0069150368 -0.0683978636  0.033051323 0.641335579 8.055660e-02    16 rPC_aa_C36_5 -0.056448129 0.2103357754 0.56889069  0     0     1        0     0      FALSE
# 9     NM17 -0.0372107817 -0.1085932974  0.007760875 0.926078543 6.438721e-03    16 rPC_aa_C36_6 -0.108593297 0.0182107122 0.12383284  0     0     1        0     0      FALSE
# 10    NM18 -0.0046502815 -0.0416027729  0.066596052 0.996147505 4.375051e-02    64 rPC_aa_C38_0  0.011121993 0.8060742208 0.88894248  0     0     1        0     0      FALSE
# 11    NM19 -0.0189724882 -0.0285986028 -0.004753614 0.915558761 5.136120e-01     8 rPC_aa_C38_1 -0.018486582 0.6885835201 0.83613713  0     0     1        0     0      FALSE
# 12     NM2 -0.0070233383 -0.0888416935  0.033432047 0.991703569 2.656991e-02    64 rPC_aa_C30_0 -0.088841693 0.0490937285 0.27819779  0     0     1        0     0      FALSE
# 13    NM20 -0.0007314511 -0.0371947129  0.011060290 0.982355189 4.151118e-01    16 rPC_aa_C38_3 -0.037194713 0.4151117708 0.70569001  0     0     1        0     0      FALSE
# 14    NM21  0.0033401094 -0.0345908573  0.017839180 0.883789268 3.679643e-01     8 rPC_aa_C38_4 -0.034590857 0.4482566727 0.72115708  0     0     1        0     0      FALSE
# 15    NM22 -0.0183878938 -0.0546585272  0.002823182 0.851350105 2.342491e-01     4 rPC_aa_C38_5 -0.054658527 0.2342491088 0.56889069  0     0     1        0     0      FALSE
# 16    NM23 -0.0321111325 -0.1016160987  0.001251200 0.986142581 2.709050e-03    32 rPC_aa_C38_6 -0.081671739 0.0659458167 0.27843283  0     0     1        0     0      FALSE
# 17    NM24  0.0005820456 -0.0086654868  0.016300332 0.983853721 7.081491e-01     8 rPC_aa_C40_2  0.016300332 0.7235499410 0.84829993  0     0     1        0     0      FALSE
# 18    NM25  0.0318922931  0.0118102469  0.082340461 0.734804240 1.301641e-02    32 rPC_aa_C40_3  0.027821405 0.5459851907 0.74253986  0     0     1        0     0      FALSE
# 19    NM26 -0.0047690596 -0.0487211497  0.030564342 0.962349959 2.804648e-01    16 rPC_aa_C40_4 -0.048721150 0.2804648359 0.57312914  0     0     1        0     0      FALSE
# 20    NM27 -0.0005126620 -0.0489290704  0.014919495 0.922186392 2.873096e-01    16 rPC_aa_C40_5 -0.048929070 0.2873095834 0.57312914  0     0     1        0     0      FALSE
# 21    NM28 -0.0323349399 -0.0779796744 -0.003194923 0.854653055 8.713226e-02     4 rPC_aa_C40_6 -0.077979674 0.0871322579 0.27843283  0     0     1        0     0      FALSE
# 22    NM29  0.0584777464  0.0380329962  0.078384864 0.332257595 6.215013e-03    16 rPC_aa_C42_0  0.078384864 0.0900812086 0.27843283  0     0     1        0     0      FALSE
# 23     NM3  0.0016009036 -0.0502705931  0.032317176 0.990743022 1.041128e-01    32 rPC_aa_C32_0 -0.046805431 0.3034213069 0.57312914  0     0     1        0     0      FALSE
# 24    NM30 -0.0116481230 -0.0434046164  0.022450337 0.844822628 7.442513e-02    16 rPC_aa_C42_1  0.022450337 0.6241172173 0.78592538  0     0     1        0     0      FALSE
# 25    NM31  0.0167052411 -0.0006175696  0.052620968 0.999833517 1.767983e-01    64 rPC_aa_C42_2  0.024321764 0.5930700115 0.77555309  0     0     1        0     0      FALSE
# 26    NM32  0.0119297811 -0.0097574262  0.032531339 0.926142716 3.384697e-01    16 rPC_aa_C42_4  0.004326435 0.9261427157 0.93482610  0     0     1        0     0      FALSE
# 27    NM33 -0.0597927977 -0.0977533231 -0.023326212 0.506123143 3.195701e-03     8 rPC_aa_C42_5 -0.078394345 0.0890675227 0.27843283  0     0     1        0     0      FALSE
# 28    NM34  0.0129112795 -0.0433320631  0.052874923 0.750051437 4.920608e-02    16 rPC_aa_C42_6 -0.028966300 0.5301608676 0.74253986  0     0     1        0     0      FALSE
# 29     NM4 -0.1029311164 -0.1583227827 -0.072300673 0.004100293 6.800801e-05     8 rPC_aa_C32_1 -0.158322783 0.0003506137 0.01192087  2     2     1        2     2      FALSE
# 30     NM5 -0.0414341039 -0.1327900048 -0.013125900 0.555200710 5.237759e-03     8 rPC_aa_C32_2 -0.132790005 0.0052377592 0.05936127  0     2     1        2     0       TRUE
# 31     NM6  0.0422135892 -0.0035997681  0.091874345 0.928533025 4.850752e-03     8 rPC_aa_C32_3  0.010705827 0.8182761699 0.88894248  0     0     1        0     0      FALSE
# 32     NM7 -0.0067320424 -0.0789324766  0.057536808 0.267322399 1.251074e-02    16 rPC_aa_C34_1 -0.078932477 0.0761425279 0.27843283  0     0     1        0     0      FALSE
# 33     NM8  0.0194471325 -0.0271564655  0.048865073 0.949482752 7.312927e-02    16 rPC_aa_C34_2  0.003805232 0.9348260977 0.93482610  0     0     1        0     0      FALSE
# 34     NM9 -0.0642663199 -0.1199641283 -0.008188597 0.764397736 1.785879e-04     8 rPC_aa_C34_3 -0.112367032 0.0150428440 0.12383284  0     0     1        0     0      FALSE
