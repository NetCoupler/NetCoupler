
#' Rename feature names in order to avoid clash with glm.
#'
#' @param dat Dataset of metabolites (columns) by samples (rows)
#'
#' @return Outputs a list of renamed column variables.
#' @export
#'
rename_met <- function(dat) {
  # dat: samples x metabolites data matrix
  Ll <- paste("NM", c(1:dim(dat)[2]), sep = "") # generate shorter metabolite names

  names_mapping <- cbind(colnames(dat), Ll) # mapping of old and new metabolite names
  colnames(names_mapping) <- c("Metabolite", "Outcome")

  data_renamed <- dat
  colnames(data_renamed) <- Ll # is character!

  return(list(data_renamed = data_renamed, names_mapping = names_mapping))
}


# NetCoupler Out code -----------------------------------------------------

#' Estimating the outcome from a NetCoupler DAG.
#'
#' This algorithm estimates direct effect of a predefined exposure
#' (network-variable) on time-to-event for all causal models that agree with the
#' input-network: Cox prop. hazards regression models are used to estimate the
#' efect of all network-variables on survival time adjusted for all possible
#' combinations of direct neighbors (adjacency set) -> Output is a multiset of
#' possible causal effects.
#'
#' @param graph_skel Estimated DAG skeleton of samples x metabolites data matrix
#' @param dat Renamed samples x metabolites data matrix
#' @param adjustment_data Exposure/phenotype data
#' @param DE Indicator if direct effects were already identified
#' @param survival_obj "survival" object
#'
#' @return Outputs a list with model details and outcome estimates.
#' @export
#'
net_coupler_out <- function(graph_skel, dat, adjustment_data, DE, survival_obj) {
    # TODO: DE is variable given?


  # TODO: get metabolic variable names from graph.
  node_names <- colnames(dat)
  # always_set: fixed set of covariates always included in model
  always_set <- paste0(names(adjustment_data), collapse = " + ")

  model_details_all <- purrr::imap(node_names, ~ NULL)
  names(model_details_all) <- node_names

  # for loop with already identified direct effects:
  if (is.null(DE) == FALSE) {
    always_set <- paste0(always_set, " + ", paste0(DE, collapse = " + "), collapse = " + ")
    node_names <- setdiff(node_names, DE)
    adjustment_data <- dplyr::bind_cols(adjustment_data, dplyr::select(dat, noquote(c(DE))))
  }

  node_names_test <- node_names

  for (i in node_names_test) { # net.coupler.out loop over all metabolite nodes

    ########################################## prepare separate datasets for each metabolite:########################################################

    exposure_metabolite <- i

    # check if exposure_metabolite is character string:
    # if (is.character(exposure_metabolite) == FALSE)
    #   stop("'exposure_metabolite' is not a character string as required")

    # select data on exposure_metabolite within samples x metabolites data matrix, store in "exposure_metabolite_data":
    exposure_metabolite_data <- dat[, exposure_metabolite]

    # check if exposure_metabolite_data is numeric:
    # if (is.numeric(exposure_metabolite_data) == FALSE)
    #   stop("'exposure_metabolite_data' is not a numeric vector as required")

    # create vector with integers indicating adjacent variables, i.e. metabolites in skeleton:
    edge_list <- graph_skel@graph@edgeL # extract edge-list from skeleton

    adjset <- edge_list[[exposure_metabolite]][[1]] # extract adjacency set for selected node/metabolite

    match_nodes <- tibble::tibble(
        Names = graph_skel@graph@nodes,
        Index = 1:length(Names)
    )
    # extract indices of adjacency set
    adjset <-
      adjset[adjset != as.numeric(match_nodes %>% dplyr::filter(Names == exposure_metabolite) %>% dplyr::select(Index))]

    if (is.null(DE) == FALSE) {
      DE_index <- match_nodes %>%
        dplyr::filter(Names %in% DE) %>%
        dplyr::select(Index)
      adjset <- adjset[!adjset %in% DE_index$Index]
    }

    # select data on adjacency set, store in "adjset_data":
    if (rlang::is_empty(adjset) == FALSE) {
      adjset_data <- subset(dat, select = c(adjset))


      # check if adjset_data is numeric:
      # if (is.numeric(adjset_data[1, 1]) == FALSE)
      #   stop("'adjset_data' is not a matrix of numeric variables as required")

      # create vector with names of adjacency set as strings:
      adjset_names <- colnames(adjset_data)[!colnames(adjset_data) %in% "ident"]

      # combine data on exposure (including already identified direct effect metabolites), metabolite outcome, and adjacency set, store as dataframe ("modeldata"), input data for glmulti:
      modeldata <- data.frame(cbind(exposure_metabolite_data, adjustment_data, adjset_data))
      colnames(modeldata)[1] <- exposure_metabolite
    }
    else {
      message("No direct neighbors left after removing DE")
      adjset_names <- c("Age")
      modeldata <- data.frame(cbind(exposure_metabolite_data, adjustment_data))
      colnames(modeldata)[1] <- exposure_metabolite
    }

    ########################################################### estimate multimodel coefficients#####################################################
   # message("Is fixed set of covariables okay?")


    coxph.redefined <- function(formula, data, always = "", ...) {
      survival::coxph(as.formula(paste(deparse(formula), always)), data = data, method = c("efron"), ...)
    }

    modeldata$Age <- round(modeldata$Age, digits = 0)
    modeldata$Age
    if (!requireNamespace("glmulti", quietly = TRUE)) {
      stop("This function requires glmulti, please install it.")
    }
    glmulti_obj <- glmulti::glmulti(
      y = survival_obj,
      xr = adjset_names,
      data = modeldata,
      level = 1,
      minsize = 0,
      maxsize = -1,
      minK = 0,
      maxK = -1,
      confsetsize = 2^length(adjset_names),
      method = "h",
      fitfunc = coxph.redefined,
      always = paste0(" + ",
                      exposure_metabolite,
                      " + ",
                      always_set,
                      " + survival::strata(Age)",
                      collapse = " + "),
      plotty = FALSE,
      includeobjects = TRUE,
      report = FALSE,
      intercept = TRUE
    )

    # output summaries:
    glmulti_obj_objects <- glmulti_obj@objects # collect glmulti objects, i.e. all fitted models (?)

    nbmds <- 1:glmulti_obj@nbmods # number of fitted models of glmulti-function

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

    model_details_all[[exposure_metabolite]] <- list(
      # collect all model details for different adjacency subsets for specific metabolite in complete model list for all metabolites
      Model_summaries = model_details,
      Number_of_Models = length(nbmds),
      Adj_set = paste(adjset_names, collapse = ", "),
      Outcome_metabolite = exposure_metabolite
    )
  }

  return(list(
    model_details_all = model_details_all,
    outcome = list(survival_obj = survival_obj, class = class(survival_obj))
  ))
}

#################################################### Extract exposure coefficients per outcome, i.e. metabolite#######################################################

#' Extract the exposure coefficients per outcome from the netcoupler outcome estimation.
#'
#' @param object Output of net.coupler.out
#' @param metabolite Specific metabolite to evaluate
#' @param DE Direct effect...
#'
#' @return Outputs a data frame of estimates.
#' @export
#'
getExp.coef.permetabolite <- function(object, metabolite, DE = NA) {

  # create vector containing integers from 1 to number of metabolite-specific models:
  nbm <- c(1:length(object$model_details_all[[metabolite]]$Model_summaries))
  Nbmds <- max(nbm)

  # define empty dataframes for output:
  mm_coef_temp1 <- structure(list(
    Model = as.character(),
    Nbmds = as.numeric(),
    Exposure = as.character(),
    Neighbors = as.character(),
    Adjustment = as.character(),
    N_Adjustment = as.numeric(),
    HR = as.numeric(),
    LCL = as.numeric(),
    UCL = as.numeric(),
    Beta = as.numeric(),
    rSE = as.numeric(),
    P = as.numeric()
  ),
  class = "data.frame"
  )

  mm_coef_temp2 <- structure(list(
    Model = as.character(),
    Nbmds = as.numeric(),
    Exposure = as.character(),
    Neighbors = as.character(),
    Adjustment = as.character(),
    N_Adjustment = as.numeric(),
    HR = as.numeric(),
    LCL = as.numeric(),
    UCL = as.numeric(),
    Beta = as.numeric(),
    rSE = as.numeric(),
    P = as.numeric()
  ),
  class = "data.frame"
  )

  # loop along number of metabolite-specific models:
  for (i in seq(along = nbm)) {

    # get metabolite-effect estimates (hazard-ratio, lower and upper confidence limit, beta, robust SE, and p-value) from single model and write into dataframe
    SUM <- summary(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model_summary)
    Cov <- c()
    Cov1 <- c()

    All_Neighbors <- c()
    All_Neighbors <- strsplit(object$model_details_all[[metabolite]]$Adj_set, ", ")[[1]] %>%
      str_replace("^[T]*", "") %>%
      str_replace("[D]", "dh") %>%
      str_replace("CER*", "") %>%
      str_replace("a_*", "") %>%
      str_replace("_", ":") %>%
      str_replace("_", "")
    adj_Neighbors <- c()
    adj_Neighbors <-
      rownames(SUM$coefficients)[rownames(SUM$coefficients) %in% strsplit(object$model_details_all[[metabolite]]$Adj_set, ", ")[[1]]] %>%
      str_replace("^[T]*", "") %>%
      str_replace("[D]", "dh") %>%
      str_replace("CER*", "") %>%
      str_replace("a_*", "") %>%
      str_replace("_", ":") %>%
      str_replace("_", "")

    mm_coef_temp2 <- data.frame(
      Model = as.character(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model),
      Nbmds = Nbmds,
      Exposure = metabolite %>%
        str_replace("^[T]*", "") %>%
        str_replace("[D]", "dh") %>%
        str_replace("CER*", "") %>%
        str_replace("a_*", "") %>%
        str_replace("_", ":") %>%
        str_replace("_", ""),
      Neighbors = paste0(All_Neighbors, collapse = " "),
      Adjustment = paste0(adj_Neighbors, collapse = " "),
      N_Adjustment = as.numeric(length(adj_Neighbors)),
      HR = as.numeric(round(SUM$conf.int[[metabolite, 1]], digits = 2)),
      LCL = as.numeric(round(SUM$conf.int[[metabolite, 3]], digits = 2)),
      UCL = as.numeric(round(SUM$conf.int[[metabolite, 4]], digits = 2)),
      Beta = as.numeric(format(SUM$coefficients[[metabolite, 1]], digits = 2, nsmall = 2)),
      rSE = as.numeric(format(SUM$coefficients[[metabolite, 3]], digits = 2, nsmall = 2)),
      P = as.numeric(format(SUM$coefficients[[metabolite, 5]], digits = 3, nsmall = 3))
    )

    # bind information to an Exposure-specific dataframe:
    mm_coef_temp1 <- dplyr::bind_rows(mm_coef_temp1, mm_coef_temp2)
  }

  return(mm_coef_temp1)
}

#################################################### Extract exposure coefficients for metabolite(s) on time-to-incidence###############################################

#' Extract coefficients of effect estimates for outcomes.
#'
#' This function produces a table of effect estimates of all (some)
#' network-variables on an outcome (time-to-event) for all possibles causal
#' models based on conditional independence criteria encoded in the
#' input-network => MULTISET OF POSSIBLE EFFECTS PER EXPOSURE-OUTCOME PAIR;
#' network-variables of interest are selected by indicating the variable-names
#' as character vector.
#'
#' @param object Output of net.coupler
#' @param metabolite Specific metabolite to evaluate
#' @param DE Direct effect ...
#'
#' @return Outputs a data frame of effect estimates.
#' @export
#'
getExp.coef.out <- function(object, metabolite, DE = NA) {

  # define empty dataframes for output:

  mm_coef <- structure(list(
    Model = as.character(),
    Nbmds = as.numeric(),
    Exposure = as.character(),
    Neighbors = as.character(),
    Adjustment = as.character(),
    N_Adjustment = as.numeric(),
    HR = as.numeric(),
    LCL = as.numeric(),
    UCL = as.numeric(),
    Beta = as.numeric(),
    rSE = as.numeric(),
    P = as.numeric()
  ),
  class = "data.frame"
  )

  mm_coef_temp <- structure(list(
    Model = as.character(),
    Nbmds = as.numeric(),
    Exposure = as.character(),
    Neighbors = as.character(),
    Adjustment = as.character(),
    N_Adjustment = as.numeric(),
    HR = as.numeric(),
    LCL = as.numeric(),
    UCL = as.numeric(),
    Beta = as.numeric(),
    rSE = as.numeric(),
    P = as.numeric()
  ),
  class = "data.frame"
  )

  # for each metabolite: get metabolite coefficients from all possible models and write to table:

  for (j in metabolite) {
    mm_coef_temp <- data.frame(getExp.coef.permetabolite(object = object, metabolite = j, DE = DE))
    mm_coef <- dplyr::bind_rows(mm_coef, mm_coef_temp)
  }

  return(mm_coef)
}


##################### Summary statistics and multiple testing adjustment for net.coupler.out with survival object############################################################################

mult.stat.surv <-
  function(sum_netout,
           adjust_method,
           rule_1 = 12,
           rule_1_cut = 0.1,
           rule_2 = 5,
           rule_2_cut = 0.05,
           rule_3 = 15,
           rule_4 = 14,
           ass_rule1 = 16,
           round_number) {

  # sum_netout: output of getExp.coef.out after adding original metabolite names
  # adjust_method: select adjustment method(s), e.g. "fdr" or c("fdr","bonferroni"), see p.adjust.methods for available options
  # rule= "set DE=1 or =2 if ..."
  # rule_1: determines first argument of rule (e.g. adjusted p-value of full model < rule_1_cut)
  # rule_2: determines second argument of rule (e.g. upperP < rule_2_cut)
  # rule_3: determines first variable of third argument of rule (e.g. highEst)
  # rule_4: determines second variable of third argument of rule (e.g. lowEst)
  # possible arguments for rule_...: adjusted p-values of full model=12=rule_1
  #                                 largest p-value across all adjacency models=upperP=5=rule_2
  #                                 highest estimated beta across all adjacency models=highEst=15=rule_3
  #                                 lowest estimated beta across all adjacency models=lowEst=14=rule_4
  # default rule: if(adjusted p< 0.1 & upperP<0.05 & (sign(highEst)==sign(lowEst)) & lowEst>0){DE=1}
  #              if(adjusted p< 0.1 & upperP<0.05 & (sign(highEst)==sign(lowEst)) & lowEst<0){DE=2}else{DE=0}
  # ass_rule1: determines first variable of second argument of association rule (e.g. bestGuess=16)
  # default association rule: if(adjusted p<0.1 & bestGuess>0){Assoc=1}
  #                          if(adjusted p<0.1 & bestGuess<0){Assoc=2}else{Assoc=0}
  # round_number: actual round number in ambiguous metabolites loop

  # return: for each metabolite: maximum number of models, mean hazard ratio, maximum HR, minimum HR, maximum p-value, minimum p-value, as well as number of covariates, HR and p-value for model with full adjustment:
  sum_netout_sum <-
    dplyr::left_join(
      sum_netout %>%
        dplyr::group_by(Exposure) %>%
        dplyr::summarise(
          avgHR = mean(HR),
          minHR = min(HR),
          maxHR = max(HR),
          upperP = max(P),
          lowerP = min(P),
          Nbmds = max(Nbmds)
        ),
      sum_netout %>%
        dplyr::group_by(Exposure) %>%
        dplyr::filter(nchar(Covariables) == max(nchar(Covariables))) %>%
        dplyr::select(
          NCov = NCov,
          Exposure = Exposure,
          ind_HR = HR,
          ind_P = P
        ),
      by = "Exposure"
    )

  # multiple testing correction:
  p_adjust_M <- p.adjust.methods[p.adjust.methods %in% adjust_method] # select multiple testing correction method(s)
  p_adj <- sapply(p_adjust_M, function(meth) {
    stats::p.adjust(sum_netout_sum$ind_P, meth)
  }) # calculate adjusted p-values for ind_P

  sum_netout_sum_adj <- bind_cols(sum_netout_sum, data.frame(p_adj))
  colnames(sum_netout_sum_adj)[12] <- adjust_method

  # add beta-coefficients:
  sum_netout_sum_adj <-
    sum_netout_sum_adj %>% dplyr::mutate(
      avgEst = log(avgHR),
      lowEst = log(minHR),
      highEst = log(maxHR),
      bestGuess = log(ind_HR)
    )

  # is there an effect of specified metabolite on time-to-event?: determine effect-indicator (DE!=0: direct effect of metabolite on time-to-event, DE=0: ambiguous if adjusted ind_p<0.1):
  sum_netout_sum_adj <- mutate(sum_netout_sum_adj,
                               DE = mosaic::derivedFactor(
                                 "1" = (
                                   sum_netout_sum_adj[, rule_1] < rule_1_cut &
                                     sum_netout_sum_adj[, rule_2] < rule_2_cut &
                                     (sign(sum_netout_sum_adj[, rule_3]) == sign(sum_netout_sum_adj[, rule_4])) &
                                     sum_netout_sum_adj[, rule_4] > 0
                                 ),
                                 "2" = (
                                   sum_netout_sum_adj[, rule_1] < rule_1_cut &
                                     sum_netout_sum_adj[, rule_2] < rule_2_cut &
                                     (sign(sum_netout_sum_adj[, rule_3]) == sign(sum_netout_sum_adj[, rule_4])) &
                                     sum_netout_sum_adj[, rule_4] < 0
                                 ),
                                 .method = "first",
                                 .default = 0
                               ))

  # is there any effect of specified metabolite on time-to-event?: determine association-indicator (Assoc!=0: there is an effect (no differenciation between direct or indirect), Assoc=0: there is no effect)
  sum_netout_sum_adj <- mutate(sum_netout_sum_adj,
    Assoc = mosaic::derivedFactor(
      "1" = (sum_netout_sum_adj[, rule_1] < rule_1_cut & sum_netout_sum_adj[, ass_rule1] > 0),
      "2" = (sum_netout_sum_adj[, rule_1] < rule_1_cut & sum_netout_sum_adj[, ass_rule1] < 0),
      .method = "first", .default = 0
    )
  )

  # collect ambiguous effects:
  amb <- sum_netout_sum_adj %>% dplyr::filter(Assoc != 0 & DE == 0)

  # collect direct effects:
  direct <- sum_netout_sum_adj %>% dplyr::filter(Assoc != 0 & DE != 0)

  # collect summary statistics including round information
  sum_netout_sum_adj_FV <-
    sum_netout_sum_adj %>% dplyr::mutate(
      round = as.numeric(round_number),
      Assoc_FV = as.character(Assoc),
      DE_FV = as.character(DE)
    )

  return(
    list(
      sum_netout_sum_adj = sum_netout_sum_adj,
      amb = amb,
      direct = direct,
      sum_netout_sum_adj_FV = sum_netout_sum_adj_FV
    )
  )
}
