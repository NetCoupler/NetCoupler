# return: estimate of direct effects of each metabolite on survival time, models are adjusted for all possible combinations of direct neighboring metabolites and all covariates

net.coupler.out <- function(graph_skel, dat, dat_compl, exp_dat, DE, survival_obj, always_set, name_surv_obj) {

  # graph_skel: estimated DAG skeleton of samples x metabolites data matrix
  # dat: renamed samples x metabolites data matrix
  # exp_dat: exposure/phenotype data
  # DE: indicator if direct effects were already identified
  # survival_obj: "survival" object
  # always_set: fixed set of covariates always included in model

  cat("*****************************************************************************************************\n")
  cat("This algorithm estimates direct effect of a predefined exposure (network-variable) on time-to-event  \n")
  cat("for all causal models that agree with the input-network: Cox prop. hazards regression models are     \n")
  cat("used to estimate the efect of all network-variables on survival time adjusted for all possible       \n")
  cat("combinations of direct neighbors (adjacency set) -> Output is a multiset of possible causal effects \n")
  cat("*****************************************************************************************************")

  model_details_all <- list(NULL) # prepare empty list to store output
  node_names <- colnames(dat) # create vector "nodes" with node names as strings

  for (i in seq(along = node_names)) { # create empty list with slots for each network-variable, i.e. metabolite

    model_details <- list(NULL)
    model_details_all[[i]] <- model_details
  }
  names(model_details_all) <- node_names

  for (i in node_names) { # net.coupler.out loop over all metabolite nodes

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

    # combine data on exposure (including already identified direct effect metabolites), metabolite outcome, and adjacency set, store as dataframe ("modeldata"), input data for glmulti:
    modeldata <- data.frame(cbind(met_out_data, exp_dat, adjset_data))
    colnames(modeldata)[1] <- met_outcome

    # estimate multimodel coefficients#

    # modify fitting function of coxph so that it always includes one predefined set of variables (here: exposure + covariates), while subsetting adjacent metabolite set:
    coxph.redefined <- function(formula, data, always="", ...) {
      coxph(as.formula(paste(deparse(formula), always)), data = data, ...)
    }

    # fit all possible causal models using glmulti: iterate over all possible combinations of adjacent metabolites, fixed exposure + covariates set

    glmulti_obj <- glmulti(
      y = name_surv_obj, # response variable to be predicted: SURV-object
      xr = c(adjset_names[1:length(adjset_names)]), # predictor variables to be tested in iteration
      data = modeldata, # dataframe containing data
      level = 1, # only main effects (terms of order 1) are used to build the candidate set
      confsetsize = 2^(length(adjset_names)), # number of models to be looked for, dependent on size of adjacency set
      fitfunction = coxph.redefined, # fitting function to be used: here: coxph-function with fixed predictor set "always"
      always = paste0("+", met_outcome, "+", always_set), # defines "always"-set: complete set of exposures and covariates and current metabolite
      includeobjects = TRUE, # fitted models are included as objects
      plotty = FALSE, # no progress of the IC profile is plotted
      report = FALSE # report about progress at run time is given
    )

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
  }

  return(list(model_details_all = model_details_all, outcome = list(survival_obj = survival_obj, class = class(survival_obj))))
}
