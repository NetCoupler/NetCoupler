
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
