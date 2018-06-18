
net.coupler <- function(Graph, Graph_data, Exp_data_SC, exposure) {
  cat("*************************************************************************************************** \n")
  cat("This algorithm estimates direct effects of a predefined exposure on each network-variable for all   \n")
  cat("causal models that agree with the input-network: models are adjusted for all possible combinations  \n")
  cat("of direct neighbors (==variables in the adjecency set) -> Output ist a multiset of possible effects \n")
  cat("***************************************************************************************************")

  mod_coeff_q_new <- c()
  mod_coeff_q <- c()
  Nodes <- c()
  i <- c()

  # PREPARE AN EMPTY LIST TO STORE THE OUTPUT
  Model_details_all <- list(NULL)
  Nodes <- as.character(Graph@graph@nodes) # Create vector "Nodes" with node names as strings
  CC1 <- setdiff(CC$Outcome, DE_exp$Redmeat)
  DE1 <- intersect(CC$Outcome, DE_exp$Redmeat)
  EXP <- colnames(Exp_data_SC)

  for (i in 1:length(CC1)) # Create an empty list with slots each network-variable
  {
    Model_details <- list(NULL)
    Model_details_all[[i]] <- Model_details
  }
  names(Model_details_all) <- CC1

  for (i in CC1)
  {
    # PREPARE DATASETS
    outcome <- i

    if (is.character(outcome) == FALSE) {
      stop("'outcome' is not a character string as required")
    }

    # Select data on Outcome within network, store in "out"
    out <- Graph_data[, outcome]
    if (is.numeric(out) == FALSE) {
      stop("'out' is not a numeric vector as required")
    }

    # create vector with integers indicating adjecent variables
    edgeList <- slot(Graph@graph, "edgeL")
    adjset <- c(edgeList[[outcome]])
    adjset <- c(adjset[[1]])

    # adjdata1<-subset(Graph_data, select=c(adjset)) #Select data on adjecency set, store in adjdata
    # adjdata2<-subset(Graph_data, select=adj_plus$Outcome,sep=","))) #Select data on adjecency set, store in adjdata
    if (is.numeric(out) == FALSE) stop("'adjdata' is not a matrix of numeric variables as required")

    # Create vector with names of adjecency set as strings
    adjset_char <- array(Nodes[adjset])
    adjset_char <- setdiff(adjset_char, DE_exp$Redmeat)
    print(adjset_char)
    if (is.character(adjset_char) == FALSE) {
      stop("'adjset_char' is not a is not a vector of character strings as required")
    }

    # Combine data on Exposure, Outcome and Adjacency set, store as dataframe (modeldata)
    modeldata <- data.frame(cbind(out, Exp_data_SC, Graph_data))
    # print(modeldata[1:10,])
    # modeldata<-rename(modeldata,exp=Exp_data_SC)
    # print(modeldata[1:10,])
    # if(is.(modeldata)==FALSE)
    #  stop("'modeldata' is not a df of numeric variables as required")

    # ESTIMATE MULTIMODEL COEFFICIENTS
    # Modify the fitting function to always include one set of variables while shuffling another set
    glm.redefined <- function(formula, data, always="", ...) {
      glm(as.formula(paste(deparse(formula), always)), data = data, ...)
    }
    # Fit all possible causal models using glmulti
    glmulti_obj <- glmulti(
      y = "out",
      xr = c(adjset_char[1:length(adjset_char)]),
      data = modeldata,
      level = 1,
      confsetsize = 512,
      fitfunc = glm.redefined,
      always = paste0("+", paste0(DE1, collapse = " + "), "+", paste0(EXP, collapse = "+ ")),
      includeobjects = TRUE,
      intercept = FALSE,
      plotty = F,
      report = F
    )

    # OUTPUT SUMMARIES
    # avg_coef<-as.data.frame(coef.glmulti(glmulti_obj, icmethod="Lukacs", alphaIC=0.05));mm_coef<-tbl_df(mm_coef)
    coffee <- glmulti_obj@objects
    Nbmds <- c(1:glmulti_obj@nbmods)
    Model_details <- list(NULL)
    for (j in seq(along = Nbmds))
    {
      Model_details[[j]] <- list(NULL)
    }
    names(Model_details) <- paste("Model", Nbmds, sep = "_")
    for (j in seq(along = Nbmds))
    {
      model_summary <- coffee[[j]]
      Model_details[[j]] <- list(
        Model = paste("Model", j, "of", length(Nbmds)),
        Model_summary = model_summary
      )
      # Model_details[[j]]$Model_summary$data<-NULL
      # Model_details[[j]]$Model_summary$qr$qr<-NULL
      # Model_details[[j]]$Model_summary$model<-NULL
      # Model_details[[j]]$Model_summary$residuals<-NULL
      # Model_details[[j]]$Model_summary$fitted.values<-NULL
      # Model_details[[j]]$Model_summary$effects<-NULL
      # Model_details[[j]]$Model_summary$linear.predictors<-NULL
      # Model_details[[j]]$Model_summary$weights<-NULL
      # Model_details[[j]]$Model_summary$prior.weights<-NULL
      # Model_details[[j]]$Model_summary$y<-NULL
    }
    Model_details_all[[outcome]] <- list(
      Model_summaries = Model_details,
      Number_of_Models = length(Nbmds),
      Adj_set = paste(adjset_char, collapse = ", "),
      Outcome = outcome
    )
  }
  OUT1 <- list(Outcomes = Model_details_all, Exposure = exposure) # list(Model_details_all,mmcf)
  OUT1
}
