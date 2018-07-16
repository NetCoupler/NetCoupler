
net.coupler.out <- function(Graph, data) {
  cat("*****************************************************************************************************\n")
  cat("This algorithm estimates direct effect of a predefined exposure (network-variable) on time-to-event  \n")
  cat("for all causal models that agree with the input-network: Cox prop. hazards regression models are     \n")
  cat("used to estimate the efect of all network-variables on survival time adjusted for all possible       \n")
  cat("combinations of direct neighbors (adjecency set) -> Output ist a multiset of possible causal effects \n")
  cat("*****************************************************************************************************")
  mod_coeff_q_new <- c()
  mod_coeff_q <- c()
  Nodes <- c()
  i <- c()

  Model_details_all <- list(NULL)
  Nodes <- as.character(Graph@graph@nodes) # Create vector "Nodes" with node names as strings
  for (i in seq(along = Nodes)) # Create an empty list with slots each network-variable
  {
    Model_details <- list(NULL)
    Model_details_all[[i]] <- Model_details
  }
  names(Model_details_all) <- Nodes

  # Modify the fitting function to Cox phreg + always include one set of variables while shuffling another set
  coxph.redefined <- function(formula, data, always="", ...) {
    coxph(as.formula(paste(deparse(formula), always)), data = data, ...)
  }
  # PREPARE AN EMPTY LIST TO STORE THE OUTPUT

  for (i in seq(along = Nodes))
  {
    Model_details <- list(NULL)
    Model_details_all[[i]] <- Model_details
  }

  names(Model_details_all) <- Nodes

  # ST<-Surv(time = data[["fup_time"]], event=data[["cens"]])
  Nodes <- as.character(Graph@graph@nodes) # Create vector "Nodes" with node names as strings
  for (i in Nodes)
  {

    # PREPARE DATASETS
    exposure <- i
    edgeList <- slot(Graph@graph, "edgeL")
    adjset <- c(edgeList[[i]])
    adjset <- c(adjset[[1]])

    # Create vector with names of adjecency set as strings
    adjset_char <- array(Nodes[adjset])

    glmulti_obj <- glmulti(
      y = "SURV",
      xr = c(adjset_char[1:length(adjset_char)]),
      data = data,
      level = 1,
      confsetsize = 512,
      method = "h",
      fitfunc = coxph.redefined,
      always = paste0(paste0("+", i, collapse = ""), sep = " + ", "WGBperMJ + TMperMJ + CofCup + g2perMJ + g4perMJ + g5perMJ + g6perMJ + g9perMJ +
                                       g10perMJ + g12perMJ + g13perMJ + g15perMJ + g16perMJ + g17perMJ + g18perMJ + g19perMJ  +  g20perMJ + g22perMJ + g24perMJ +
                                       g25perMJ + g26perMJ +  g27perMJ +  g28perMJ + g31perMJ + g34perMJ + g36perMJ + g39perMJ +   g40perMJ + g41perMJ +
                                       g42perMJ + g43perMJ + g45perMJ + g49perMJ +   fasting + GJ_c + bike + sport + alk_1 + alk_2 + alk_3 + alk_4 + alk_5 + alk_6 +
                                       smk1 + smk2 + smk3  +  educ1 +  educ2 +  educ3 + Med_HLipid + Med_Hypert + cluster(ID) + strata(age_years)"),
      plotty = FALSE,
      includeobjects = TRUE
    )

    # OUTPUT SUMMARIES
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
    }

    Model_details_all[[exposure]] <- list(
      Model_summaries = Model_details,
      Number_of_Models = length(Nbmds),
      Adj_set = paste(adjset_char, collapse = ", "),
      Exposure = exposure
    )
  }

  OUT1 <- list(Exposures = Model_details_all, Outcome = Outcome <- list(data = SURV, class = class(SURV)))
  OUT1
}
