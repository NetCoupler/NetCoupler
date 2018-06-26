# Extract exposure coefficients per outcome, i.e. metabolite#

getExp.coef.permetabolite <- function(object, metabolite) {

  # object: output of net.coupler.out
  # metabolite: specific metabolite to evaluate

  # create vector containing integers from 1 to number of metabolite-specific models:
  nbm <- c(1:length(object$model_details_all[[metabolite]]$Model_summaries))
  Nbmds <- max(nbm)

  # define empty dataframes for output:
  mm_coef_temp1 <- structure(list(
    Model = as.character(),
    Nbmds = as.numeric(),
    Outcome = as.character(),
    Covariables = as.character(),
    # NCov=as.numeric(),
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
    Outcome = as.character(),
    Covariables = as.character(),
    NCov = as.numeric(),
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
    Cov <- rownames(SUM$coefficients)
    # Cov<-dplyr::filter(data.frame(row.names(SUM$coefficients)),row.names(SUM$coefficients)!=exposure)                     #record specific covariates of this model

    mm_coef_temp2 <- data.frame(
      Model = as.character(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model),
      Nbmds = Nbmds,
      Outcome = metabolite,
      Covariables = as.character(paste(Cov, collapse = ", ")),
      NCov = max(unlist(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model_summary$assign)) - 1,
      HR = as.numeric(SUM$conf.int[[metabolite, 1]]),
      LCL = as.numeric(SUM$conf.int[[metabolite, 3]]),
      UCL = as.numeric(SUM$conf.int[[metabolite, 4]]),
      Beta = as.numeric(SUM$coefficients[[metabolite, 1]]),
      rSE = as.numeric(SUM$coefficients[[metabolite, 4]]),
      P = as.numeric(SUM$coefficients[[metabolite, 6]])
    )

    # bind information to an outcome-specific dataframe:
    mm_coef_temp1 <- bind_rows(mm_coef_temp1, mm_coef_temp2)
  }

  return(mm_coef_temp1)
}

