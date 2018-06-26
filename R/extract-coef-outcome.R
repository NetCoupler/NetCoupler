
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

# Extract exposure coefficients per outcome, i.e. metabolite#

# old version? there are a lot of differences between this function and the other three
getExp.coef.peroutcome_old <- function(object, outcome, exposure) {

  # Create a vector containing integers from 1 to number of outcome-specific models
  nbm <- c(1:length(object$Outcomes[[outcome]]$Model_summaries))
  Nbmds <- max(nbm)
  exposure <- exposure
  # Create emty data-frames for the output
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

  # loop along number of outcome-specific models
  for (i in seq(along = nbm))
  { # get exposure-effect estimates (beta coefficient, SE, t-value, p-value) from single model and write into a dataframe
    SUM <- summary(object$Outcomes[[outcome]]$Model_summaries[[i]]$Model_summary)
    Cov <- (dplyr::filter(data.frame(row.names(SUM$coefficients)), row.names(SUM$coefficients) != "exp"))

    mm_coef_temp2 <- data.frame(
      Model = as.character(object$Outcomes[[outcome]]$Model_summaries[[i]]$Model),
      Nbmds = Nbmds,
      Outcome = as.character(object$Outcomes[[outcome]]$Outcome),
      Exposure = exposure,
      Covariables = as.character(paste(Cov[, 1], collapse = ", ")),
      Estimate = as.numeric(SUM$coefficients[paste(exposure), 1]),
      SE = as.numeric(SUM$coefficients[paste(exposure), 2]),
      tval = as.numeric(SUM$coefficients[paste(exposure), 3]),
      P = as.numeric(SUM$coefficients[paste(exposure), 4])
    )

    # bind information to a outcome-specific dataframe
    mm_coef_temp1 <- bind_rows(mm_coef_temp1, mm_coef_temp2)
  }
  mm_coef_temp1
}

