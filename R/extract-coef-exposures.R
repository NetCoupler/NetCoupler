
# 1.3 GETTERS (OUTPUT SUMMARY FUNCTIONS)
# 1.3.1 For specified exposure-outcome pairs: Get Exposure coefficients from all possible models
# A. Get exposure coefficients for a single outcome
getExp.coef.perexposure <- function(object, exposure) {
  # Create a vector containing integers from 1 to number of exposure-specific models
  nbm <- c(1:length(object$Exposures[[exposure]]$Model_summaries))
  Nbmds <- max(nbm)

  # Create emty data-frames for the output
  mm_coef_temp1 <- structure(list(
    Model = as.character(),
    Nbmds = as.numeric(),
    Exposure = as.character(),
    Covariables = as.character(),
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
  # loop along number of exposure-specific models

  for (i in seq(along = nbm))
  {
    SUM <- summary(object$Exposures[[exposure]]$Model_summaries[[i]]$Model_summary)
    Cov <- (dplyr::filter(data.frame(row.names(SUM$coefficients)), row.names(SUM$coefficients) != exposure))

    mm_coef_temp2 <- data.frame(
      Model = as.character(object$Exposures[[exposure]]$Model_summaries[[i]]$Model),
      Nbmds = Nbmds,
      Exposure = as.character(exposure),
      Covariables = as.character(paste(Cov[, 1], collapse = ", ")),
      NCov = max(unlist(object$Exposures[[exposure]]$Model_summaries[[i]]$Model_summary$assign)) - 1,
      HR = as.numeric(SUM$conf.int[[exposure, 1]]),
      LCL = as.numeric(SUM$conf.int[[exposure, 3]]),
      UCL = as.numeric(SUM$conf.int[[exposure, 4]]),
      Beta = as.numeric(SUM$coefficients[[exposure, 1]]),
      rSE = as.numeric(SUM$coefficients[[exposure, 4]]),
      P = as.numeric(SUM$coefficients[[exposure, 6]])
    )

    # bind information to a exposure-specific dataframe
    mm_coef_temp1 <- bind_rows(mm_coef_temp1, mm_coef_temp2)
  }
  mm_coef_temp1
}
