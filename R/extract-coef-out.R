
#
# B. Get exposure coefficients for a a group of exposures (e.g. all network-variables)  on time-to-incidence
getExp.coef.out <- function(object, exposure) {
  cat("*************************************************************************************************** \n")
  cat("This function produces a table of effect estimates of all (some) network-variables on an outcome    \n")
  cat("(time-to-event) for all possibles causal models based on conditional independence criteria encoded  \n")
  cat("in the input-network => MULTISET OF POSSIBLE EFFECTS PER EXPOSURE-OUTCOME PAIR; network-variables   \n")
  cat("of interest are selected by indicating the variable-names as character vector                       \n")
  cat("***************************************************************************************************")

  # Define empty dataframes
  mm_coef <- structure(list(
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
  mm_coef_temp <- structure(list(
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
  # For each exposure-outcome pair: Get Exposure coefficients from all possible models and write to a table
  for (j in exposure)
  {
    mm_coef_temp <- data.frame(getExp.coef.perexposure(object = object, exposure = j))
    mm_coef <- bind_rows(mm_coef, mm_coef_temp)
  }
  # Merge result tables for all exposure-outcome pairs
  mm_coef
}
