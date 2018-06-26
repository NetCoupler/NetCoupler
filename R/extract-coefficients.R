
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

# # B. Get exposure coefficients for a a group of outcomes: e.g. all network-variables
#
# getExp.coef <- function(object, outcome, exposure) {
#   cat("*************************************************************************************************** \n")
#   cat("This function produces a table of effect estimates of the exposure (betas, SE, t-value, p-value) on \n")
#   cat("a selected (set of) outcome(s) [specified as character-vector of variable-names] for all possible   \n")
#   cat("causal models based on conditional independence criteria encoded in the input-network => MULTISET OF\n")
#   cat("POSSIBLE EFFECTS PER EXPOSURE-OUTCOME PAIRS: NUMBER OF POSSIBLE EFFECTS ~ SIZE OF THE ADJECENCY SET!\n")
#   cat("***************************************************************************************************")
#   exposure <- exposure
#   # Define empty dataframes
#   mm_coef <- structure(list(
#     Model = as.character(),
#     Nbmds = as.numeric(),
#     Outcome = as.character(),
#     Exposure = as.character(),
#     Covariables = as.character(),
#     Estimate = as.numeric(),
#     SE = as.numeric(),
#     tval = as.numeric(),
#     P = as.numeric()
#   ),
#   class = "data.frame"
#   )
#   mm_coef_temp <- structure(list(
#     Model = as.character(),
#     Nbmds = as.numeric(),
#     Outcome = as.character(),
#     Exposure = as.character(),
#     Covariables = as.character(),
#     Estimate = as.numeric(),
#     SE = as.numeric(),
#     tval = as.numeric(),
#     P = as.numeric()
#   ),
#   class = "data.frame"
#   )
#   # For each exposure-outcome pair: Get Exposure coefficients from all possible models and write to a table
#   for (j in outcome)
#   {
#     mm_coef_temp <- data.frame(getExp.coef.peroutcome(object = object, outcome = j, exposure = exposure))
#     mm_coef <- bind_rows(mm_coef, mm_coef_temp)
#   }
#   # Merge result tables for all exposure-outcome pairs
#   mm_coef
# }

# old?
getExp.coef_old <- function(object, exposure) {
  SUM <- data.frame(NULL)
  for (j in 1:length(names(object$Outcomes)))
  {
    print(j)
    for (i in 1:PC_IN_CC1$Outcomes[[j]]$Number_of_Models)
    # for (i in 1:PC_IN_CC2$Outcomes[[j]]$Number_of_Models)
    # for (i in 1:PC_IN_CC3$Outcomes[[j]]$Number_of_Models)
    { # get exposure-effect estimates (beta coefficient, SE, t-value, p-value) from single model and write into a dataframe

      SUM1 <- data.frame(NULL)
      SUM1 <- data.frame(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]) %>%
        draw_rownames_Exp() %>%
        dplyr::filter(exposure == Exp) %>%
        dplyr::rename(SE = Std..Error, tval = t.value, P = Pr...t..) %>%
        dplyr::mutate(Outcome = PC_IN_CC1$Outcomes[[j]]$Outcome, Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
        # dplyr::mutate(Outcome = PC_IN_CC2$Outcomes[[j]]$Outcome, Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
        # dplyr::mutate(Outcome = outcome[[j]], Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
        # dplyr::mutate(Outcome = PC_IN_CC3$Outcomes[[j]]$Outcome, Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
      SUM <- dplyr::bind_rows(SUM, SUM1)
    }
  }
  SUM
}
