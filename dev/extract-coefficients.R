#
# # Extract exposure coefficients for specific exposure and outcomes, i.e. metabolites (?)#
# getExp.coef <- function(object, outcome, exposure) {
#
#   # object: output of net.coupler
#   # outcome: list of metabolites to evaluate, e.g. all
#   # exposure: name of considered exposure variable
#
#   cat("*************************************************************************************************** \n")
#   cat("This function produces a table of effect estimates of the exposure (betas, SE, t-value, p-value) on \n")
#   cat("a selected (set of) outcome(s) [specified as character-vector of variable-names] for all possible   \n")
#   cat("causal models based on conditional independence criteria encoded in the input-network => MULTISET OF\n")
#   cat("POSSIBLE EFFECTS PER EXPOSURE-OUTCOME PAIRS: NUMBER OF POSSIBLE EFFECTS ~ SIZE OF THE ADJECENCY SET!\n")
#   cat("***************************************************************************************************")
#
#   # define empty dataframes for output:
#
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
#
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
#
#   # for each exposure-outcome pair: get exposure coefficients from all possible models and write to table:
#
#   for (j in outcome) {
#     mm_coef_temp <- data.frame(getExp.coef.peroutcome(object = object, outcome = j, exposure = exposure))
#     mm_coef <- bind_rows(mm_coef, mm_coef_temp)
#   }
#
#   return(mm_coef)
# }
#
# # # B. Get exposure coefficients for a a group of outcomes: e.g. all network-variables
# #
# # getExp.coef <- function(object, outcome, exposure) {
# #   cat("*************************************************************************************************** \n")
# #   cat("This function produces a table of effect estimates of the exposure (betas, SE, t-value, p-value) on \n")
# #   cat("a selected (set of) outcome(s) [specified as character-vector of variable-names] for all possible   \n")
# #   cat("causal models based on conditional independence criteria encoded in the input-network => MULTISET OF\n")
# #   cat("POSSIBLE EFFECTS PER EXPOSURE-OUTCOME PAIRS: NUMBER OF POSSIBLE EFFECTS ~ SIZE OF THE ADJECENCY SET!\n")
# #   cat("***************************************************************************************************")
# #   exposure <- exposure
# #   # Define empty dataframes
# #   mm_coef <- structure(list(
# #     Model = as.character(),
# #     Nbmds = as.numeric(),
# #     Outcome = as.character(),
# #     Exposure = as.character(),
# #     Covariables = as.character(),
# #     Estimate = as.numeric(),
# #     SE = as.numeric(),
# #     tval = as.numeric(),
# #     P = as.numeric()
# #   ),
# #   class = "data.frame"
# #   )
# #   mm_coef_temp <- structure(list(
# #     Model = as.character(),
# #     Nbmds = as.numeric(),
# #     Outcome = as.character(),
# #     Exposure = as.character(),
# #     Covariables = as.character(),
# #     Estimate = as.numeric(),
# #     SE = as.numeric(),
# #     tval = as.numeric(),
# #     P = as.numeric()
# #   ),
# #   class = "data.frame"
# #   )
# #   # For each exposure-outcome pair: Get Exposure coefficients from all possible models and write to a table
# #   for (j in outcome)
# #   {
# #     mm_coef_temp <- data.frame(getExp.coef.peroutcome(object = object, outcome = j, exposure = exposure))
# #     mm_coef <- bind_rows(mm_coef, mm_coef_temp)
# #   }
# #   # Merge result tables for all exposure-outcome pairs
# #   mm_coef
# # }
#
# # old?
# getExp.coef_old <- function(object, exposure) {
#   SUM <- data.frame(NULL)
#   for (j in 1:length(names(object$Outcomes)))
#   {
#     print(j)
#     for (i in 1:PC_IN_CC1$Outcomes[[j]]$Number_of_Models)
#     # for (i in 1:PC_IN_CC2$Outcomes[[j]]$Number_of_Models)
#     # for (i in 1:PC_IN_CC3$Outcomes[[j]]$Number_of_Models)
#     { # get exposure-effect estimates (beta coefficient, SE, t-value, p-value) from single model and write into a dataframe
#
#       SUM1 <- data.frame(NULL)
#       SUM1 <- data.frame(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]) %>%
#         draw_rownames_Exp() %>%
#         dplyr::filter(exposure == Exp) %>%
#         dplyr::rename(SE = Std..Error, tval = t.value, P = Pr...t..) %>%
#         dplyr::mutate(Outcome = PC_IN_CC1$Outcomes[[j]]$Outcome, Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
#         # dplyr::mutate(Outcome = PC_IN_CC2$Outcomes[[j]]$Outcome, Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
#         # dplyr::mutate(Outcome = outcome[[j]], Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
#         # dplyr::mutate(Outcome = PC_IN_CC3$Outcomes[[j]]$Outcome, Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
#       SUM <- dplyr::bind_rows(SUM, SUM1)
#     }
#   }
#   SUM
# }
#
# # 1.3 GETTERS (OUTPUT SUMMARY FUNCTIONS)
# # 1.3.1 For specified exposure-outcome pairs: Get Exposure coefficients from all possible models
# # A. Get exposure coefficients for a single outcome
# getExp.coef.perexposure <- function(object, exposure) {
#   # Create a vector containing integers from 1 to number of exposure-specific models
#   nbm <- c(1:length(object$Exposures[[exposure]]$Model_summaries))
#   Nbmds <- max(nbm)
#
#   # Create emty data-frames for the output
#   mm_coef_temp1 <- structure(list(
#     Model = as.character(),
#     Nbmds = as.numeric(),
#     Exposure = as.character(),
#     Covariables = as.character(),
#     HR = as.numeric(),
#     LCL = as.numeric(),
#     UCL = as.numeric(),
#     Beta = as.numeric(),
#     rSE = as.numeric(),
#     P = as.numeric()
#   ),
#   class = "data.frame"
#   )
#
#   mm_coef_temp2 <- structure(list(
#     Model = as.character(),
#     Nbmds = as.numeric(),
#     Exposure = as.character(),
#     Covariables = as.character(),
#     NCov = as.numeric(),
#     HR = as.numeric(),
#     LCL = as.numeric(),
#     UCL = as.numeric(),
#     Beta = as.numeric(),
#     rSE = as.numeric(),
#     P = as.numeric()
#   ),
#   class = "data.frame"
#   )
#   # loop along number of exposure-specific models
#
#   for (i in seq(along = nbm))
#   {
#     SUM <- summary(object$Exposures[[exposure]]$Model_summaries[[i]]$Model_summary)
#     Cov <- (dplyr::filter(data.frame(row.names(SUM$coefficients)), row.names(SUM$coefficients) != exposure))
#
#     mm_coef_temp2 <- data.frame(
#       Model = as.character(object$Exposures[[exposure]]$Model_summaries[[i]]$Model),
#       Nbmds = Nbmds,
#       Exposure = as.character(exposure),
#       Covariables = as.character(paste(Cov[, 1], collapse = ", ")),
#       NCov = max(unlist(object$Exposures[[exposure]]$Model_summaries[[i]]$Model_summary$assign)) - 1,
#       HR = as.numeric(SUM$conf.int[[exposure, 1]]),
#       LCL = as.numeric(SUM$conf.int[[exposure, 3]]),
#       UCL = as.numeric(SUM$conf.int[[exposure, 4]]),
#       Beta = as.numeric(SUM$coefficients[[exposure, 1]]),
#       rSE = as.numeric(SUM$coefficients[[exposure, 4]]),
#       P = as.numeric(SUM$coefficients[[exposure, 6]])
#     )
#
#     # bind information to a exposure-specific dataframe
#     mm_coef_temp1 <- bind_rows(mm_coef_temp1, mm_coef_temp2)
#   }
#   mm_coef_temp1
# }
#
# #
# # B. Get exposure coefficients for a a group of exposures (e.g. all network-variables)  on time-to-incidence
# getExp.coef.out <- function(object, exposure) {
#   cat("*************************************************************************************************** \n")
#   cat("This function produces a table of effect estimates of all (some) network-variables on an outcome    \n")
#   cat("(time-to-event) for all possibles causal models based on conditional independence criteria encoded  \n")
#   cat("in the input-network => MULTISET OF POSSIBLE EFFECTS PER EXPOSURE-OUTCOME PAIR; network-variables   \n")
#   cat("of interest are selected by indicating the variable-names as character vector                       \n")
#   cat("***************************************************************************************************")
#
#   # Define empty dataframes
#   mm_coef <- structure(list(
#     Model = as.character(),
#     Nbmds = as.numeric(),
#     Exposure = as.character(),
#     Covariables = as.character(),
#     NCov = as.numeric(),
#     HR = as.numeric(),
#     LCL = as.numeric(),
#     UCL = as.numeric(),
#     Beta = as.numeric(),
#     rSE = as.numeric(),
#     P = as.numeric()
#   ),
#   class = "data.frame"
#   )
#   mm_coef_temp <- structure(list(
#     Model = as.character(),
#     Nbmds = as.numeric(),
#     Exposure = as.character(),
#     Covariables = as.character(),
#     NCov = as.numeric(),
#     HR = as.numeric(),
#     LCL = as.numeric(),
#     UCL = as.numeric(),
#     Beta = as.numeric(),
#     rSE = as.numeric(),
#     P = as.numeric()
#   ),
#   class = "data.frame"
#   )
#   # For each exposure-outcome pair: Get Exposure coefficients from all possible models and write to a table
#   for (j in exposure)
#   {
#     mm_coef_temp <- data.frame(getExp.coef.perexposure(object = object, exposure = j))
#     mm_coef <- bind_rows(mm_coef, mm_coef_temp)
#   }
#   # Merge result tables for all exposure-outcome pairs
#   mm_coef
# }
#
# # Extract exposure coefficients per outcome, i.e. metabolite#
#
# getExp.coef.peroutcome <- function(object, outcome, exposure) {
#
#   # object: output of net.coupler
#   # outcome: specific outcome, i.e. metabolite, to evaluate
#   # exposure: name of considered exposure variable
#
#   # create vector containing integers from 1 to number of outcome-specific models:
#   nbm <- c(1:length(object[[outcome]]$Model_summaries))
#   Nbmds <- max(nbm)
#
#   # define empty dataframes for output:
#   mm_coef_temp1 <- structure(list(
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
#
#   mm_coef_temp2 <- structure(list(
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
#
#   # loop along number of outcome(metabolite)-specific models:
#   for (i in seq(along = nbm)) {
#     # get exposure-effect estimates (beta coefficient, SE, t-value, p-value) from single model and write into dataframe
#     SUM <- summary(object[[outcome]]$Model_summaries[[i]]$Model_summary)
#     Cov <- rownames(SUM$coefficients)
#     # Cov<-dplyr::filter(data.frame(row.names(SUM$coefficients)),row.names(SUM$coefficients)!=exposure)                     #record specific covariates of this model
#
#     mm_coef_temp2 <- data.frame(
#       Model = as.character(object[[outcome]]$Model_summaries[[i]]$Model),
#       Nbmds = Nbmds,
#       Outcome = as.character(object[[outcome]]$Outcome_metabolite),
#       Exposure = exposure,
#       Covariables = as.character(paste(Cov, collapse = ", ")),
#       Estimate = as.numeric(SUM$coefficients[paste(exposure), 1]),
#       SE = as.numeric(SUM$coefficients[paste(exposure), 2]),
#       tval = as.numeric(SUM$coefficients[paste(exposure), 3]),
#       P = as.numeric(SUM$coefficients[paste(exposure), 4])
#     )
#
#     # bind information to an outcome-specific dataframe:
#     mm_coef_temp1 <- bind_rows(mm_coef_temp1, mm_coef_temp2)
#   }
#
#   return(mm_coef_temp1)
# }
#
# # Extract exposure coefficients per outcome, i.e. metabolite#
#
# # old version? there are a lot of differences between this function and the other three
# getExp.coef.peroutcome_old <- function(object, outcome, exposure) {
#
#   # Create a vector containing integers from 1 to number of outcome-specific models
#   nbm <- c(1:length(object$Outcomes[[outcome]]$Model_summaries))
#   Nbmds <- max(nbm)
#   exposure <- exposure
#   # Create emty data-frames for the output
#   mm_coef_temp1 <- structure(list(
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
#   mm_coef_temp2 <- structure(list(
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
#
#   # loop along number of outcome-specific models
#   for (i in seq(along = nbm))
#   { # get exposure-effect estimates (beta coefficient, SE, t-value, p-value) from single model and write into a dataframe
#     SUM <- summary(object$Outcomes[[outcome]]$Model_summaries[[i]]$Model_summary)
#     Cov <- (dplyr::filter(data.frame(row.names(SUM$coefficients)), row.names(SUM$coefficients) != "exp"))
#
#     mm_coef_temp2 <- data.frame(
#       Model = as.character(object$Outcomes[[outcome]]$Model_summaries[[i]]$Model),
#       Nbmds = Nbmds,
#       Outcome = as.character(object$Outcomes[[outcome]]$Outcome),
#       Exposure = exposure,
#       Covariables = as.character(paste(Cov[, 1], collapse = ", ")),
#       Estimate = as.numeric(SUM$coefficients[paste(exposure), 1]),
#       SE = as.numeric(SUM$coefficients[paste(exposure), 2]),
#       tval = as.numeric(SUM$coefficients[paste(exposure), 3]),
#       P = as.numeric(SUM$coefficients[paste(exposure), 4])
#     )
#
#     # bind information to a outcome-specific dataframe
#     mm_coef_temp1 <- bind_rows(mm_coef_temp1, mm_coef_temp2)
#   }
#   mm_coef_temp1
# }
#
# # Extract exposure coefficients per outcome, i.e. metabolite#
#
# # getExp.coef.permetabolite <- function(object, metabolite) {
# #
# #   # object: output of net.coupler.out
# #   # metabolite: specific metabolite to evaluate
# #
# #   # create vector containing integers from 1 to number of metabolite-specific models:
# #   nbm <- c(1:length(object$model_details_all[[metabolite]]$Model_summaries))
# #   Nbmds <- max(nbm)
# #
# #   # define empty dataframes for output:
# #   mm_coef_temp1 <- structure(list(
# #     Model = as.character(),
# #     Nbmds = as.numeric(),
# #     Outcome = as.character(),
# #     Covariables = as.character(),
# #     # NCov=as.numeric(),
# #     HR = as.numeric(),
# #     LCL = as.numeric(),
# #     UCL = as.numeric(),
# #     Beta = as.numeric(),
# #     rSE = as.numeric(),
# #     P = as.numeric()
# #   ),
# #   class = "data.frame"
# #   )
# #
# #   mm_coef_temp2 <- structure(list(
# #     Model = as.character(),
# #     Nbmds = as.numeric(),
# #     Outcome = as.character(),
# #     Covariables = as.character(),
# #     NCov = as.numeric(),
# #     HR = as.numeric(),
# #     LCL = as.numeric(),
# #     UCL = as.numeric(),
# #     Beta = as.numeric(),
# #     rSE = as.numeric(),
# #     P = as.numeric()
# #   ),
# #   class = "data.frame"
# #   )
# #
# #   # loop along number of metabolite-specific models:
# #   for (i in seq(along = nbm)) {
# #     # get metabolite-effect estimates (hazard-ratio, lower and upper confidence limit, beta, robust SE, and p-value) from single model and write into dataframe
# #     SUM <- summary(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model_summary)
# #     Cov <- rownames(SUM$coefficients)
# #     # Cov<-dplyr::filter(data.frame(row.names(SUM$coefficients)),row.names(SUM$coefficients)!=exposure)                     #record specific covariates of this model
# #
# #     mm_coef_temp2 <- data.frame(
# #       Model = as.character(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model),
# #       Nbmds = Nbmds,
# #       Outcome = metabolite,
# #       Covariables = as.character(paste(Cov, collapse = ", ")),
# #       NCov = max(unlist(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model_summary$assign)) - 1,
# #       HR = as.numeric(SUM$conf.int[[metabolite, 1]]),
# #       LCL = as.numeric(SUM$conf.int[[metabolite, 3]]),
# #       UCL = as.numeric(SUM$conf.int[[metabolite, 4]]),
# #       Beta = as.numeric(SUM$coefficients[[metabolite, 1]]),
# #       rSE = as.numeric(SUM$coefficients[[metabolite, 4]]),
# #       P = as.numeric(SUM$coefficients[[metabolite, 6]])
# #     )
# #
# #     # bind information to an outcome-specific dataframe:
# #     mm_coef_temp1 <- bind_rows(mm_coef_temp1, mm_coef_temp2)
# #   }
# #
# #   return(mm_coef_temp1)
# # }
#
