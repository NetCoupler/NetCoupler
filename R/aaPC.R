#
# NETcoupler.Out 1
#

memory.limit(size = 250000)
# PC
a <- data.frame(l = letters[1:26])
b <- bind_rows(data.frame(l = letters[2:26]), data.frame(l = letters[1]))
c <- bind_rows(data.frame(l = letters[3:26]), data.frame(l = letters[1:2]))
d <- bind_rows(data.frame(l = letters[4:26]), data.frame(l = letters[1:3]))
A1 <- data.frame(Let = LETTERS[1:26])
A2 <- A1
A3 <- A1
A4 <- A1
LA <- dplyr::bind_rows(A1, A2, A3, A4)
l_abcd <- dplyr::bind_rows(a, b, c, d)

Ll <- dplyr::bind_cols(LA, l_abcd) %>% dplyr::transmute(Ll = paste0(Let, l))

load("H:/Metabolomics/DISS I/R_obj/T2D/T2D_data")
# load("H:/Metabolomics/DISS I/R_obj/T2D/STR_GPL")
GPL_data <- readRDS("H:/Metabolomics/DISS I/R_obj/T2D/STR_GPL.rds")
GPL_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/GPL/STR_GPL_SC.rds")
EXP_data <- readRDS("H:/Metabolomics/DISS I/R_obj/T2D/EXP_DATA.rds")
PC_data <- dplyr::select(GPL_data, contains("_aa_"))
PC_data_SC <- dplyr::select(GPL_data_SC, contains("_aa_"))

#
# Using Gaussian Data
#
# Load predefined data
n <- nrow(PC_data_SC)
V <- colnames(data.frame(PC_data_SC)) # labels aka node names
# estimate CPDAG
# PC_skel <-  skeleton(suffStat = list(C = cor(PC_data), n = n),
#                     indepTest = gaussCItest,labels = V, method ="stable", # indep.test: partial correlations
#                     alpha=0.05, fixedGaps = NULL, fixedEdges = NULL,verbose = FALSE)
# saveRDS(PC_skel,"H:/Metabolomics/DISS I/R_Obj/T2D/PC_skel.rds")

# PC_DAG <-  pc(suffStat = list(C = cor(PC_data_SC), n = n),
#              indepTest = gaussCItest,labels = V, skel.method ="stable", # indep.test: partial correlations
#              alpha=0.05, fixedGaps = NULL, fixedEdges = NULL,verbose = FALSE,maj.rule = F, solve.confl = F)
# saveRDS(PC_DAG,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/PC_DAG.rds")

PC_ren <- c(Ll$Ll[1:length(colnames(PC_data_SC))])
RENAME_PC <- dplyr::bind_cols(data.frame(colnames(PC_data_SC)), data.frame(PC_ren))
colnames(RENAME_PC) <- c("Metabolite", "Exposure")
colnames(PC_data_SC) <- PC_ren
colnames(PC_data) <- PC_ren

LABEL <- data.frame(Metabolite = RENAME_PC$Metabolite, Met_label = sapply(strsplit(as.character(RENAME_PC$Metabolite), split = "_", fixed = TRUE), function(x) (paste0(x[3], sep = "/", x[4]))))
LABEL$Metabolite <- as.character(LABEL$Metabolite)

# T2D_data <- data.frame(T2D_data)
# SURV<-data.frame(SURV)
PCT2D_data <- dplyr::bind_cols(PC_data, T2D_data, EXP_data)
SURV <- Surv(T2D_data$sta_time, T2D_data$sto_time, T2D_data$fall)
print(is.Surv(SURV))
#
# Using Gaussian Data
#
# Load predefined data
n <- nrow(PC_data_SC)
V <- colnames(data.frame(PC_data_SC)) # labels aka node names
PC_skel <- skeleton(
  suffStat = list(C = cor(PC_data_SC), n = n),
  indepTest = gaussCItest, labels = V, method = "stable", # indep.test: partial correlations
  alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE
)
# saveRDS(PC_skel,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/PC_skel.rds")

#
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

#
#
# 2.          MULTIMODEL INFERENCE

PC_OUT <- net.coupler.out(Graph = PC_skel, data = PCT2D_data)
# save(PC_OUT,file="D:/R_work/Clemens/PC/PC_OUT")
PC_T2D <- getExp.coef.out(object = PC_OUT, exposure = V)
PC_T2D$Exposure <- as.character(PC_T2D$Exposure)
RENAME_PC$Exposure <- as.character(RENAME_PC$Exposure)
PC_T2D <- dplyr::full_join(PC_T2D, RENAME_PC, by = "Exposure")
PC_T2D$Metabolite <- as.character(PC_T2D$Metabolite)

PC_T2Dsum <- left_join(PC_T2D %>% dplyr::group_by(Exposure) %>% dplyr::summarise(Nbmds = max(Nbmds), avgHR = mean(HR), minHR = min(HR), maxHR = max(HR), upperP = max(P), lowerP = min(P)),
  PC_T2D %>% dplyr::group_by(Exposure) %>% dplyr::filter(nchar(Covariables) == max(nchar(Covariables))) %>% dplyr::select(NCov = NCov, Exposure = Exposure, Metabolite = Metabolite, ind_HR = HR, ind_P = P)
  ,
  by = "Exposure"
)

#
# Adjust p-Values (controlling false discovery rate):
p.adjust.M <- p.adjust.methods[p.adjust.methods == "fdr"]
p.adj <- sapply(p.adjust.M, function(meth) p.adjust(PC_T2Dsum$ind_P, meth))
# fdr<-noquote(apply(p.adj, 2, format.pval, digits = 3))
PC_T2Dsum <- bind_cols(PC_T2Dsum, data.frame(p.adj))
rm(p.adj)
PC_T2Dsum <- PC_T2Dsum %>% dplyr::mutate(avgEst = log(avgHR), lowEst = log(minHR), highEst = log(maxHR), bestGuess = log(ind_HR))
PC_T2Dsum <- mutate(PC_T2Dsum,
  DE = derivedFactor(
    "1" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst > 0),
    "2" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst < 0),
    .method = "first",
    .default = 0
  )
)
# marginally associated (Assoc) if fdr-controlled p-Value below threshold
# "1" if positive betas, "2" if negative betas
PC_T2Dsum <- mutate(PC_T2Dsum,
  Assoc = derivedFactor(
    "1" = (fdr < 0.1 & bestGuess > 0),
    "2" = (fdr < 0.1 & bestGuess < 0),
    .method = "first",
    .default = 0
  )
)
# saveRDS(PC_T2D,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/PC_T2D.rds")
# saveRDS(PC_T2Dsum,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/PC_T2Dsum.rds")
# rm(PC_OUT)

AMB1_PC_Diab <- PC_T2Dsum %>% filter(Assoc != 0 & DE == 0)
DE1_PC_Diab <- PC_T2Dsum %>% filter(Assoc != 0 & DE != 0)

#
#
# NETcoupler.In 1
#
#
# Exposure
#
#
# PC
memory.limit(size = 250000)
Ll <- dplyr::bind_cols(LA, l_abcd) %>% dplyr::transmute(Ll = paste0(Let, l))
a <- data.frame(l = letters[1:26])
b <- dplyr::bind_rows(data.frame(l = letters[2:26]), data.frame(l = letters[1]))
c <- dplyr::bind_rows(data.frame(l = letters[3:26]), data.frame(l = letters[1:2]))
d <- dplyr::bind_rows(data.frame(l = letters[4:26]), data.frame(l = letters[1:3]))
A1 <- data.frame(Let = LETTERS[1:26])
A2 <- A1
A3 <- A1
A4 <- A1
LA <- dplyr::bind_rows(A1, A2, A3, A4)
l_abcd <- dplyr::bind_rows(a, b, c, d)

Ll <- dplyr::bind_cols(LA, l_abcd)
Ll <- Ll %>% dplyr::transmute(Ll = paste0(Let, l))
rm(a, b, c, d, A1, A2, A3, A4, LA, l_abcd)
GPL_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/GPL/STR_GPL_SC.rds")
Exp_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/GPL/EXP_SC.rds")
PC_data_SC <- dplyr::select(GPL_data_SC, contains("_aa_"))
#
# Using Gaussian Data
#
# Load predefined data
PC_ren_SC <- c(Ll$Ll[1:length(colnames(PC_data_SC))])
RENAME_PC_SC <- dplyr::bind_cols(data.frame(colnames(PC_data_SC)), data.frame(PC_ren_SC))
colnames(RENAME_PC_SC) <- c("Metabolite", "Outcome")
colnames(PC_data_SC) <- PC_ren_SC

n <- nrow(PC_data_SC)
V <- colnames(data.frame(PC_data_SC)) # labels aka node names

# estimate CPDAG
PC_skel <- skeleton(
  suffStat = list(C = cor(PC_data_SC), n = n),
  indepTest = gaussCItest, labels = V, method = "stable", # indep.test: partial correlations
  alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE
)

PC_DAG <- pc(
  suffStat = list(C = cor(PC_data_SC), n = n),
  indepTest = gaussCItest, labels = V, skel.method = "stable", # indep.test: partial correlations
  alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, maj.rule = F, solve.confl = F
)

# saveRDS(PC_skel,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/aaPC/FV/Results/Networks/PC_skel.rds")
# saveRDS(PC_DAG,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/aaPC/FV/Results/Networks/PC_DAG.rds")

getExp.coef.peroutcome <- function(object, outcome, exposure) {

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

# B. Get exposure coefficients for a a group of outcomes: e.g. all network-variables
getExp.coef <- function(object, outcome, exposure) {
  cat("*************************************************************************************************** \n")
  cat("This function produces a table of effect estimates of the exposure (betas, SE, t-value, p-value) on \n")
  cat("a selected (set of) outcome(s) [specified as character-vector of variable-names] for all possible   \n")
  cat("causal models based on conditional independence criteria encoded in the input-network => MULTISET OF\n")
  cat("POSSIBLE EFFECTS PER EXPOSURE-OUTCOME PAIRS: NUMBER OF POSSIBLE EFFECTS ~ SIZE OF THE ADJECENCY SET!\n")
  cat("***************************************************************************************************")
  exposure <- exposure
  # Define empty dataframes
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
  # For each exposure-outcome pair: Get Exposure coefficients from all possible models and write to a table
  for (j in outcome)
  {
    mm_coef_temp <- data.frame(getExp.coef.peroutcome(object = object, outcome = j, exposure = exposure))
    mm_coef <- bind_rows(mm_coef, mm_coef_temp)
  }
  # Merge result tables for all exposure-outcome pairs
  mm_coef
}
# 1.4.1 Round numeric columns in dataframes
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame
  # digits: number of digits to round
  numeric_columns <- sapply(x, class) == "numeric"
  x[numeric_columns] <- round(x[numeric_columns], digits)
  x
}

# 1.4.2 get object-names as string-variables
name.as.string <- function(v1) {
  deparse(substitute(v1))
}

MinMod <- c(paste(colnames(Exp_data_SC), collapse = ", "))
# saveRDS(MinMod,"Covs.rds")
print(MinMod)

#
# WGB
PC_IN <- net.coupler(Graph = PC_skel, Graph_data = PC_data_SC, Exp_data_SC = Exp_data_SC, exposure = "WGBperMJ")
SC_PC_WGB <- getExp.coef(object = PC_IN, outcome = colnames(data.frame(PC_data_SC)), exposure = "WGBperMJ")
# saveRDS(PC_IN,"D:/R_work/Clemens/PC/MULTI_PC_IN.rds")
RENAME_PC_SC$Outcome <- as.character(RENAME_PC_SC$Outcome)
SC_PC_WGB <- merge(SC_PC_WGB, RENAME_PC_SC, by = "Outcome")
SC_PCsum_WGB <- dplyr::left_join(SC_PC_WGB %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgbeta = mean(Estimate), minbeta = min(Estimate), maxbeta = max(Estimate), upperP = max(P), lowerP = min(P), Nbmds = mean(Nbmds)),
  SC_PC_WGB %>% dplyr::filter(Covariables == MinMod) %>% dplyr::select(Outcome = Outcome, Metabolite = Metabolite, marg_Est = Estimate, marg_P = P)
  ,
  by = "Outcome"
)
p.adjust.M <- p.adjust.methods[p.adjust.methods == "fdr"]
p.adj <- sapply(p.adjust.M, function(meth) p.adjust(SC_PCsum_WGB$marg_P, meth))
# fdr<-noquote(apply(p.adj, 2, format.pval, digits = 3))
SC_PCsum_WGB <- bind_cols(SC_PCsum_WGB, data.frame(p.adj))
rm(p.adj)
SC_PCsum_WGB <- SC_PCsum_WGB %>% dplyr::rename(avgEst = avgbeta, lowEst = minbeta, highEst = maxbeta, bestGuess = marg_Est)

SC_PCsum_WGB <- mutate(SC_PCsum_WGB,
  DE = derivedFactor(
    "1" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst > 0),
    "2" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst < 0),
    .method = "first",
    .default = 0
  )
)
# marginally associated (Assoc) if fdr-controlled p-Value below threshold
# "1" if positive betas, "2" if negative betas
SC_PCsum_WGB <- mutate(SC_PCsum_WGB,
  Assoc = derivedFactor(
    "1" = (fdr < 0.1 & bestGuess > 0),
    "2" = (fdr < 0.1 & bestGuess < 0),
    .method = "first",
    .default = 0
  )
)

AMB1_PC_WGB <- SC_PCsum_WGB %>% filter(Assoc != 0 & DE == 0)
DE1_PC_WGB <- SC_PCsum_WGB %>% filter(Assoc != 0 & DE != 0)

SC_PCsum_WGB_FV <- SC_PCsum_WGB %>% dplyr::mutate(round = 1, Assoc_FV = as.character(Assoc), DE_FV = as.character(DE))
SC_PCsum_WGB_FV <- dplyr::left_join(SC_PCsum_WGB_FV, LABEL, by = "Metabolite")
SC_PCsum_WGB_FV$Met_label <- as.character(SC_PCsum_WGB_FV$Met_label)
SC_PCsum_WGB_FV <- SC_PCsum_WGB_FV %>% dplyr::rowwise() %>% dplyr::mutate(Est_range_FV = paste0(formatC(bestGuess, format = "f", digits = 2), " (", formatC(lowEst, format = "f", digits = 2), ", ", formatC(highEst, format = "f", digits = 2), ")"))
# WGBCC<-CC.WGB.PC%>%dplyr::select(Metabolite,contains("CC"))
# SC_PCsum_WGB_FV<-dplyr::left_join(SC_PCsum_WGB_FV,WGBCC,by="Metabolite")

#
# Redmeat
# PC_IN<-net.coupler(Graph=PC_skel , Graph_data=PC_data_SC , Exp_data_SC=Exp_data_SC ,exposure="TMperMJ")
SC_PC_Redmeat <- getExp.coef(object = PC_IN, outcome = colnames(data.frame(PC_data_SC)), exposure = "TMperMJ")
# saveRDS(PC_IN,"D:/R_work/Clemens/PC/MULTI_PC_IN.rds")
RENAME_PC_SC$Outcome <- as.character(RENAME_PC_SC$Outcome)
SC_PC_Redmeat <- merge(SC_PC_Redmeat, RENAME_PC_SC, by = "Outcome")
SC_PCsum_Redmeat <- dplyr::left_join(SC_PC_Redmeat %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgbeta = mean(Estimate), minbeta = min(Estimate), maxbeta = max(Estimate), upperP = max(P), lowerP = min(P), Nbmds = mean(Nbmds)),
  SC_PC_Redmeat %>% dplyr::filter(Covariables == MinMod) %>% dplyr::select(Outcome = Outcome, Metabolite = Metabolite, marg_Est = Estimate, marg_P = P)
  ,
  by = "Outcome"
)
p.adjust.M <- p.adjust.methods[p.adjust.methods == "fdr"]
p.adj <- sapply(p.adjust.M, function(meth) p.adjust(SC_PCsum_Redmeat$marg_P, meth))
# fdr<-noquote(apply(p.adj, 2, format.pval, digits = 3))
SC_PCsum_Redmeat <- bind_cols(SC_PCsum_Redmeat, data.frame(p.adj))
rm(p.adj)
SC_PCsum_Redmeat <- SC_PCsum_Redmeat %>% dplyr::rename(avgEst = avgbeta, lowEst = minbeta, highEst = maxbeta, bestGuess = marg_Est)

# directly exposure-affected (DE) if upper P from all models > threshold AND all betas in same direction
# "1" if positive betas, "2" if negative betas
SC_PCsum_Redmeat <- mutate(SC_PCsum_Redmeat,
  DE = derivedFactor(
    "1" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst > 0),
    "2" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst < 0),
    .method = "first",
    .default = 0
  )
)
# marginally associated (Assoc) if fdr-controlled p-Value below threshold
# "1" if positive betas, "2" if negative betas
SC_PCsum_Redmeat <- mutate(SC_PCsum_Redmeat,
  Assoc = derivedFactor(
    "1" = (fdr < 0.1 & bestGuess > 0),
    "2" = (fdr < 0.1 & bestGuess < 0),
    .method = "first",
    .default = 0
  )
)

AMB1_PC_Redmeat <- SC_PCsum_Redmeat %>% filter(Assoc != 0 & DE == 0)
DE1_PC_Redmeat <- SC_PCsum_Redmeat %>% filter(Assoc != 0 & DE != 0)

# SC_PCsum_Redmeat_FV<- SC_PCsum_Redmeat%>%dplyr::mutate(round=1,Assoc_FV=as.character(Assoc),DE_FV=as.character(DE))

#
# Coffee
# PC_IN<-net.coupler(Graph=PC_skel , Graph_data=PC_data_SC , Exp_data_SC=Exp_data_SC ,exposure="CofCup")
SC_PC_Coffee <- getExp.coef(object = PC_IN, outcome = colnames(data.frame(PC_data_SC)), exposure = "CofCup")
# saveRDS(PC_IN,"D:/R_work/Clemens/PC/MULTI_PC_IN.rds")
RENAME_PC_SC$Outcome <- as.character(RENAME_PC_SC$Outcome)
SC_PC_Coffee <- merge(SC_PC_Coffee, RENAME_PC_SC, by = "Outcome")
SC_PCsum_Coffee <- dplyr::left_join(SC_PC_Coffee %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgbeta = mean(Estimate), minbeta = min(Estimate), maxbeta = max(Estimate), upperP = max(P), lowerP = min(P), Nbmds = mean(Nbmds)),
  SC_PC_Coffee %>% dplyr::filter(Covariables == MinMod) %>% dplyr::select(Outcome = Outcome, Metabolite = Metabolite, marg_Est = Estimate, marg_P = P)
  ,
  by = "Outcome"
)
p.adjust.M <- p.adjust.methods[p.adjust.methods == "fdr"]
p.adj <- sapply(p.adjust.M, function(meth) p.adjust(SC_PCsum_Coffee$marg_P, meth))
# fdr<-noquote(apply(p.adj, 2, format.pval, digits = 3))
SC_PCsum_Coffee <- bind_cols(SC_PCsum_Coffee, data.frame(p.adj))
rm(p.adj)
SC_PCsum_Coffee <- SC_PCsum_Coffee %>% dplyr::rename(avgEst = avgbeta, lowEst = minbeta, highEst = maxbeta, bestGuess = marg_Est)

# directly exposure-affected (DE) if upper P from all models > threshold AND all betas in same direction
# "1" if positive betas, "2" if negative betas
SC_PCsum_Coffee <- mutate(SC_PCsum_Coffee,
  DE = derivedFactor(
    "1" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst > 0),
    "2" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst < 0),
    .method = "first",
    .default = 0
  )
)
# marginally associated (Assoc) if fdr-controlled p-Value below threshold
# "1" if positive betas, "2" if negative betas
SC_PCsum_Coffee <- mutate(SC_PCsum_Coffee,
  Assoc = derivedFactor(
    "1" = (fdr < 0.1 & bestGuess > 0),
    "2" = (fdr < 0.1 & bestGuess < 0),
    .method = "first",
    .default = 0
  )
)

AMB1_PC_Coffee <- SC_PCsum_Coffee %>% filter(Assoc != 0 & DE == 0)
DE1_PC_Coffee <- SC_PCsum_Coffee %>% filter(Assoc != 0 & DE != 0)

SC_PCsum_Coffee_FV <- SC_PCsum_Coffee %>% dplyr::mutate(round = 1, Assoc_FV = as.character(Assoc), DE_FV = as.character(DE))
SC_PCsum_Coffee_FV <- dplyr::left_join(SC_PCsum_Coffee_FV, LABEL, by = "Metabolite")
SC_PCsum_Coffee_FV$Met_label <- as.character(SC_PCsum_Coffee_FV$Met_label)
SC_PCsum_Coffee_FV <- SC_PCsum_Coffee_FV %>% dplyr::rowwise() %>% dplyr::mutate(Est_range_FV = paste0(formatC(bestGuess, format = "f", digits = 2), " (", formatC(lowEst, format = "f", digits = 2), ", ", formatC(highEst, format = "f", digits = 2), ")"))
# CoffeeCC<-CC.Coffee.PC%>%dplyr::select(Metabolite,contains("CC"))
# SC_PCsum_Coffee_FV<-dplyr::left_join(SC_PCsum_Coffee_FV,CoffeeCC,by="Metabolite")

rm(PC_IN)

#
#
# Generate ConnectedComponents
#
#
# deleteAllWindows(CytoscapeConnection())
memory.limit(size = 250000)
#
# load("H:/Metabolomics/DISS I/R_obj/PC/MULTI_PC_IN")
# show(CHECK_IN)

#
# define function to get rownames as variable in dplyr
draw_rownames <- function(.data) .data %>% do(mutate(., Metabolite = rownames(.)))
draw_rownames_out <- function(.data) .data %>% do(mutate(., Metabolite = as.factor(rownames(.))))

#
# import multi-model estimates for Exposure effects on Acylcarnitines:
GPL_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/GPL/STR_GPL_SC.rds")
PC_data_SC <- dplyr::select(GPL_data_SC, contains("_aa_"))
pCor_PC <- ppcor::pcor(PC_data_SC)
colnames(pCor_PC$estimate) <- colnames(PC_data_SC)
rownames(pCor_PC$estimate) <- colnames(PC_data_SC)
is.matrix(pCor_PC$estimate)

#
# Generate networks: DAG; Skeleton; Adjacency matrix
n <- nrow(PC_data_SC)
V <- colnames(data.frame(PC_data_SC)) # labels aka node names
PC_skel <- skeleton(
  suffStat = list(C = cor(PC_data_SC), n = n),
  indepTest = gaussCItest, labels = V, method = "stable", # indep.test: partial correlations
  alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE
)

PC_DAG <- pc(
  suffStat = list(C = cor(PC_data_SC), n = n),
  indepTest = gaussCItest, labels = V, skel.method = "stable", # indep.test: partial correlations
  alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, maj.rule = F, solve.confl = F
)

# get adjacency matrix (transform GraphNEL in Igraph)
adjM_PC <- get.adjacency(graphNEL2igraph(PC_skel@graph))

#
# Write multi-model summaries into named list for looped data-processing:

Exp <- c("Redmeat", "WGB", "Coffee")
EXP_PC <- list(Redmeat = SC_PCsum_Redmeat, WGB = SC_PCsum_WGB, Coffee = SC_PCsum_Coffee)
names(EXP_PC)

Out <- c("Diabetes")
OUT_PC <- list(Diabetes = PC_T2Dsum)
names(OUT_PC)

#
# get estimated network + vector with metabolic-network node-names

NW_PC_EPIC <- PC_DAG@graph

NW_PC_EPIC <- initEdgeAttribute(NW_PC_EPIC,
  attribute.name = "pCor",
  attribute.type = "numeric",
  default.value = 1
)
NW_PC_EPIC <- initEdgeAttribute(NW_PC_EPIC,
  attribute.name = "Est_range",
  attribute.type = "char",
  default.value = 1
)

for (i in 1:length(NW_PC_EPIC@edgeData@data))
{
  Y <- unlist(strsplit(attributes(NW_PC_EPIC@edgeData@data)[[1]][[i]], "[|]"))
  print(Y)
  edgeData(NW_PC_EPIC, from = Y[[1]], to = Y[[2]], "pCor") <- pCor_PC$estimate[Y[[1]], Y[[2]]]
}
for (i in 1:length(NW_PC_EPIC@edgeData@data))
{
  Y <- unlist(strsplit(attributes(NW_PC_EPIC@edgeData@data)[[1]][[i]], "[|]"))
  print(Y)
  edgeData(NW_PC_EPIC, from = Y[[1]], to = Y[[2]], "Est_range") <- as.character(format(round(pCor_PC$estimate[Y[[1]], Y[[2]]], 2), nsmall = 2))
}
Mets <- NW_PC_EPIC@nodes

#
# GET CONNECTED COMPONENTS PER EXPOSURE
# Connected component(CC):
# CC Cluster(>=2) of metabolites where
# Met~Exp in non-metabolite adjusted model
NC.PC_res <- list()
nam <- c()
ALL_CC <- list()
# For each exposure do
for (k in Exp)
{
  print(k)
  EXP_PC[[k]]$Metabolite <- as.character(EXP_PC[[k]]$Metabolite)
  # delete edges from adjacency matrix whenever non-exposure associated node is part of the node-pair
  adjM2CoCo <- adjM_PC
  nonAssociated_Nodes <- c(unlist(EXP_PC[[k]] %>% dplyr::filter(Assoc == 0) %>% dplyr::select(Metabolite)))
  Associated_Nodes <- c(unlist(EXP_PC[[k]] %>% dplyr::filter(Assoc != 0) %>% dplyr::select(Metabolite)))
  adjM2CoCo[nonAssociated_Nodes, ] <- 0
  adjM2CoCo[, nonAssociated_Nodes] <- 0

  # extract connected components
  c <- clusters(graph_from_adjacency_matrix(adjM2CoCo))
  c$membership
  n_cc <- max(c$membership)
  CC_PC <- list()
  nam <- c()
  j <- c(0)
  # For each connected component (CC) assign values to Nodes
  # "True" if in CC, "False" otherwise
  for (i in 1:n_cc)
  {
    # who's in there?
    idx <- c$membership == i

    if (sum(idx) > 1) {
      print(i)
      j <- j + 1
      # print(j)
      # no singleton
      # do something with this one
      # checking algorithm
      # C = B[idx,idx]
      CC_PC[j] <- list(assign(paste0(k, "CC", sep = "_", j), idx)) # adment Output to  list
      nam <- c(nam, paste0(k, "CC", sep = "_", j)) # define name for list object
      # cat(idx)
      # cat('\n')
      # cat('=====\n')
    }
    names(CC_PC) <- nam # assign names to listed objects
    # transform into dataframe
    CC <- data.frame(CC_PC[1]) %>% draw_rownames() %>% dplyr::select(Metabolite)
    for (i in 1:length(CC_PC))
    {
      CC_temp <- data.frame()
      CC_temp <- data.frame(CC_PC[i]) %>% draw_rownames()
      CC_temp$Metabolite <- as.character(CC_temp$Metabolite)
      CC <- dplyr::left_join(CC, CC_temp, by = "Metabolite")
    }
    # Combine with multi_model output
    assign(paste0("CC", sep = ".", k, sep = ".", "PC"), dplyr::left_join(EXP_PC[[k]], CC, by = "Metabolite"))
    NC.PC_res[k] <- list(assign(paste0("CC", sep = ".", k, sep = ".", "PC"), dplyr::left_join(EXP_PC[[k]], CC, by = "Metabolite")))

    rm(CC)
  }
  # Keep a list of CC as well
  assign(paste0("CC", sep = "_", k, sep = "_", "PC"), CC_PC)
  ALL_CC <- append(ALL_CC, assign(paste0("CC", sep = "_", k, sep = "_", "PC"), CC_PC))
  rm(CC_PC)
}
#
# GET CONNECTED COMPONENTS PER Outcome
# Connected component(CC):
# CC Cluster(>=2) of metabolites where
# Met~Exp in non-metabolite adjusted model

# For each Outcome do
for (k in Out)
{
  # delete edges from adjacency matrix whenever non-exposure associated node is part of the node-pair
  adjM2CoCo <- adjM_PC
  nonAssociated_Nodes <- c(unlist(OUT_PC[[k]] %>% filter(Assoc == 0) %>% dplyr::select(Metabolite)))
  adjM2CoCo[nonAssociated_Nodes, ] <- 0
  adjM2CoCo[, nonAssociated_Nodes] <- 0

  # extract connected components
  c <- clusters(graph_from_adjacency_matrix(adjM2CoCo))
  c$membership
  n_cc <- max(c$membership)

  CC_PC <- list()
  nam <- c()
  j <- c(0)
  # For each connected component (CC) assign values to Nodes
  # "True" if in CC, "False" otherwise
  for (i in 1:n_cc)
  {
    # who's in there?
    idx <- c$membership == i
    # print(idx)

    if (sum(idx) > 1) {
      j <- j + 1
      # print(j)
      # no singleton
      # do something with this one
      # checking algorithm
      # C = B[idx,idx]
      CC_PC[j] <- list(assign(paste0(k, "CC", sep = "_", j), idx)) # adment Output to  list
      nam <- c(nam, paste0(k, "CC", sep = "_", j)) # define name for list object
      # cat(idx)
      # cat('\n')
      # cat('=====\n')
    }
    names(CC_PC) <- nam # assign names to listed objects
    # print(nam)
    # transform into dataframe
    CC <- data.frame(CC_PC[1]) %>% draw_rownames_out() %>% dplyr::select(Metabolite)
    CC$Metabolite <- as.character(CC$Metabolite)

    for (i in 1:length(CC_PC))
    {
      CC_temp <- data.frame()
      CC_temp <- data.frame(CC_PC[i]) %>% draw_rownames_out()
      CC_temp$Metabolite <- as.character(CC_temp$Metabolite)
      CC <- left_join(CC, CC_temp, by = "Metabolite")
    }
    # Combine with multi_model output
    assign(paste0("CC", sep = ".", k, sep = ".", "PC"), left_join(OUT_PC[[k]], CC, by = "Metabolite"))
    NC.PC_res[k] <- list(assign(paste0("CC", sep = ".", k, sep = ".", "PC"), left_join(OUT_PC[[k]], CC, by = "Metabolite")))
    rm(CC)
  }
  # Keep a list of CC as well
  assign(paste0("CC", sep = "_", k, sep = "_", "PC"), CC_PC)
  ALL_CC <- append(ALL_CC, assign(paste0("CC", sep = "_", k, sep = "_", "PC"), CC_PC))
  rm(CC_PC)
}

NC.PC <- list()

for (i in 1:length(NC.PC_res))
{
  NC.PC[[i]] <- NC.PC_res[[i]] %>% dplyr::mutate(AssEf = paste0(Assoc, sep = "_", DE))
}
nam <- names(NC.PC_res)
names(NC.PC) <- sprintf(nam, seq_along(nam))
CC.labels <- paste0("CC", sep = ".", k, sep = ".", "PC")

#
#
# NETcoupler.In 1
#
#
# Exposure
#

WGBCC <- CC.WGB.PC %>% dplyr::select(Metabolite, contains("CC"))
SC_PCsum_WGB_FV <- dplyr::left_join(SC_PCsum_WGB_FV, WGBCC, by = "Metabolite")
CoffeeCC <- CC.Coffee.PC %>% dplyr::select(Metabolite, contains("CC"))
SC_PCsum_Coffee_FV <- dplyr::left_join(SC_PCsum_Coffee_FV, CoffeeCC, by = "Metabolite")

# saveRDS(SC_PC_Coffee,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/aaPC/FV/Results/Coffee/SC_PC_Coffee.rds")
# saveRDS(SC_PCsum_Coffee,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/aaPC/FV/Results/Coffee/SC_PCsum_Coffee.rds")
# saveRDS(SC_PCsum_Coffee_FV,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/aaPC/FV/Results/Coffee/SC_PCsum_Coffee_FV.rds")
SC_PCsum_Coffee_FV <- SC_PCsum_Coffee_FV %>% dplyr::rowwise() %>% dplyr::mutate(Est_range_FV = paste0(formatC(bestGuess, format = "f", digits = 2), " (", formatC(lowEst, format = "f", digits = 2), ", ", formatC(highEst, format = "f", digits = 2), ")"))
rm(PC_IN)

#
# PC
memory.limit(size = 250000)
Ll <- dplyr::bind_cols(LA, l_abcd) %>% dplyr::transmute(Ll = paste0(Let, l))
a <- data.frame(l = letters[1:26])
b <- dplyr::bind_rows(data.frame(l = letters[2:26]), data.frame(l = letters[1]))
c <- dplyr::bind_rows(data.frame(l = letters[3:26]), data.frame(l = letters[1:2]))
d <- dplyr::bind_rows(data.frame(l = letters[4:26]), data.frame(l = letters[1:3]))
A1 <- data.frame(Let = LETTERS[1:26])
A2 <- A1
A3 <- A1
A4 <- A1
LA <- dplyr::bind_rows(A1, A2, A3, A4)
l_abcd <- dplyr::bind_rows(a, b, c, d)

Ll <- dplyr::bind_cols(LA, l_abcd)
Ll <- Ll %>% dplyr::transmute(Ll = paste0(Let, l))
rm(a, b, c, d, A1, A2, A3, A4, LA, l_abcd)
GPL_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/GPL/STR_GPL_SC.rds")
Exp_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/GPL/EXP_SC.rds")
PC_data_SC <- dplyr::select(GPL_data_SC, contains("_aa_"))
#
# Using Gaussian Data
#
# Load predefined data
PC_ren_SC <- c(Ll$Ll[1:length(colnames(PC_data_SC))])
RENAME_PC_SC <- dplyr::bind_cols(data.frame(colnames(PC_data_SC)), data.frame(PC_ren_SC))
colnames(RENAME_PC_SC) <- c("Metabolite", "Outcome")
colnames(PC_data_SC) <- PC_ren_SC

#
# deleteAllWindows(CytoscapeConnection())
memory.limit(size = 250000)
#
# load("H:/Metabolomics/DISS I/R_obj/PC/MULTI_PC_IN")
# show(CHECK_IN)
draw_rownames_Exp <- function(.data) .data %>% do(mutate(., Exp = rownames(.)))
draw_rownames <- function(.data) .data %>% do(mutate(., Outcome = rownames(.)))
draw_rownames_out <- function(.data) .data %>% do(mutate(., Metabolite = as.factor(rownames(.))))
#
# Write multi-model summaries into named list for looped data-processing:
# PC_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/PC/STR_PC_SC.rds")
# Exp_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/PC/EXP_SC.rds")
# Exp<-c("Redmeat","WGB","Coffee")
# EXP_PC<-list(Redmeat=SC_PCsum_Redmeat,WGB=SC_PCsum_WGB,Coffee=SC_PCsum_Coffee)
# names(EXP_PC)
LABEL <- data.frame(Metabolite = RENAME_PC$Metabolite, Met_label = sapply(strsplit(as.character(RENAME_PC$Metabolite), split = "_", fixed = TRUE), function(x) (paste0(x[3], sep = "/", x[4]))))
# PC_ren_SC<-c(Ll$Ll[1:length(colnames(PC_data_SC))])
# RENAME_PC_SC<-dplyr::bind_cols(data.frame(colnames(PC_data_SC)),data.frame(PC_ren_SC))
# colnames(RENAME_PC_SC)<-c("Metabolite","Outcome")
# colnames(PC_data_SC)<-PC_ren_SC
PC_skel <- readRDS("L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/PC_skel.rds")
# PREPARE AN EMPTY LISTS TO STORE DE and AMB per exposure

AMB_exp <- list()
DE_exp <- list()
for (i in Exp)
{
  AMB <- list(NULL)
  AMB_exp[[i]] <- AMB

  DE <- list(NULL)
  DE_exp[[i]] <- AMB
}

names(AMB_exp) <- Exp
names(DE_exp) <- Exp

for (i in Exp)
{
  AMB_exp[i] <- assign(paste0("AMB1_PC", sep = "_", i), EXP_PC[[i]] %>% dplyr::filter(Assoc != 0 & DE == 0))
  DE_exp[i] <- assign(paste0("DE1_PC", sep = "_", i), EXP_PC[[i]] %>% dplyr::filter(Assoc != 0 & DE != 0))
}

# CC<-data.frame()
# CC<-data.frame(CC_Redmeat_PC[2])%>%draw_rownames_out()%>%filter(CC_Redmeat_PC[[2]]==TRUE)%>%dplyr::select(Metabolite)
# CC$Metabolite<-as.character(CC$Metabolite)
# RENAME_PC_SC$Metabolite<-as.character(RENAME_PC_SC$Metabolite)
# RENAME_PC_SC$Outcome<-as.character(RENAME_PC_SC$Outcome)
# CC <-dplyr::inner_join(CC,RENAME_PC_SC,by="Metabolite")

# adj_plus<-CC%>%dplyr::filter(Outcome!=i)%>%dplyr::select(Outcome)
# MinMod<-paste0(paste0(adj_plus$Outcome, collapse=", "),", ",paste0(colnames(Exp_data_SC), collapse = ", "))
# print(MinMod)
# PC_IN_CC2<-net.coupler(Graph=PC_skel , Graph_data=PC_data_SC , Exp_data_SC=Exp_data_SC ,exposure="TMperMJ")
#
getExp.coef <- function(object, exposure) {
  SUM <- data.frame(NULL)
  for (j in 1:length(names(object$Outcomes)))
  {
    print(j)
    for (i in 1:PC_IN_CC2$Outcomes[[j]]$Number_of_Models)
    { # get exposure-effect estimates (beta coefficient, SE, t-value, p-value) from single model and write into a dataframe

      SUM1 <- data.frame(NULL)
      SUM1 <- data.frame(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]) %>%
        draw_rownames_Exp() %>%
        dplyr::filter(exposure == Exp) %>%
        dplyr::rename(SE = Std..Error, tval = t.value, P = Pr...t..) %>%
        dplyr::mutate(Outcome = outcome[[j]], Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
      SUM <- dplyr::bind_rows(SUM, SUM1)
    }
  }
  SUM
}

# 1.4.2 get object-names as string-variables
name.as.string <- function(v1) {
  deparse(substitute(v1))
}
#
#   Connected Component 1
#
EXP <- colnames(Exp_data_SC)
CC <- data.frame()
CC <- data.frame(CC_Redmeat_PC[1]) %>% draw_rownames_out() %>% filter(CC_Redmeat_PC[[1]] == TRUE) %>% dplyr::select(Metabolite)
CC$Metabolite <- as.character(CC$Metabolite)
RENAME_PC_SC$Metabolite <- as.character(RENAME_PC_SC$Metabolite)
RENAME_PC_SC$Outcome <- as.character(RENAME_PC_SC$Outcome)
CC <- dplyr::inner_join(CC, RENAME_PC_SC, by = "Metabolite")

DE1 <- intersect(CC$Outcome, DE_exp$Redmeat)
EXP <- colnames(Exp_data_SC)
MinMod <- paste0(EXP, collapse = ", ")
print(MinMod)
PC_IN_CC1 <- net.coupler(Graph = PC_skel, Graph_data = PC_data_SC, Exp_data_SC = Exp_data_SC, exposure = "TMperMJ")

getExp.coef <- function(object, exposure) {
  SUM <- data.frame(NULL)
  for (j in 1:length(names(object$Outcomes)))
  {
    print(j)
    for (i in 1:PC_IN_CC1$Outcomes[[j]]$Number_of_Models)
    { # get exposure-effect estimates (beta coefficient, SE, t-value, p-value) from single model and write into a dataframe

      SUM1 <- data.frame(NULL)
      SUM1 <- data.frame(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]) %>%
        draw_rownames_Exp() %>%
        dplyr::filter(exposure == Exp) %>%
        dplyr::rename(SE = Std..Error, tval = t.value, P = Pr...t..) %>%
        dplyr::mutate(Outcome = PC_IN_CC1$Outcomes[[j]]$Outcome, Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
      SUM <- dplyr::bind_rows(SUM, SUM1)
    }
  }
  SUM
}
SC_PC_Redmeat_CC1 <- getExp.coef(object = PC_IN_CC1, exposure = "TMperMJ")
# SC_PC_Redmeat$Outcome<-as.character(SC_PC_Redmeat$Outcome)

# RENAME_PC_SC$Metabolite<-as.character(RENAME_PC_SC$Metabolite)
# SC_PC_Redmeat<-merge(SC_PC_Redmeat, RENAME_PC_SC,by="Outcome")
# SC_PC_Redmeat$Metabolite<-as.character(SC_PC_Redmeat$Metabolite)
SC_PC_Redmeat_CC1 <- dplyr::inner_join(SC_PC_Redmeat_CC1, RENAME_PC_SC, by = "Outcome")

SC_PCsum_Redmeat_CC1 <- SC_PC_Redmeat_CC1 %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgbeta = mean(Estimate), minbeta = min(Estimate), maxbeta = max(Estimate), upperP = max(P), lowerP = min(P), Nbmds = mean(Nbmds), Metabolite = first(Metabolite))
SC_PCsum_Redmeat_CC1 <- dplyr::left_join(SC_PC_Redmeat_CC1 %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgbeta = mean(Estimate), minbeta = min(Estimate), maxbeta = max(Estimate), upperP = max(P), lowerP = min(P), Nbmds = mean(Nbmds), Metabolite = first(Metabolite)),
  SC_PC_Redmeat_CC1 %>% dplyr::filter(Covariables == MinMod) %>% dplyr::select(Outcome = first(Outcome), marg_Est = Estimate, marg_P = P)
  ,
  by = "Outcome"
)
p.adjust.M <- p.adjust.methods[p.adjust.methods == "fdr"]
p.adj <- sapply(p.adjust.M, function(meth) p.adjust(SC_PCsum_Redmeat_CC1$marg_P, meth))
# fdr<-noquote(apply(p.adj, 2, format.pval, digits = 3))
SC_PCsum_Redmeat_CC1 <- bind_cols(SC_PCsum_Redmeat_CC1, data.frame(p.adj))
rm(p.adj)
SC_PCsum_Redmeat_CC1 <- SC_PCsum_Redmeat_CC1 %>% dplyr::rename(avgEst = avgbeta, lowEst = minbeta, highEst = maxbeta, bestGuess = marg_Est)

SC_PCsum_Redmeat_CC1 <- mutate(SC_PCsum_Redmeat_CC1,
  DE = derivedFactor(
    "1" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst > 0),
    "2" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst < 0),
    .method = "first",
    .default = 0
  )
)
SC_PCsum_Redmeat_CC1 <- mutate(SC_PCsum_Redmeat_CC1,
  Assoc = derivedFactor(
    "1" = (fdr < 0.1 & bestGuess > 0),
    "2" = (fdr < 0.1 & bestGuess < 0),
    .method = "first",
    .default = 0
  )
)

rmPC_new_NEM <- filter(SC_PCsum_Redmeat_CC1, Assoc == 0)
rmPC_new_DE <- filter(SC_PCsum_Redmeat_CC1, Assoc != 0 & DE != 0)
rmPC_still_AMB <- filter(SC_PCsum_Redmeat_CC1, Assoc != 0 & DE == 0)
# saveRDS(SC_PC_Redmeat_CC1,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/SC_PC_Redmeat_CC1.rds")
# saveRDS(SC_PCsum_Redmeat_CC1,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/SC_PCsum_Redmeat_CC1.rds")

#
#   Connected Component 2
#
EXP <- colnames(Exp_data_SC)
CC <- data.frame()
CC <- data.frame(CC_Redmeat_PC[2]) %>% draw_rownames_out() %>% filter(CC_Redmeat_PC[[2]] == TRUE) %>% dplyr::select(Metabolite)
CC$Metabolite <- as.character(CC$Metabolite)
RENAME_PC_SC$Metabolite <- as.character(RENAME_PC_SC$Metabolite)
RENAME_PC_SC$Outcome <- as.character(RENAME_PC_SC$Outcome)
CC <- dplyr::inner_join(CC, RENAME_PC_SC, by = "Metabolite")

DE1 <- intersect(CC$Outcome, DE_exp$Redmeat)
EXP <- colnames(Exp_data_SC)
MinMod <- paste0(paste0(DE1, collapse = ", "), ", ", paste0(EXP, collapse = ", "))
print(MinMod)
PC_IN_CC2 <- net.coupler(Graph = PC_skel, Graph_data = PC_data_SC, Exp_data_SC = Exp_data_SC, exposure = "TMperMJ")
# saveRDS(PC_IN,"D:/R_work/Clemens/PC/MULTI_PC_IN.rds")
getExp.coef <- function(object, exposure) {
  SUM <- data.frame(NULL)
  for (j in 1:length(names(object$Outcomes)))
  {
    print(j)
    for (i in 1:PC_IN_CC2$Outcomes[[j]]$Number_of_Models)
    { # get exposure-effect estimates (beta coefficient, SE, t-value, p-value) from single model and write into a dataframe

      SUM1 <- data.frame(NULL)
      SUM1 <- data.frame(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]) %>%
        draw_rownames_Exp() %>%
        dplyr::filter(exposure == Exp) %>%
        dplyr::rename(SE = Std..Error, tval = t.value, P = Pr...t..) %>%
        dplyr::mutate(Outcome = PC_IN_CC2$Outcomes[[j]]$Outcome, Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
      SUM <- dplyr::bind_rows(SUM, SUM1)
    }
  }
  SUM
}
SC_PC_Redmeat_CC2 <- getExp.coef(object = PC_IN_CC2, exposure = "TMperMJ")
# SC_PC_Redmeat$Outcome<-as.character(SC_PC_Redmeat$Outcome)

# RENAME_PC_SC$Metabolite<-as.character(RENAME_PC_SC$Metabolite)
# SC_PC_Redmeat<-merge(SC_PC_Redmeat, RENAME_PC_SC,by="Outcome")
# SC_PC_Redmeat$Metabolite<-as.character(SC_PC_Redmeat$Metabolite)
SC_PC_Redmeat_CC2 <- dplyr::inner_join(SC_PC_Redmeat_CC2, RENAME_PC_SC, by = "Outcome")

SC_PCsum_Redmeat_CC2 <- SC_PC_Redmeat_CC2 %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgbeta = mean(Estimate), minbeta = min(Estimate), maxbeta = max(Estimate), upperP = max(P), lowerP = min(P), Nbmds = mean(Nbmds), Metabolite = first(Metabolite))
SC_PCsum_Redmeat_CC2 <- dplyr::left_join(SC_PC_Redmeat_CC2 %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgbeta = mean(Estimate), minbeta = min(Estimate), maxbeta = max(Estimate), upperP = max(P), lowerP = min(P), Nbmds = mean(Nbmds), Metabolite = first(Metabolite)),
  SC_PC_Redmeat_CC2 %>% dplyr::filter(Covariables == MinMod) %>% dplyr::select(Outcome = first(Outcome), marg_Est = Estimate, marg_P = P)
  ,
  by = "Outcome"
)
p.adjust.M <- p.adjust.methods[p.adjust.methods == "fdr"]
p.adj <- sapply(p.adjust.M, function(meth) p.adjust(SC_PCsum_Redmeat_CC2$marg_P, meth))
# fdr<-noquote(apply(p.adj, 2, format.pval, digits = 3))
SC_PCsum_Redmeat_CC2 <- bind_cols(SC_PCsum_Redmeat_CC2, data.frame(p.adj))
rm(p.adj)
SC_PCsum_Redmeat_CC2 <- SC_PCsum_Redmeat_CC2 %>% dplyr::rename(avgEst = avgbeta, lowEst = minbeta, highEst = maxbeta, bestGuess = marg_Est)

SC_PCsum_Redmeat_CC2 <- mutate(SC_PCsum_Redmeat_CC2,
  DE = derivedFactor(
    "1" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst > 0),
    "2" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst < 0),
    .method = "first",
    .default = 0
  )
)
SC_PCsum_Redmeat_CC2 <- mutate(SC_PCsum_Redmeat_CC2,
  Assoc = derivedFactor(
    "1" = (fdr < 0.1 & bestGuess > 0),
    "2" = (fdr < 0.1 & bestGuess < 0),
    .method = "first",
    .default = 0
  )
)

rmPC_new_NEM_temp <- filter(SC_PCsum_Redmeat_CC2, Assoc == 0)
rmPC_new_NEM <- bind_rows(rmPC_new_NEM, rmPC_new_NEM_temp)
rmPC_new_DE_temp <- filter(SC_PCsum_Redmeat_CC2, Assoc != 0 & DE != 0)
rmPC_new_DE <- bind_rows(rmPC_new_DE, rmPC_new_DE_temp)
rmPC_still_AMB_temp <- filter(SC_PCsum_Redmeat_CC2, Assoc != 0 & DE == 0)
rmPC_still_AMB <- bind_rows(rmPC_still_AMB, rmPC_still_AMB_temp)
# saveRDS(SC_PC_Redmeat_CC2,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/SC_PC_Redmeat_CC2.rds")
# saveRDS(SC_PCsum_Redmeat_CC2,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/SC_PCsum_Redmeat_CC2.rds")

#

#
#   Connected Component 3
#
EXP <- colnames(Exp_data_SC)
CC <- data.frame()
CC <- data.frame(CC_Redmeat_PC[3]) %>% draw_rownames_out() %>% filter(CC_Redmeat_PC[[3]] == TRUE) %>% dplyr::select(Metabolite)
CC$Metabolite <- as.character(CC$Metabolite)
RENAME_PC_SC$Metabolite <- as.character(RENAME_PC_SC$Metabolite)
RENAME_PC_SC$Outcome <- as.character(RENAME_PC_SC$Outcome)
CC <- dplyr::inner_join(CC, RENAME_PC_SC, by = "Metabolite")

DE1 <- intersect(CC$Outcome, DE_exp$Redmeat)
EXP <- colnames(Exp_data_SC)
MinMod <- paste0(paste0(DE1, collapse = ", "), ", ", paste0(EXP, collapse = ", "))
print(MinMod)
PC_IN_CC3 <- net.coupler(Graph = PC_skel, Graph_data = PC_data_SC, Exp_data_SC = Exp_data_SC, exposure = "TMperMJ")
# saveRDS(PC_IN,"D:/R_work/Clemens/PC/MULTI_PC_IN.rds")
getExp.coef <- function(object, exposure) {
  SUM <- data.frame(NULL)
  for (j in 1:length(names(object$Outcomes)))
  {
    print(j)
    for (i in 1:PC_IN_CC3$Outcomes[[j]]$Number_of_Models)
    { # get exposure-effect estimates (beta coefficient, SE, t-value, p-value) from single model and write into a dataframe

      SUM1 <- data.frame(NULL)
      SUM1 <- data.frame(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]) %>%
        draw_rownames_Exp() %>%
        dplyr::filter(exposure == Exp) %>%
        dplyr::rename(SE = Std..Error, tval = t.value, P = Pr...t..) %>%
        dplyr::mutate(Outcome = PC_IN_CC3$Outcomes[[j]]$Outcome, Nbmds = length(object$Outcomes[[j]]$Model_summaries), Exposure = exposure, Covariables = paste0(row.names(summary(object$Outcomes[[j]]$Model_summaries[[i]]$Model_summary)[["coefficients"]]), collapse = ", "))
      SUM <- dplyr::bind_rows(SUM, SUM1)
    }
  }
  SUM
}
SC_PC_Redmeat_CC3 <- getExp.coef(object = PC_IN_CC3, exposure = "TMperMJ")

SC_PC_Redmeat$Outcome <- as.character(SC_PC_Redmeat$Outcome)
SC_PC_Redmeat$Metabolite <- as.character(SC_PC_Redmeat$Metabolite)

# RENAME_PC_SC$Metabolite<-as.character(RENAME_PC_SC$Metabolite)
# SC_PC_Redmeat<-merge(SC_PC_Redmeat, RENAME_PC_SC,by="Outcome")
# SC_PC_Redmeat$Metabolite<-as.character(SC_PC_Redmeat$Metabolite)
SC_PC_Redmeat_CC3 <- dplyr::inner_join(SC_PC_Redmeat_CC3, RENAME_PC_SC, by = "Outcome")

SC_PC_Redmeat_FV <- dplyr::bind_rows(dplyr::select(dplyr::setdiff(SC_PC_Redmeat, dplyr::filter(SC_PC_Redmeat, Outcome %in% SC_PC_Redmeat_CC2$Outcome | Outcome %in% SC_PC_Redmeat_CC3$Outcome)), Metabolite, Outcome, Nbmds, Exposure, Covariables, Estimate, SE, P), dplyr::select(SC_PC_Redmeat_CC2, Metabolite, Outcome, Nbmds, Exposure, Covariables, Estimate, SE, P), dplyr::select(SC_PC_Redmeat_CC3, Metabolite, Outcome, Nbmds, Exposure, Covariables, Estimate, SE, P))
SC_PC_Redmeat_FV <- left_join(SC_PC_Redmeat_FV %>% group_by(Metabolite), LABEL, by = "Metabolite")
SC_PC_Redmeat_FV$Met_label <- as.character(SC_PC_Redmeat_FV$Met_label)
SC_PC_Redmeat_FV <- tbl_df(SC_PC_Redmeat_FV)

SC_PCsum_Redmeat_CC3 <- SC_PC_Redmeat_CC3 %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgbeta = mean(Estimate), minbeta = min(Estimate), maxbeta = max(Estimate), upperP = max(P), lowerP = min(P), Nbmds = mean(Nbmds), Metabolite = first(Metabolite))
SC_PCsum_Redmeat_CC3 <- dplyr::left_join(SC_PC_Redmeat_CC3 %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgbeta = mean(Estimate), minbeta = min(Estimate), maxbeta = max(Estimate), upperP = max(P), lowerP = min(P), Nbmds = mean(Nbmds), Metabolite = first(Metabolite)),
  SC_PC_Redmeat_CC3 %>% dplyr::filter(Covariables == MinMod) %>% dplyr::select(Outcome = first(Outcome), marg_Est = Estimate, marg_P = P)
  ,
  by = "Outcome"
)
# p.adjust.M <- p.adjust.methods[p.adjust.methods == "fdr"]
# p.adj    <- sapply(p.adjust.M, function(meth) p.adjust(SC_PCsum_Redmeat_CC3$marg_P, meth))
# fdr<-noquote(apply(p.adj, 2, format.pval, digits = 3))
SC_PCsum_Redmeat_CC3 <- SC_PCsum_Redmeat_CC3 %>% dplyr::mutate(fdr = marg_P)
# rm(p.adj )
SC_PCsum_Redmeat_CC3 <- SC_PCsum_Redmeat_CC3 %>% dplyr::rename(avgEst = avgbeta, lowEst = minbeta, highEst = maxbeta, bestGuess = marg_Est)

SC_PCsum_Redmeat_CC3 <- mutate(SC_PCsum_Redmeat_CC3,
  DE = derivedFactor(
    "1" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst > 0),
    "2" = (fdr < 0.1 & upperP < 0.05 & (sign(highEst) == sign(lowEst)) & lowEst < 0),
    .method = "first",
    .default = 0
  )
)
SC_PCsum_Redmeat_CC3 <- mutate(SC_PCsum_Redmeat_CC3,
  Assoc = derivedFactor(
    "1" = (fdr < 0.1 & bestGuess > 0),
    "2" = (fdr < 0.1 & bestGuess < 0),
    .method = "first",
    .default = 0
  )
)
rmPC_new_NEM_temp <- filter(SC_PCsum_Redmeat_CC3, Assoc == 0)
rmPC_new_NEM <- bind_rows(rmPC_new_NEM, rmPC_new_NEM_temp)
rmPC_new_DE_temp <- filter(SC_PCsum_Redmeat_CC3, Assoc != 0 & DE != 0)
rmPC_new_DE <- bind_rows(rmPC_new_DE, rmPC_new_DE_temp)
rmPC_still_AMB_temp <- filter(SC_PCsum_Redmeat_CC3, Assoc != 0 & DE == 0)
rmPC_still_AMB <- bind_rows(rmPC_still_AMB, rmPC_still_AMB_temp)
# saveRDS(SC_PC_Redmeat_CC3,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/SC_PC_Redmeat_CC3.rds")
# saveRDS(SC_PCsum_Redmeat_CC3,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/SC_PCsum_Redmeat_CC3.rds")
# rm(SC_PCsum_Redmeat_FV)
SC_PCsum_Redmeat_FV <- bind_rows(
  DE1_PC_Redmeat %>% dplyr::mutate(round = 1, Assoc_FV = as.character(Assoc), DE_FV = as.character(DE)),
  # rmPC_new_DE%>%dplyr::mutate(round=1,Assoc_FV=as.character(Assoc),DE_FV=as.character(DE)),
  rmPC_still_AMB %>% dplyr::mutate(round = 2, Assoc_FV = as.character(Assoc), DE_FV = as.character(DE)),
  SC_PCsum_Redmeat %>% filter(Metabolite %in% rmPC_new_NEM$Metabolite) %>% dplyr::mutate(round = 2, Assoc_FV = as.character(0), DE_FV = as.character(0)),
  SC_PCsum_Redmeat %>% filter(Assoc == 0) %>% dplyr::mutate(round = 1, Assoc_FV = as.character(Assoc), DE_FV = as.character(DE))
)

# RedmeatCC<-CC.Redmeat.PC%>%dplyr::select(Metabolite,contains("CC"))
# RedmeatCC<-dplyr::mutate(RedmeatCC,rmCCf = derivedFactor(
#  "1"= (RedmeatCC_1=="TRUE" ) ,
#  "2"= (RedmeatCC_2=="TRUE" )  ,"3"= (RedmeatCC_3=="TRUE" )  ,
#  .method = "first" ,
#  .default = 0  ))
# RedmeatCC<-RedmeatCC%>%dplyr::select(Metabolite,rmCCf)
# SC_PCsum_Redmeat_FV<-dplyr::left_join(SC_PCsum_Redmeat_FV,RedmeatCC,by="Metabolite")
# rm(RedmeatCC)

SC_PCsum_Redmeat_FV <- dplyr::left_join(SC_PCsum_Redmeat_FV, LABEL, by = "Metabolite")
SC_PCsum_Redmeat_FV$Met_label <- as.character(SC_PCsum_Redmeat_FV$Met_label)
#

rmPC_new_NEM_round2 <- bind_rows(filter(SC_PCsum_Redmeat_CC2, Assoc == 0), filter(SC_PCsum_Redmeat_CC3, Assoc == 0))
rmPC_new_NEM_temp1 <- SC_PCsum_Redmeat %>% filter(Metabolite %in% rmPC_new_NEM_round2$Metabolite) %>% dplyr::select(Metabolite, Assoc, DE) %>% dplyr::mutate(Metabolite = as.character(Metabolite))
rmPC_new_NEM_temp2 <- rmPC_new_NEM_round2 %>% dplyr::select(-Assoc, -DE) %>% dplyr::mutate(Metabolite = as.character(Metabolite))
rmPC_new_NEM_round2 <- left_join(rmPC_new_NEM_temp2, rmPC_new_NEM_temp1, by = "Metabolite") %>% mutate(round = as.character(3), Assoc_FV = as.character(0), DE_FV = as.character(0))

NEM_temp2_names <- filter(AMB1_PC_Redmeat, Metabolite %in% rmPC_new_NEM_round2$Metabolite)
NEM_temp2_names$Metabolite <- as.character(rmPC_new_NEM_round2$Metabolite)
rmPC_AMB_round2 <- setdiff(AMB1_PC_Redmeat, AMB1_PC_Redmeat %>% dplyr::filter(Outcome %in% NEM_temp2_names$Outcome))

SC_PCsum_Redmeat_FV <- dplyr::bind_rows(
  DE1_PC_Redmeat %>% dplyr::mutate(Assoc = as.character(Assoc), DE = as.character(DE), round = as.numeric(1), Assoc_FV = as.character(Assoc), DE_FV = as.character(DE)),
  rmPC_AMB_round2 %>% dplyr::mutate(Assoc = as.character(Assoc), DE = as.character(DE), round = as.numeric(1), Assoc_FV = as.character(Assoc), DE_FV = as.character(DE)),
  rmPC_new_NEM_round2 %>% dplyr::mutate(Assoc = as.character(Assoc), DE = as.character(DE), round = as.numeric(round), Assoc_FV = as.character(Assoc_FV), DE_FV = as.character(DE_FV)),
  filter(SC_PCsum_Redmeat %>% dplyr::mutate(Assoc = as.character(Assoc), DE = as.character(DE), round = as.numeric(1), Assoc_FV = as.character(Assoc), DE_FV = as.character(DE)), Assoc == 0)
)
# SC_PCsum_Redmeat_FV<-dplyr::left_join(SC_PCsum_Redmeat_FV,RedmeatCC,by="Metabolite")
# rm(RedmeatCC)

SC_PCsum_Redmeat_FV <- dplyr::left_join(SC_PCsum_Redmeat_FV %>% dplyr::mutate(Metabolite = as.character(Metabolite)), LABEL, by = "Metabolite")
SC_PCsum_Redmeat_FV$Met_label <- as.character(SC_PCsum_Redmeat_FV$Met_label)
SC_PCsum_Redmeat_FV <- SC_PCsum_Redmeat_FV %>% dplyr::rowwise() %>% dplyr::mutate(Est_range_FV = paste0(formatC(bestGuess, format = "f", digits = 2), " (", formatC(lowEst, format = "f", digits = 2), ", ", formatC(highEst, format = "f", digits = 2), ")"))
RedmeatCC <- CC.Redmeat.PC %>% dplyr::select(Metabolite, contains("CC"))
SC_PCsum_Redmeat_FV <- dplyr::left_join(SC_PCsum_Redmeat_FV, RedmeatCC, by = "Metabolite")

SC_PCsum_Redmeat_FV <- SC_PCsum_Redmeat_FV %>% dplyr::select(-Met_label.x, -Met_label.y)
SC_PCsum_Redmeat_FV <- dplyr::left_join(SC_PCsum_Redmeat_FV, LABEL, by = "Metabolite")

# saveRDS(SC_PC_Redmeat_FV,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/aaPC/FV/Results/Redmeat/aePC_Redmeat_FV.rds")
# saveRDS(SC_PCsum_Redmeat_FV,"L:/!MEP/Projekte/EPIC-Potsdam/Diabetes/Metabolomic traits/Simulation/PC/aaPC/FV/Results/Redmeat/aePCsum_Redmeat_FV.rds")
