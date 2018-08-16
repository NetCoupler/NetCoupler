
#############################################################
#1.1 DEFINE DAG SIMULATION

##############################################################
#  1.1     Simulation: Generate random DAG
#          and sample DAG-generated DATA

library(RBGL)
library(dplyr)
library(glmulti)
library(pcalg)
library(gRbase)
library(rlist)
library(GeneNet)
# source("https://bioconductor.org/biocLite.R")
# biocLite("RBGL")

################################################################################
# Modification of the randomDAG function from the pcalg package                #
# Function to fully specify a data-generating model (DAG) with one exposure    #
#  and nVert=number of nodes in the network (exposure excluded)                #
# exp.effects specifies which nodes are affected by the exposure               #
# exp.weights specifies the strenght of exposure effect for each affected node #
################################################################################

exp.randomDAG <-
    function(nVert,
             prob,
             lB = 0.1,
             uB = 1,
             V = as.character(1:nVert),
             exp.effects,
             exp.weights) {
        stopifnot(
            #sanity checks
            nVert >= 2,
            is.numeric(prob),
            length(prob) == 1,
            0 <= prob,
            prob <= 1,
            is.numeric(lB),
            is.numeric(uB),
            lB <= uB
        )

        edL <- vector("list", nVert + 1)
        names(edL) <- c(V, "EXP1")
        nmbEdges <- 0L

        for (i in seq_len(nVert)) {
            listSize <- rbinom(1, nVert - i, prob)
            nmbEdges <- nmbEdges + listSize
            edgeList <- sample(seq(i + 1, nVert), size = listSize)
            weightList <- runif(length(edgeList), min = lB, max = uB)
            edL[[i]] <- list(edges = edgeList, weights = weightList)
        }

        listSize <- rbinom(1, 1, prob)
        if (listSize > 0) {
            nmbEdges <- nmbEdges + 1
            edgeList <- nVert
            weightList <- runif(1, min = lB, max = uB)
        } else {
            edgeList <- integer(0)
            weightList <- numeric(0)
        }

        edL[[nVert - 1]] <- list(edges = edgeList, weights = weightList)
        if (nmbEdges > 0) {
            edL[[nVert]] <- list(edges = integer(0), weights = numeric(0))
            edL[["EXP1"]] <-
                list(edges = exp.effects, weights = exp.weights)
            new(
                "graphNEL",
                nodes = c(V, "EXP1"),
                edgeL = edL,
                edgemode = "directed"
            )
        } else {
            new("graphNEL",
                nodes =  c(V, "EXP1"),
                edgemode = "directed")
        }
    }





##############################################################
# 2.2         ESTIMATE CPDAG (skeleton) ON DAG-DATA

XDAG <- c()
Mat <- matrix()
DAGsample <- matrix()
DAGsample1 <- matrix()
#Generate a large sample (<< study sample) as source population
n <- 20000

#Define network size (number of variables in the data-generating model)
nVert <- 25
L <- LETTERS[1:nVert]
set.seed(1221)

#Generate a random DAG (data-generating model for network variables)

XDAG <-
    exp.randomDAG(
        nVert = 25,
        prob = 0.16,
        lB = 0.5,
        uB = 0.8,
        V = L,
        exp.effects = c(2, 6, 22),
        exp.weights = c(0.05 , 0.10, 0.05)
    )
plot(XDAG)
Mat <- graphNEL2M(XDAG, result = "matrix")
Mat <-
    Mat[order(match(row.names(Mat), tsort(XDAG))), order(match(colnames(Mat), tsort(XDAG)))]
XDAG <- M2graphNEL(Mat)

#Generate random observations; data-generating process defined by the DAG
DAGsample <- rmvDAG(n, XDAG, errDist = "normal")
DAGsample1 <-
    as.matrix(dplyr::select(data.frame(DAGsample), -EXP1)) #remove exposure variable, matrix  of observations of networkvariables
DAGsample1 <- scale(DAGsample1)
EXP1 <-
    as.matrix(dplyr::select(data.frame(DAGsample), EXP1))        #vector of observations of exposure-levels
nV <- colnames(DAGsample1)



####Estimate an undirected Network: sceleton of a DAG####
pcDAG_skel <- skeleton(
    suffStat = list(C = cor(DAGsample1), n = n),
    indepTest = gaussCItest,
    method = "stable",
    ## indep.test: partial correlations
    alpha = 0.05,
    labels = nV,
    fixedGaps = NULL,
    fixedEdges = NULL,
    verbose = FALSE
)

plot(XDAG)
plot(pcDAG_skel)



#############################################################
#1.3 DEFINE SURVIVAL TIME SIMULATION
#     (based on Gompertz-Distribution)
library(survival)
survsim.cw <-
    function(object,
             IV1 = X1,
             IV2 = X2,
             IV3 = X3,
             beta1,
             beta2,
             beta3) {
        n <- 20000
        U = round(runif(1, 0, 1), 2)  	#Parameter Gompertz Verteilung ->Zahl zwischen 0 und 1
        alpha = 0.2138				      	#Parameter Gompertz Verteilung
        lambdaEvent =  0.7 * 10 ** (-11)	#Parameter Gompertz Verteilung

        #X0<-EXP1
        X1 <- object[, IV1]
        X2 <- object[, IV2]
        X3 <- object[, IV3]
        RV1 <- rnorm(n, mean = 10, sd = 1)
        RV2 <- rnorm(n, mean = 10, sd = 1)
        #print(X1[1:10])
        #print(X2[1:10])
        #print(X3[1:10])

        expBetaX = exp(beta1 * X1 + beta2 * X2 + beta3 * X3 + RV1 + RV2)

        fup_time = (1 / alpha) * log(1 - ((alpha * log(U)) / (lambdaEvent * expBetaX)))


        #find the first 10% event times;
        sort_time = sort(fup_time)
        five <- n / 20
        #p10 = (sort_time[ten]+sort_time[ten+1])/2
        third <- n / 3
        p33 = (sort_time[five] + sort_time[five + 1]) / 2

        # calculate censoring variable
        cens = matrix(0, n, )
        for (i in 1:n) {
            if (fup_time[i] < p33) {
                cens[i] <- 1
            } else {
                fup_time[i] <- p33
            }
        }
        cens <- as.numeric(cens)
        #print(mean(fup_time))
        print(max(fup_time))
        #print(min(fup_time))
        print(sum(cens))
        simSurv <- bind_cols(data.frame(fup_time, cens))
    }

ST <- c()
DAG_SURV_EXP <- data.frame()
SURV <- c()
ST <-
    data.frame(
        survsim.cw(
            object = DAGsample1,
            IV1 = "D" ,
            IV2 = "J",
            IV3 = "G",
            beta1 = 0.5 ,
            beta2 = 0.5 ,
            beta3 = -0.5
        )
    )
DAG_SURV_EXP <- bind_cols(data.frame(DAGsample1), data.frame(EXP1), ST)
DAG_SURV_EXP <-
    DAG_SURV_EXP %>% mutate(
        EXP2 = rnorm(20000, mean = 0, sd = 1) ,
        EXP3 = rnorm(20000, mean = 0, sd = 1),
        ID = c(1:20000)
    )
DAG_SURV_EXP_SC <- DAG_SURV_EXP %>% sample_n(2000, replace = FALSE)
DAG_SURV_EXP_Case <-
    DAG_SURV_EXP %>% dplyr::filter(fup_time < max(DAG_SURV_EXP$fup_time))
DAG_SURV_EXP_extCase <-
    setdiff(DAG_SURV_EXP_Case, DAG_SURV_EXP_SC) %>% mutate(sc = 0,
                                                           start = (fup_time - 0.0005),
                                                           stop_t = fup_time)
DAG_SURV_EXP_SC <-
    DAG_SURV_EXP_SC %>% mutate(sc = 1, start = 0, stop_t = fup_time)
DAG_SURV_EXP_CC <- bind_rows(DAG_SURV_EXP_SC, DAG_SURV_EXP_extCase)
rm(DAG_SURV_EXP_SC)
rm(DAG_SURV_EXP_Case)
EXP <- c("EXP1", "EXP2", "EXP3")
SURV <-
    Surv(time = DAG_SURV_EXP[["fup_time"]], event = DAG_SURV_EXP[["cens"]])
fit <- coxph(Surv(fup_time, cens) ~ D + J + G, data = DAG_SURV_EXP)
summary(fit)

EXP_data <- as.matrix(dplyr::select(data.frame(DAG_SURV_EXP), EXP3))
EXP_data <-
    as.matrix(dplyr::select(data.frame(DAG_SURV_EXP), one_of(c(
        "EXP1", "EXP2", "EXP3"
    ))))
