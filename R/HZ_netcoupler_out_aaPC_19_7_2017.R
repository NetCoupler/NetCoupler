#####################################################HZ MGM Netcoupler out Clemens Wittenbecher#############################################################
################################################################aaPCs######################################################################################

rm(list=ls())

















#########################################################Rename feature names in order to avoid clash with glm#######################################

rename.met<-function(dat){
  
  #dat: samples x metabolites data matrix
  
  Ll<-paste("NM",c(1:dim(dat)[2]),sep="")                                #generate shorter metabolite names
  
  names_mapping<-cbind(colnames(dat),Ll)                                 #mapping of old and new metabolite names
  colnames(names_mapping)<-c("Metabolite","Outcome")               
  
  data_renamed<-dat
  colnames(data_renamed)<-Ll                                                #is character!
  
  return(list(data_renamed=data_renamed,names_mapping=names_mapping))                                           
  
}



#########################################################Obtain partial correlation matrix, DAG skeleton, DAG and adjacency matrix for DAG skeleton#########################


est.pcor.skel.DAG.adj<-function(dat){
  
  #dat: samples x metabolites data matrix
  
  
  
  
  #check if input data is gaussian
  
  pCor_mat<-ppcor::pcor(dat)$estimate                             #estimate Pearson's partial correlation coefficients
  colnames(pCor_mat)<-colnames(dat)
  rownames(pCor_mat)<-colnames(dat)
  
  n_samples<-nrow(dat)                                            #number of samples
  V_met<-colnames(dat)                                            #labels equal to node names i.e. metabolites
  
  skel_est<-skeleton(suffStat = list(C=cor(dat),n = n_samples),     #estimate order-independent "PC-stable" skeleton of DAG using PC-algorithm
                     indepTest = gaussCItest,                        #test conditional independence of Gaussians via Fisher's Z
                     labels = V_met, method = "stable", alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE)
  
  DAG_est<-pc(suffStat = list(C=cor(dat),n = n_samples),           #estimate equivalence class of DAG using PC-algorithm
              indepTest = gaussCItest, labels = V_met, skel.method = "stable",          #order-independent skeleton
              alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, maj.rule = FALSE, solve.confl = FALSE)
  
  adj_matrix<-get.adjacency(graphNEL2igraph(skel_est@graph))          #return adjacency matrix of DAG skeleton
  
  return(list(pCor_mat = pCor_mat, skel_est = skel_est, DAG_est = DAG_est, adj_matrix = adj_matrix))
  
}


###############################################Net Coupler Out function: metabolome -> time-to-type-2-diabetes-incident############################################################

#return: estimate of direct effects of each metabolite on survival time, models are adjusted for all possible combinations of direct neighboring metabolites and all covariates

net.coupler.out<-function(graph_skel,dat,dat_compl,exp_dat,DE,survival_obj,always_set){
  
  #graph_skel: estimated DAG skeleton of samples x metabolites data matrix
  #dat: renamed samples x metabolites data matrix
  #exp_dat: exposure/phenotype data
  #DE: indicator if direct effects were already identified
  #survival_obj: "survival" object
  #always_set: fixed set of covariates always included in model
  
  cat("*****************************************************************************************************\n")
  cat("This algorithm estimates direct effect of a predefined exposure (network-variable) on time-to-event  \n")
  cat("for all causal models that agree with the input-network: Cox prop. hazards regression models are     \n")
  cat("used to estimate the efect of all network-variables on survival time adjusted for all possible       \n")
  cat("combinations of direct neighbors (adjacency set) -> Output is a multiset of possible causal effects \n")
  cat("*****************************************************************************************************")
  
  model_details_all<-list(NULL)                                   #prepare empty list to store output
  node_names<-colnames(dat)                                      #create vector "nodes" with node names as strings
  
  for(i in seq(along=node_names)){                                #create empty list with slots for each network-variable, i.e. metabolite
    
    model_details<-list(NULL)
    model_details_all[[i]]<-model_details
    
  }
  names(model_details_all)<-node_names  
  
  for(i in node_names){                                          #net.coupler.out loop over all metabolite nodes
    
    ##########################################prepare separate datasets for each metabolite:######################################################## 
    
    met_outcome<-i
    
    #check if met_outcome is character string:
    if(is.character(met_outcome)==FALSE)stop("'met_outcome' is not a character string as required")
    
    #select data on met_outcome within samples x metabolites data matrix, store in "met_out_data":
    met_out_data<-dat[,met_outcome]
    
    #check if met_out_data is numeric:
    if(is.numeric(met_out_data)==FALSE)stop("'met_out_data' is not a numeric vector as required")
    
    #create vector with integers indicating adjacent variables, i.e. metabolites in skeleton:
    edge_list<-slot(graph_skel@graph,"edgeL")                     #extract edge-list from skeleton
    adjset<-c(edge_list[[met_outcome]])                           #extract adjacency set for selected node/metabolite
    adjset<-c(adjset[[1]])                                        #extract indices of adjacency set
    
    #select data on adjacency set, store in "adjset_data":
    adjset_data<-subset(dat_compl,select=c(adjset))
    
    #check if adjset_data is numeric:
    if(is.numeric(adjset_data[1,1])==FALSE)stop("'adjset_data' is not a matrix of numeric variables as required")
    
    #create vector with names of adjacency set as strings:
    adjset_names<-colnames(adjset_data)
    
    print(adjset_names)
    #for loop with already identified direct effects:
    if(is.null(DE)==FALSE){adjset_names<-setdiff(adjset_names,as.vector(DE))
    adjset_data<-adjset_data[,c(adjset_names)]}
    print(adjset_names)
    
    #check if adjset_names is vector of character strings:
    if(is.character(adjset_names)==FALSE)stop("'adjset_names' is not a vector of character strings as required")
    
    #combine data on exposure (including already identified direct effect metabolites), metabolite outcome, and adjacency set, store as dataframe ("modeldata"), input data for glmulti:
    modeldata<-data.frame(cbind(met_out_data,exp_dat,adjset_data))
    colnames(modeldata)[1]<-met_outcome
    
    ###########################################################estimate multimodel coefficients#####################################################
    
    #modify fitting function of coxph so that it always includes one predefined set of variables (here: exposure + covariates), while subsetting adjacent metabolite set:
    coxph.redefined<-function(formula,data,always="",...){
      coxph(as.formula(paste(deparse(formula),always)),data=data,...)
    }
    
    #fit all possible causal models using glmulti: iterate over all possible combinations of adjacent metabolites, fixed exposure + covariates set
    
    glmulti_obj<-glmulti(y=deparse(substitute(survival_obj)),                                      #response variable to be predicted: SURV-object
                         xr=c(adjset_names[1:length(adjset_names)]),            #predictor variables to be tested in iteration
                         data=modeldata,                                        #dataframe containing data
                         level=1,                                               #only main effects (terms of order 1) are used to build the candidate set
                         confsetsize = 2^(length(adjset_names)),                #number of models to be looked for, dependent on size of adjacency set
                         fitfunction=coxph.redefined,                           #fitting function to be used: here: coxph-function with fixed predictor set "always"
                         always=paste0("+",met_outcome,"+",always_set),         #defines "always"-set: complete set of exposures and covariates and current metabolite
                         includeobjects=TRUE,                                   #fitted models are included as objects
                         plotty=FALSE,                                          #no progress of the IC profile is plotted
                         report=FALSE                                           #report about progress at run time is given
    )
    
    #output summaries:
    
    glmulti_obj_objects<-glmulti_obj@objects                                    #collect glmulti objects, i.e. all fitted models (?)
    
    nbmds<-c(1:glmulti_obj@nbmods)                                              #number of fitted models of glmulti-function
    
    model_details<-list(NULL)
    
    for(j in seq(along=nbmds)){                                                 #generate list for model details
      model_details[[j]]<-list(NULL)
    }
    names(model_details)<-paste("Model",nbmds,sep="_")
    
    for(j in seq(along=nbmds)){                                                 #collect all model details for different adjacency subsets for specific metabolite
      model_summary<-glmulti_obj_objects[[j]]
      model_details[[j]]<-list(Model=paste("Model",j,"of",length(nbmds)),
                               Model_summary=model_summary)
    }
    
    model_details_all[[met_outcome]]<-list(Model_summaries=model_details,      #collect all model details for different adjacency subsets for specific metabolite in complete model list for all metabolites
                                           Number_of_Models=length(nbmds),
                                           Adj_set=paste(adjset_names,collapse = ", "),
                                           Outcome_metabolite=met_outcome)
    
  }  
  
  return(list(model_details_all=model_details_all,outcome=list(survival_obj=survival_obj,class=class(survival_obj))))
  
}



####################################################Extract exposure coefficients per outcome, i.e. metabolite#######################################################

getExp.coef.permetabolite<-function(object,metabolite){
  
  #object: output of net.coupler.out
  #metabolite: specific metabolite to evaluate
  
  #create vector containing integers from 1 to number of metabolite-specific models:
  nbm<-c(1:length(object$model_details_all[[metabolite]]$Model_summaries))
  Nbmds<-max(nbm)
  
  #define empty dataframes for output:
  mm_coef_temp1<-structure(list(Model=as.character(),
                                Nbmds=as.numeric(),
                                Outcome=as.character(),
                                Covariables=as.character(),
                                #NCov=as.numeric(),
                                HR=as.numeric(),
                                LCL=as.numeric(),
                                UCL=as.numeric(),
                                Beta=as.numeric(),
                                rSE=as.numeric(),
                                P=as.numeric()),
                           class="data.frame")
  
  mm_coef_temp2<-structure(list(Model=as.character(),
                                Nbmds=as.numeric(),
                                Outcome=as.character(),
                                Covariables=as.character(),
                                NCov=as.numeric(),
                                HR=as.numeric(),
                                LCL=as.numeric(),
                                UCL=as.numeric(),
                                Beta=as.numeric(),
                                rSE=as.numeric(),
                                P=as.numeric()),
                           class="data.frame")
  
  #loop along number of metabolite-specific models:
  for(i in seq(along=nbm)){
    #get metabolite-effect estimates (hazard-ratio, lower and upper confidence limit, beta, robust SE, and p-value) from single model and write into dataframe
    SUM<-summary(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model_summary)
    Cov<-rownames(SUM$coefficients)
    #Cov<-dplyr::filter(data.frame(row.names(SUM$coefficients)),row.names(SUM$coefficients)!=exposure)                     #record specific covariates of this model
    
    mm_coef_temp2<-data.frame(Model=as.character(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model),
                              Nbmds=Nbmds,
                              Outcome=metabolite,
                              Covariables=as.character(paste(Cov,collapse=", ")),
                              NCov=max(unlist(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model_summary$assign))-1,
                              HR=as.numeric(SUM$conf.int[[metabolite,1]]),
                              LCL=as.numeric(SUM$conf.int[[metabolite,3]]),
                              UCL=as.numeric(SUM$conf.int[[metabolite,4]]),
                              Beta=as.numeric(SUM$coefficients[[metabolite,1]]),
                              rSE=as.numeric(SUM$coefficients[[metabolite,4]]),
                              P=as.numeric(SUM$coefficients[[metabolite,6]]))
    
    #bind information to an outcome-specific dataframe:
    mm_coef_temp1<-bind_rows(mm_coef_temp1,mm_coef_temp2)
    
  }
  
  return(mm_coef_temp1)
  
}





####################################################Extract exposure coefficients for metabolite(s) on time-to-incidence###############################################

getExp.coef.out<-function(object,metabolite){
  
  #object: output of net.coupler
  #metabolite: specific metabolite to evaluate
  
  cat("*************************************************************************************************** \n")
  cat("This function produces a table of effect estimates of all (some) network-variables on an outcome    \n")
  cat("(time-to-event) for all possibles causal models based on conditional independence criteria encoded  \n")
  cat("in the input-network => MULTISET OF POSSIBLE EFFECTS PER EXPOSURE-OUTCOME PAIR; network-variables   \n")
  cat("of interest are selected by indicating the variable-names as character vector                       \n")
  cat("***************************************************************************************************")
  
  
  #define empty dataframes for output:
  
  mm_coef<-structure(list(Model=as.character(),
                          Nbmds=as.numeric(),
                          Outcome=as.character(),
                          Covariables=as.character(),
                          NCov=as.numeric(),
                          HR=as.numeric(),
                          LCL=as.numeric(),
                          UCL=as.numeric(),
                          Beta=as.numeric(),
                          rSE=as.numeric(),
                          P=as.numeric()),
                     class="data.frame")
  
  mm_coef_temp<-structure(list(Model=as.character(),
                               Nbmds=as.numeric(),
                               Outcome=as.character(),
                               Covariables=as.character(),
                               NCov=as.numeric(),
                               HR=as.numeric(),
                               LCL=as.numeric(),
                               UCL=as.numeric(),
                               Beta=as.numeric(),
                               rSE=as.numeric(),
                               P=as.numeric()),
                          class="data.frame")
  
  #for each metabolite: get metabolite coefficients from all possible models and write to table:
  
  for(j in metabolite){
    mm_coef_temp<-data.frame(getExp.coef.permetabolite(object=object,metabolite=j))
    mm_coef<-bind_rows(mm_coef,mm_coef_temp)
  }
  
  return(mm_coef)
  
}


#####################Summary statistics and multiple testing adjustment for net.coupler.out with survival object############################################################################


mult.stat.surv<-function(sum_netout,adjust_method,rule_1 = 12,rule_1_cut = 0.1,rule_2 = 5,rule_2_cut = 0.05,rule_3 = 15,rule_4 = 14, ass_rule1=16, round_number){
  
  #sum_netout: output of getExp.coef.out after adding original metabolite names
  #adjust_method: select adjustment method(s), e.g. "fdr" or c("fdr","bonferroni"), see p.adjust.methods for available options
  #rule= "set DE=1 or =2 if ..."
  #rule_1: determines first argument of rule (e.g. adjusted p-value of full model < rule_1_cut)
  #rule_2: determines second argument of rule (e.g. upperP < rule_2_cut)
  #rule_3: determines first variable of third argument of rule (e.g. highEst)
  #rule_4: determines second variable of third argument of rule (e.g. lowEst)
  #possible arguments for rule_...: adjusted p-values of full model=12=rule_1
  #                                 largest p-value across all adjacency models=upperP=5=rule_2
  #                                 highest estimated beta across all adjacency models=highEst=15=rule_3
  #                                 lowest estimated beta across all adjacency models=lowEst=14=rule_4
  #default rule: if(adjusted p< 0.1 & upperP<0.05 & (sign(highEst)==sign(lowEst)) & lowEst>0){DE=1}
  #              if(adjusted p< 0.1 & upperP<0.05 & (sign(highEst)==sign(lowEst)) & lowEst<0){DE=2}else{DE=0}
  #ass_rule1: determines first variable of second argument of association rule (e.g. bestGuess=16)
  #default association rule: if(adjusted p<0.1 & bestGuess>0){Assoc=1}
  #                          if(adjusted p<0.1 & bestGuess<0){Assoc=2}else{Assoc=0}
  #round_number: actual round number in ambiguous metabolites loop
  
  #return: for each metabolite: maximum number of models, mean hazard ratio, maximum HR, minimum HR, maximum p-value, minimum p-value, as well as number of covariates, HR and p-value for model with full adjustment:
  sum_netout_sum<-dplyr::left_join(sum_netout%>%dplyr::group_by(Outcome)%>%dplyr::summarise(avgHR=mean(HR),minHR=min(HR),maxHR=max(HR),upperP=max(P),lowerP=min(P),Nbmds=max(Nbmds)),
                                   sum_netout%>%dplyr::group_by(Outcome)%>%dplyr::filter(nchar(Covariables)==max(nchar(Covariables)))%>%dplyr::select(NCov=NCov,Outcome=Outcome,Metabolite=Metabolite,ind_HR=HR,ind_P=P),
                                   by="Outcome")
  # print(head(sum_netout_sum))
  # 
  # sum_netout_sum<-sum_netout_sum%>%group_by(Outcome)%>%filter(row_number(avgHR)==1)
  # 
  # print(head(sum_netout_sum))
  # 
  #multiple testing correction:
  p_adjust_M<-p.adjust.methods[p.adjust.methods %in% adjust_method]                            #select multiple testing correction method(s)
  p_adj<-sapply(p_adjust_M,function(meth){p.adjust(sum_netout_sum$ind_P,meth)})            #calculate adjusted p-values for ind_P
  
  sum_netout_sum_adj<-bind_cols(sum_netout_sum,data.frame(p_adj))
  colnames(sum_netout_sum_adj)[12]<-adjust_method
  
  #add beta-coefficients:
  sum_netout_sum_adj<-sum_netout_sum_adj%>%dplyr::mutate(avgEst=log(avgHR),lowEst=log(minHR),highEst=log(maxHR),bestGuess=log(ind_HR))
  
  #is there an effect of specified metabolite on time-to-event?: determine effect-indicator (DE!=0: direct effect of metabolite on time-to-event, DE=0: ambiguous if adjusted ind_p<0.1):
  sum_netout_sum_adj<-mutate(sum_netout_sum_adj,
                             DE=derivedFactor(
                               "1"=(sum_netout_sum_adj[,rule_1]<rule_1_cut & sum_netout_sum_adj[,rule_2]<rule_2_cut & (sign(sum_netout_sum_adj[,rule_3])==sign(sum_netout_sum_adj[,rule_4])) & sum_netout_sum_adj[,rule_4]>0),
                               "2"=(sum_netout_sum_adj[,rule_1]<rule_1_cut & sum_netout_sum_adj[,rule_2]<rule_2_cut & (sign(sum_netout_sum_adj[,rule_3])==sign(sum_netout_sum_adj[,rule_4])) & sum_netout_sum_adj[,rule_4]<0),
                               .method="first",.default=0)
  )
  
  
  #is there any effect of specified metabolite on time-to-event?: determine association-indicator (Assoc!=0: there is an effect (no differenciation between direct or indirect), Assoc=0: there is no effect)
  sum_netout_sum_adj<-mutate(sum_netout_sum_adj,
                             Assoc=derivedFactor(
                               "1"=(sum_netout_sum_adj[,rule_1]<rule_1_cut & sum_netout_sum_adj[,ass_rule1]>0),
                               "2"=(sum_netout_sum_adj[,rule_1]<rule_1_cut & sum_netout_sum_adj[,ass_rule1]<0),
                               .method="first",.default=0)
  )
  
  #collect ambiguous effects:
  amb<-sum_netout_sum_adj%>%filter(Assoc!=0 & DE==0)
  
  #collect direct effects:
  direct<-sum_netout_sum_adj%>%filter(Assoc!=0 & DE!=0)
  
  #collect summary statistics including round information
  sum_netout_sum_adj_FV<-sum_netout_sum_adj%>%dplyr::mutate(round=as.numeric(round_number),Assoc_FV=as.character(Assoc),DE_FV=as.character(DE))
  
  return(list(sum_netout_sum_adj=sum_netout_sum_adj,amb=amb,direct=direct,sum_netout_sum_adj_FV=sum_netout_sum_adj_FV))
  
}







####################################################Data import and analysis###########################################################################


#load metabolite data (e.g. aaPCs) for complete cohort:
gpl_data <- readRDS("C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/R_data/STR_GPL.rds")

met_data<-dplyr::select(gpl_data,contains("_aa_"))

dim(met_data)                             #case-cohort: subcohort + all incident diabetes cases in complete cohort
#[1] 2731   34

#metabolite data (e.g. aaPCs) for subcohort:
gpl_data_SC <- readRDS("C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/R_data/STR_GPL_SC.rds")

met_data_SC<-dplyr::select(gpl_data_SC,contains("_aa_"))

dim(met_data_SC)       #metabolite data with rows=patients and columns=metabolites
#[1] 2092   34

is.numeric(as.matrix(met_data_SC))
#[1] TRUE

#load type-2 diabetes incident times:
load("C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/R_data/T2D_data")

dim(T2D_data)
#[1] 2731   26                                    #phenotype data including time information for 2731 patients, including diet information, lifestyle, etc.

colnames(T2D_data)
# [1] "subcohort"  "SEX"        "fall"       "Med_Hypert" "Med_HLipid" "sport"      "bike"       "alk_1"      "alk_2"      "alk_3"      "alk_4"      "alk_5"      "alk_6"     
# [14] "sta_time"   "sto_time"   "WGBperMJ"   "TMperMJ"    "age_years"  "smk1"       "smk2"       "smk3"       "educ1"      "educ2"      "educ3"      "CofCup"     "ID"

#load phenotype data for complete cohort:
Exp_data<-readRDS("C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/R_data/EXP_DATA.rds")

dim(Exp_data)
#[1] 2731   56                                   #phenotype data for 2731 patients, including diet information, lifestyle, etc.


met_data_SC_rename<-rename.met(dat=met_data_SC)$data_renamed                               #rename metabolites with short names
met_mapping_SC<-rename.met(dat=met_data_SC)$names_mapping                                     #mapping information between old and new metabolite names

met_data_rename<-rename.met(dat=met_data)$data_renamed                               #rename metabolites with short names
met_mapping<-rename.met(dat=met_data)$names_mapping                                     #mapping information between old and new metabolite names

#skeleton estimation for metabolite matrix subcohort:
met_SC_estimates_norename<-est.pcor.skel.DAG.adj(dat=met_data_SC)                       
met_SC_skel_norename<-met_SC_estimates_norename$skel_est                                                   #estimate DAG skeleton for non-renamed matrix
met_SC_DAG_norename<-met_SC_estimates_norename$DAG_est                                                     #estimate DAG for non-renamed matrix
met_SC_adj_mat_norename<-met_SC_estimates_norename$adj_matrix                                              #estimate adjacency matrix for non-renamed matrix

met_SC_estimates_rename<-est.pcor.skel.DAG.adj(dat=met_data_SC_rename)
met_SC_skel_rename<-met_SC_estimates_rename$skel_est                                                       #estimate DAG skeleton for renamed matrix
met_SC_DAG_rename<-met_SC_estimates_rename$DAG_est                                                         #estimate DAG for renamed matrix
met_SC_adj_mat_rename<-met_SC_estimates_rename$adj_matrix                                                  #estimate adjacency matrix for renamed matrix

#skeleton estimation for complete cohort metabolite matrix:
met_estimates_norename<-est.pcor.skel.DAG.adj(dat=met_data)                       
met_skel_norename<-met_estimates_norename$skel_est                                                   #estimate DAG skeleton for non-renamed matrix
met_DAG_norename<-met_estimates_norename$DAG_est                                                     #estimate DAG for non-renamed matrix
met_adj_mat_norename<-met_estimates_norename$adj_matrix                                              #estimate adjacency matrix for non-renamed matrix

met_estimates_rename<-est.pcor.skel.DAG.adj(dat=met_data_rename)
met_skel_rename<-met_estimates_rename$skel_est                                                       #estimate DAG skeleton for renamed matrix
met_DAG_rename<-met_estimates_rename$DAG_est                                                         #estimate DAG for renamed matrix
met_adj_mat_rename<-met_estimates_rename$adj_matrix                                                  #estimate adjacency matrix for renamed matrix

#create "survival" object:
t2d_surv<-Surv(T2D_data$sta_time,T2D_data$sto_time,T2D_data$fall)

#define "always"-set:
always_set<-paste0(colnames(Exp_data)[-(which(colnames(Exp_data) %in% c("SEX","subcohort","ID","age","age_years")))],collapse = " + ")
always_set<-paste0(always_set," + cluster(ID) + strata(age_years)")                                   #cluster(ID) specific for case-cohort design, strata(age_years) stratification of baseline risk according to age in years

#initial net.coupler.out run:
net_coupler_out_1<-net.coupler.out(graph_skel=met_SC_skel_rename,dat=met_data_rename,dat_compl = met_data_rename,exp_dat = Exp_data,DE=NULL,survival_obj = t2d_surv,always_set = always_set)

#return results:
sum_netout<-getExp.coef.out(object=net_coupler_out_1,metabolite=colnames(met_data_rename))

#get original metabolite names back:
sum_netout<-merge(sum_netout,as.data.frame(met_mapping),by="Outcome")

#summary statistics:
sum_netout_stat<-mult.stat.surv(sum_netout=sum_netout,adjust_method = "fdr",round_number = 1)

netout_amb<-sum_netout_stat$amb                                            #summary statistics for metabolites which have ambiguous effect on time-to-diabetes-incident
netout_direct<-sum_netout_stat$direct                                      #summary statistics for metabolites which have direct effect on time-to-diabetes-incident
netout_sum<-sum_netout_stat$sum_netout_sum_adj_FV                          #summary statistics for all metabolites including round-number information

netout_amb
# A tibble: 10 x 18
# Outcome     avgHR     minHR     maxHR     upperP       lowerP Nbmds  NCov   Metabolite    ind_HR        ind_P          fdr     avgEst       lowEst      highEst  bestGuess     DE  Assoc
# <chr>     <dbl>     <dbl>     <dbl>      <dbl>        <dbl> <dbl> <dbl>       <fctr>     <dbl>        <dbl>        <dbl>      <dbl>        <dbl>        <dbl>      <dbl> <fctr> <fctr>
#   1    NM17 0.7556374 0.5291614 0.9699910 0.56107622 7.225901e-05    16    55 rPC_aa_C36_6 0.5291614 7.225901e-05 0.0012284032 -0.2801937 -0.636461791 -0.030468480 -0.6364618      0      2
# 2    NM22 0.8281669 0.4999506 1.1166445 0.45824343 2.670005e-05     4    53 rPC_aa_C38_5 0.4999506 2.670005e-05 0.0009078018 -0.1885406 -0.693245948  0.110328246 -0.6932459      0      2
# 3    NM24 1.1598871 0.9102137 1.4114873 0.58403875 2.131131e-03     8    54 rPC_aa_C40_2 1.3138229 1.722860e-02 0.0544139432  0.1483227 -0.094075843  0.344643978  0.2729411      0      1
# 4    NM25 0.8431124 0.6872579 1.0694483 0.92305197 1.598933e-05    32    56 rPC_aa_C40_3 0.7402252 4.059156e-03 0.0191980104 -0.1706550 -0.375045720  0.067142877 -0.3008008      0      2
# 5    NM28 1.3518927 1.0041347 1.7245665 0.95573685 1.941793e-06     4    53 rPC_aa_C40_6 1.5323859 2.781241e-03 0.0157603634  0.3015056  0.004126136  0.544975687  0.4268259      0      1
# 6    NM29 0.7433856 0.6474714 0.8645700 0.13146880 5.081305e-10    16    55 rPC_aa_C42_0 0.8001438 3.103473e-02 0.0879317218 -0.2965404 -0.434680606 -0.145523015 -0.2229638      0      2
# 7    NM30 0.7097091 0.5849164 0.8461914 0.09128099 2.364775e-14    16    55 rPC_aa_C42_1 0.7264600 4.517179e-03 0.0191980104 -0.3429002 -0.536286344 -0.167009664 -0.3195719      0      2
# 8    NM31 0.7949370 0.6670219 1.0083521 0.90631060 5.127145e-07    64    57 rPC_aa_C42_2 0.7998989 1.760451e-02 0.0544139432 -0.2294924 -0.404932373  0.008317406 -0.2232699      0      2
# 9     NM6 0.8465499 0.7448836 0.9554412 0.47563214 1.720361e-04     8    54 rPC_aa_C32_3 0.7734183 2.636267e-03 0.0157603634 -0.1665862 -0.294527374 -0.045582049 -0.2569353      0      2
# 10     NM7 0.9036881 0.6694385 1.1860704 0.93535740 1.428366e-03    16    55 rPC_aa_C34_1 0.6694385 1.412586e-02 0.0533643545 -0.1012710 -0.401316035  0.170645696 -0.4013160      0      2

netout_direct
# A tibble: 2 x 18
# Outcome    avgHR    minHR    maxHR      upperP       lowerP Nbmds  NCov   Metabolite   ind_HR        ind_P         fdr    avgEst    lowEst   highEst bestGuess     DE  Assoc
# <chr>    <dbl>    <dbl>    <dbl>       <dbl>        <dbl> <dbl> <dbl>       <fctr>    <dbl>        <dbl>       <dbl>     <dbl>     <dbl>     <dbl>     <dbl> <fctr> <fctr>
#   1    NM20 1.461412 1.331850 1.539315 0.001463337 2.543242e-07    16    55 rPC_aa_C38_3 1.539315 0.0005648558 0.004801274 0.3794032 0.2865691 0.4313374 0.4313374      1      1
# 2    NM33 1.305813 1.108655 1.477687 0.045878920 1.160349e-05     8    54 rPC_aa_C42_5 1.477687 0.0001603423 0.001817212 0.2668261 0.1031473 0.3904780 0.3904780      1      1


#there are direct and ambiguous effect metabolites on time-to-diabetes-incident:

#define new metabolite data matrix only including ambiguous effect metabolites:
met_data_rename_amb1_netout<-as.matrix(met_data_rename[,c(which(colnames(met_data_rename) %in% netout_amb$Outcome)),drop=FALSE])                     

DE<-netout_direct$Outcome
DE_1<-netout_direct$Outcome

#add direct effect metabolites to Exp_data:
Exp_data_DE1_netout<-cbind(met_data_rename[,c(which(colnames(met_data_rename) %in% DE_1))],Exp_data)
colnames(Exp_data_DE1_netout)<-c(DE_1,colnames(Exp_data))

#add direct effect metabolites to always set:
always_set_DE1<-paste0(paste0(DE_1,collapse=" + "),sep=" + ",always_set)

#repeat net.coupler.out with new input and new always set now including all metabolites previously identified as direct effect metabolites:
net_coupler_out_2<-net.coupler.out(graph_skel=met_SC_skel_rename,dat=met_data_rename_amb1_netout,dat_compl = met_data_rename,exp_dat = Exp_data_DE1_netout,DE=DE,survival_obj = t2d_surv,always_set = always_set_DE1)

sum_netout_2<-getExp.coef.out(object=net_coupler_out_2,metabolite=colnames(met_data_rename_amb1_netout))

#get original metabolite names back:
sum_netout_2<-merge(sum_netout_2,as.data.frame(met_mapping),by="Outcome")

#summary statistics:
sum_netout_stat_2<-mult.stat.surv(sum_netout=sum_netout_2,adjust_method = "fdr",round_number = 2)

netout_amb_2<-sum_netout_stat_2$amb                                            #summary statistics for metabolites which have ambiguous effect on time-to-diabetes-incident
netout_direct_2<-sum_netout_stat_2$direct                                      #summary statistics for metabolites which have direct effect on time-to-diabetes-incident
netout_sum_2<-sum_netout_stat_2$sum_netout_sum_adj_FV                          #summary statistics for all metabolites including round-number information

netout_amb_2
# A tibble: 6 x 18
# Outcome     avgHR     minHR     maxHR     upperP       lowerP Nbmds  NCov   Metabolite    ind_HR       ind_P        fdr      avgEst      lowEst     highEst  bestGuess     DE  Assoc
# <chr>     <dbl>     <dbl>     <dbl>      <dbl>        <dbl> <dbl> <dbl>       <fctr>     <dbl>       <dbl>      <dbl>       <dbl>       <dbl>       <dbl>      <dbl> <fctr> <fctr>
#   1    NM17 0.7069553 0.6179231 0.7999431 0.05491281 8.773698e-05    16    57 rPC_aa_C36_6 0.6711242 0.018792788 0.03643585 -0.34678783 -0.48139122 -0.22321466 -0.3988011      0      2
# 2    NM24 1.0758483 0.8366412 1.3039004 0.61792278 1.752766e-02     8    56 rPC_aa_C40_2 1.2962680 0.025505094 0.03643585  0.07310942 -0.17835995  0.26536006  0.2594894      0      1
# 3    NM25 0.7918270 0.6902090 0.9126297 0.29621698 3.282691e-05    16    57 rPC_aa_C40_3 0.7789921 0.020513652 0.03643585 -0.23341231 -0.37076085 -0.09142504 -0.2497543      0      2
# 4    NM28 1.1886150 0.9495486 1.4881860 0.73113704 5.552096e-03     4    55 rPC_aa_C40_6 1.4881860 0.005552096 0.01850699  0.17278880 -0.05176853  0.39755793  0.3975579      0      1
# 5    NM29 0.7560749 0.6637714 0.8577500 0.12462379 9.368506e-12    16    57 rPC_aa_C42_0 0.8326282 0.084087987 0.09343110 -0.27961487 -0.40981751 -0.15344264 -0.1831680      0      2
# 6    NM31 0.7747592 0.6604121 0.9053642 0.25216965 2.021956e-08    64    59 rPC_aa_C42_2 0.8299504 0.048311345 0.06038918 -0.25520297 -0.41489124 -0.09941800 -0.1863894      0      2

netout_direct_2
# A tibble: 3 x 18
# Outcome     avgHR     minHR     maxHR       upperP       lowerP Nbmds  NCov   Metabolite    ind_HR        ind_P          fdr     avgEst     lowEst    highEst  bestGuess     DE  Assoc
# <chr>     <dbl>     <dbl>     <dbl>        <dbl>        <dbl> <dbl> <dbl>       <fctr>     <dbl>        <dbl>        <dbl>      <dbl>      <dbl>      <dbl>      <dbl> <fctr> <fctr>
#   1    NM22 0.6029201 0.4688673 0.6966110 0.0002019488 6.909858e-06     4    55 rPC_aa_C38_5 0.4688673 6.909858e-06 6.909858e-05 -0.5059705 -0.7574356 -0.3615282 -0.7574356      2      2
# 2    NM30 0.6854735 0.6126198 0.7669771 0.0150326314 3.380629e-13    16    57 rPC_aa_C42_1 0.7289033 5.043250e-03 1.850699e-02 -0.3776455 -0.4900107 -0.2652984 -0.3162142      2      2
# 3     NM6 0.7999596 0.7485887 0.8481255 0.0247436411 7.576171e-05     8    56 rPC_aa_C32_3 0.8218709 2.474364e-02 3.643585e-02 -0.2231941 -0.2895656 -0.1647267 -0.1961719      2      2

#there are additional direct and ambiguous effect metabolites on time-to-diabetes-incident:

#define new metabolite data matrix only including ambiguous effect metabolites:
met_data_rename_amb2_netout<-as.matrix(met_data_rename[,c(which(colnames(met_data_rename) %in% netout_amb_2$Outcome)),drop=FALSE])                     

DE_2<-netout_direct_2$Outcome
DE_1_2<-c(DE_1,DE_2)

#add direct effect metabolites to Exp_data:
Exp_data_DE2_netout<-cbind(met_data_rename[,c(which(colnames(met_data_rename) %in% DE_2))],Exp_data_DE1_netout)
colnames(Exp_data_DE2_netout)<-c(DE_2,colnames(Exp_data_DE1_netout))

#add direct effect metabolites to always set:
always_set_DE2<-paste0(paste0(DE_2,collapse=" + "),sep=" + ",always_set_DE1)

#repeat net.coupler.out with new input and new always set now including all metabolites previously identified as direct effect metabolites:
net_coupler_out_3<-net.coupler.out(graph_skel=met_SC_skel_rename,dat=met_data_rename_amb2_netout,dat_compl = met_data_rename,exp_dat = Exp_data_DE2_netout,DE=DE_1_2,survival_obj = t2d_surv,always_set = always_set_DE2)

sum_netout_3<-getExp.coef.out(object=net_coupler_out_3,metabolite=colnames(met_data_rename_amb2_netout))

#get original metabolite names back:
sum_netout_3<-merge(sum_netout_3,as.data.frame(met_mapping),by="Outcome")

#summary statistics:
sum_netout_stat_3<-mult.stat.surv(sum_netout=sum_netout_3,adjust_method = "fdr",round_number = 3)

netout_amb_3<-sum_netout_stat_3$amb                                            #summary statistics for metabolites which have ambiguous effect on time-to-diabetes-incident
netout_direct_3<-sum_netout_stat_3$direct                                      #summary statistics for metabolites which have direct effect on time-to-diabetes-incident
netout_sum_3<-sum_netout_stat_3$sum_netout_sum_adj_FV                          #summary statistics for all metabolites including round-number information

netout_amb_3
# A tibble: 5 x 18
# Outcome     avgHR     minHR     maxHR    upperP      lowerP Nbmds  NCov   Metabolite    ind_HR       ind_P        fdr     avgEst      lowEst     highEst  bestGuess     DE  Assoc
# <chr>     <dbl>     <dbl>     <dbl>     <dbl>       <dbl> <dbl> <dbl>       <fctr>     <dbl>       <dbl>      <dbl>      <dbl>       <dbl>       <dbl>      <dbl> <fctr> <fctr>
#   1    NM17 0.8065255 0.6110597 0.9657696 0.7690416 0.005109276    16    60 rPC_aa_C36_6 0.6294707 0.009917459 0.02975238 -0.2150198 -0.49256057 -0.03483002 -0.4628759      0      2
# 2    NM24 1.1859421 1.0342496 1.3427333 0.6812973 0.008297387     8    59 rPC_aa_C40_2 1.3060686 0.018227520 0.03470137  0.1705375  0.03367615  0.29470734  0.2670216      0      1
# 3    NM25 0.8159863 0.7386868 0.9051769 0.2450431 0.005109764     8    59 rPC_aa_C40_3 0.7405413 0.005569742 0.02975238 -0.2033578 -0.30288131 -0.09962484 -0.3003738      0      2
# 4    NM29 0.8289516 0.8009908 0.8703262 0.1695697 0.036505825     8    59 rPC_aa_C42_0 0.8152448 0.056753682 0.06810442 -0.1875935 -0.22190579 -0.13888723 -0.2042669      0      2
# 5    NM31 0.8471528 0.7909966 0.9189214 0.3306723 0.014083141    32    61 rPC_aa_C42_2 0.8036937 0.023134249 0.03470137 -0.1658742 -0.23446156 -0.08455471 -0.2185371      0      2

netout_direct_3
# A tibble: 0 x 18
# ... with 18 variables: Outcome <chr>, avgHR <dbl>, minHR <dbl>, maxHR <dbl>, upperP <dbl>, lowerP <dbl>, Nbmds <dbl>, NCov <dbl>, Metabolite <fctr>, ind_HR <dbl>, ind_P <dbl>, fdr <dbl>, avgEst <dbl>,
#   lowEst <dbl>, highEst <dbl>, bestGuess <dbl>, DE <fctr>, Assoc <fctr>



#final results:

netout_sum_2_3<-merge(netout_sum_2,netout_sum_3,by="Metabolite",suffixes=c("_2","_3"))

netout_sum_2<-netout_sum_2[,c(9,1:8,10:21)]
colnames(netout_sum_2)[2:21]<-paste0(colnames(netout_sum_2)[2:21],sep="_","2")

netout_sum_2_3_final<-bind_rows(netout_sum_2[-which(netout_sum_2$Outcome_2 %in% netout_sum_2_3$Outcome_2),],netout_sum_2_3)

netout_sum_1_2_3<-merge(netout_sum,netout_sum_2_3_final,by="Metabolite")
colnames(netout_sum_1_2_3)[2:21]<-paste0(colnames(netout_sum_1_2_3)[2:21],sep="_","1")

netout_sum_1_2_3_final<-bind_rows(netout_sum[-which(netout_sum$Outcome %in% netout_sum_1_2_3$Outcome_1),],netout_sum_1_2_3)



write.xlsx(netout_sum_1_2_3_final,file="C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/HZ_netout_aaPC_20_7_2017.xls")












#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
##############################################################shorter analysis version with loop for ambiguous metabolites#######################################################################################################
#################################################################################################################################################################################################################################

rm(list=ls())

















#########################################################Rename feature names in order to avoid clash with glm#######################################

rename.met<-function(dat){
  
  #dat: samples x metabolites data matrix
  
  Ll<-paste("NM",c(1:dim(dat)[2]),sep="")                                #generate shorter metabolite names
  
  names_mapping<-cbind(colnames(dat),Ll)                                 #mapping of old and new metabolite names
  colnames(names_mapping)<-c("Metabolite","Outcome")               
  
  data_renamed<-dat
  colnames(data_renamed)<-Ll                                                #is character!
  
  return(list(data_renamed=data_renamed,names_mapping=names_mapping))                                           
  
}



#########################################################Obtain partial correlation matrix, DAG skeleton, DAG and adjacency matrix for DAG skeleton#########################


est.pcor.skel.DAG.adj<-function(dat){
  
  #dat: samples x metabolites data matrix
  
  
  
  
  #check if input data is gaussian
  
  pCor_mat<-ppcor::pcor(dat)$estimate                             #estimate Pearson's partial correlation coefficients
  colnames(pCor_mat)<-colnames(dat)
  rownames(pCor_mat)<-colnames(dat)
  
  n_samples<-nrow(dat)                                            #number of samples
  V_met<-colnames(dat)                                            #labels equal to node names i.e. metabolites
  
  skel_est<-skeleton(suffStat = list(C=cor(dat),n = n_samples),     #estimate order-independent "PC-stable" skeleton of DAG using PC-algorithm
                     indepTest = gaussCItest,                        #test conditional independence of Gaussians via Fisher's Z
                     labels = V_met, method = "stable", alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE)
  
  DAG_est<-pc(suffStat = list(C=cor(dat),n = n_samples),           #estimate equivalence class of DAG using PC-algorithm
              indepTest = gaussCItest, labels = V_met, skel.method = "stable",          #order-independent skeleton
              alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE, maj.rule = FALSE, solve.confl = FALSE)
  
  adj_matrix<-get.adjacency(graphNEL2igraph(skel_est@graph))          #return adjacency matrix of DAG skeleton
  
  return(list(pCor_mat = pCor_mat, skel_est = skel_est, DAG_est = DAG_est, adj_matrix = adj_matrix))
  
}


###############################################Net Coupler Out function: metabolome -> time-to-type-2-diabetes-incident############################################################

#return: estimate of direct effects of each metabolite on survival time, models are adjusted for all possible combinations of direct neighboring metabolites and all covariates

net.coupler.out<-function(graph_skel,dat,dat_compl,exp_dat,DE,survival_obj,always_set,name_surv_obj){
  
  #graph_skel: estimated DAG skeleton of samples x metabolites data matrix
  #dat: renamed samples x metabolites data matrix
  #exp_dat: exposure/phenotype data
  #DE: indicator if direct effects were already identified
  #survival_obj: "survival" object
  #always_set: fixed set of covariates always included in model
  
  cat("*****************************************************************************************************\n")
  cat("This algorithm estimates direct effect of a predefined exposure (network-variable) on time-to-event  \n")
  cat("for all causal models that agree with the input-network: Cox prop. hazards regression models are     \n")
  cat("used to estimate the efect of all network-variables on survival time adjusted for all possible       \n")
  cat("combinations of direct neighbors (adjacency set) -> Output is a multiset of possible causal effects \n")
  cat("*****************************************************************************************************")
  
  model_details_all<-list(NULL)                                   #prepare empty list to store output
  node_names<-colnames(dat)                                      #create vector "nodes" with node names as strings
  
  for(i in seq(along=node_names)){                                #create empty list with slots for each network-variable, i.e. metabolite
    
    model_details<-list(NULL)
    model_details_all[[i]]<-model_details
    
  }
  names(model_details_all)<-node_names  
  
  for(i in node_names){                                          #net.coupler.out loop over all metabolite nodes
    
    ##########################################prepare separate datasets for each metabolite:######################################################## 
    
    met_outcome<-i
    
    #check if met_outcome is character string:
    if(is.character(met_outcome)==FALSE)stop("'met_outcome' is not a character string as required")
    
    #select data on met_outcome within samples x metabolites data matrix, store in "met_out_data":
    met_out_data<-dat[,met_outcome]
    
    #check if met_out_data is numeric:
    if(is.numeric(met_out_data)==FALSE)stop("'met_out_data' is not a numeric vector as required")
    
    #create vector with integers indicating adjacent variables, i.e. metabolites in skeleton:
    edge_list<-slot(graph_skel@graph,"edgeL")                     #extract edge-list from skeleton
    adjset<-c(edge_list[[met_outcome]])                           #extract adjacency set for selected node/metabolite
    adjset<-c(adjset[[1]])                                        #extract indices of adjacency set
    
    #select data on adjacency set, store in "adjset_data":
    adjset_data<-subset(dat_compl,select=c(adjset))
    
    #check if adjset_data is numeric:
    if(is.numeric(adjset_data[1,1])==FALSE)stop("'adjset_data' is not a matrix of numeric variables as required")
    
    #create vector with names of adjacency set as strings:
    adjset_names<-colnames(adjset_data)
    
    #for loop with already identified direct effects:
    if(is.null(DE)==FALSE){adjset_names<-setdiff(adjset_names,as.vector(DE))
    adjset_data<-adjset_data[,c(adjset_names)]}
    
    #check if adjset_names is vector of character strings:
    if(is.character(adjset_names)==FALSE)stop("'adjset_names' is not a vector of character strings as required")
    
    #combine data on exposure (including already identified direct effect metabolites), metabolite outcome, and adjacency set, store as dataframe ("modeldata"), input data for glmulti:
    modeldata<-data.frame(cbind(met_out_data,exp_dat,adjset_data))
    colnames(modeldata)[1]<-met_outcome
    
    ###########################################################estimate multimodel coefficients#####################################################
    
    #modify fitting function of coxph so that it always includes one predefined set of variables (here: exposure + covariates), while subsetting adjacent metabolite set:
    coxph.redefined<-function(formula,data,always="",...){
      coxph(as.formula(paste(deparse(formula),always)),data=data,...)
    }
    
    #fit all possible causal models using glmulti: iterate over all possible combinations of adjacent metabolites, fixed exposure + covariates set
    
    glmulti_obj<-glmulti(y=name_surv_obj,                                      #response variable to be predicted: SURV-object
                         xr=c(adjset_names[1:length(adjset_names)]),            #predictor variables to be tested in iteration
                         data=modeldata,                                        #dataframe containing data
                         level=1,                                               #only main effects (terms of order 1) are used to build the candidate set
                         confsetsize = 2^(length(adjset_names)),                #number of models to be looked for, dependent on size of adjacency set
                         fitfunction=coxph.redefined,                           #fitting function to be used: here: coxph-function with fixed predictor set "always"
                         always=paste0("+",met_outcome,"+",always_set),         #defines "always"-set: complete set of exposures and covariates and current metabolite
                         includeobjects=TRUE,                                   #fitted models are included as objects
                         plotty=FALSE,                                          #no progress of the IC profile is plotted
                         report=FALSE                                           #report about progress at run time is given
    )
    
    #output summaries:
    
    glmulti_obj_objects<-glmulti_obj@objects                                    #collect glmulti objects, i.e. all fitted models (?)
    
    nbmds<-c(1:glmulti_obj@nbmods)                                              #number of fitted models of glmulti-function
    
    model_details<-list(NULL)
    
    for(j in seq(along=nbmds)){                                                 #generate list for model details
      model_details[[j]]<-list(NULL)
    }
    names(model_details)<-paste("Model",nbmds,sep="_")
    
    for(j in seq(along=nbmds)){                                                 #collect all model details for different adjacency subsets for specific metabolite
      model_summary<-glmulti_obj_objects[[j]]
      model_details[[j]]<-list(Model=paste("Model",j,"of",length(nbmds)),
                               Model_summary=model_summary)
    }
    
    model_details_all[[met_outcome]]<-list(Model_summaries=model_details,      #collect all model details for different adjacency subsets for specific metabolite in complete model list for all metabolites
                                           Number_of_Models=length(nbmds),
                                           Adj_set=paste(adjset_names,collapse = ", "),
                                           Outcome_metabolite=met_outcome)
    
  }  
  
  return(list(model_details_all=model_details_all,outcome=list(survival_obj=survival_obj,class=class(survival_obj))))
  
}



####################################################Extract exposure coefficients per outcome, i.e. metabolite#######################################################

getExp.coef.permetabolite<-function(object,metabolite){
  
  #object: output of net.coupler.out
  #metabolite: specific metabolite to evaluate
  
  #create vector containing integers from 1 to number of metabolite-specific models:
  nbm<-c(1:length(object$model_details_all[[metabolite]]$Model_summaries))
  Nbmds<-max(nbm)
  
  #define empty dataframes for output:
  mm_coef_temp1<-structure(list(Model=as.character(),
                                Nbmds=as.numeric(),
                                Outcome=as.character(),
                                Covariables=as.character(),
                                #NCov=as.numeric(),
                                HR=as.numeric(),
                                LCL=as.numeric(),
                                UCL=as.numeric(),
                                Beta=as.numeric(),
                                rSE=as.numeric(),
                                P=as.numeric()),
                           class="data.frame")
  
  mm_coef_temp2<-structure(list(Model=as.character(),
                                Nbmds=as.numeric(),
                                Outcome=as.character(),
                                Covariables=as.character(),
                                NCov=as.numeric(),
                                HR=as.numeric(),
                                LCL=as.numeric(),
                                UCL=as.numeric(),
                                Beta=as.numeric(),
                                rSE=as.numeric(),
                                P=as.numeric()),
                           class="data.frame")
  
  #loop along number of metabolite-specific models:
  for(i in seq(along=nbm)){
    #get metabolite-effect estimates (hazard-ratio, lower and upper confidence limit, beta, robust SE, and p-value) from single model and write into dataframe
    SUM<-summary(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model_summary)
    Cov<-rownames(SUM$coefficients)
    #Cov<-dplyr::filter(data.frame(row.names(SUM$coefficients)),row.names(SUM$coefficients)!=exposure)                     #record specific covariates of this model
    
    mm_coef_temp2<-data.frame(Model=as.character(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model),
                              Nbmds=Nbmds,
                              Outcome=metabolite,
                              Covariables=as.character(paste(Cov,collapse=", ")),
                              NCov=max(unlist(object$model_details_all[[metabolite]]$Model_summaries[[i]]$Model_summary$assign))-1,
                              HR=as.numeric(SUM$conf.int[[metabolite,1]]),
                              LCL=as.numeric(SUM$conf.int[[metabolite,3]]),
                              UCL=as.numeric(SUM$conf.int[[metabolite,4]]),
                              Beta=as.numeric(SUM$coefficients[[metabolite,1]]),
                              rSE=as.numeric(SUM$coefficients[[metabolite,4]]),
                              P=as.numeric(SUM$coefficients[[metabolite,6]]))
    
    #bind information to an outcome-specific dataframe:
    mm_coef_temp1<-bind_rows(mm_coef_temp1,mm_coef_temp2)
    
  }
  
  return(mm_coef_temp1)
  
}





####################################################Extract exposure coefficients for metabolite(s) on time-to-incidence###############################################

getExp.coef.out<-function(object,metabolite){
  
  #object: output of net.coupler
  #metabolite: specific metabolite to evaluate
  
  cat("*************************************************************************************************** \n")
  cat("This function produces a table of effect estimates of all (some) network-variables on an outcome    \n")
  cat("(time-to-event) for all possibles causal models based on conditional independence criteria encoded  \n")
  cat("in the input-network => MULTISET OF POSSIBLE EFFECTS PER EXPOSURE-OUTCOME PAIR; network-variables   \n")
  cat("of interest are selected by indicating the variable-names as character vector                       \n")
  cat("***************************************************************************************************")
  
  
  #define empty dataframes for output:
  
  mm_coef<-structure(list(Model=as.character(),
                          Nbmds=as.numeric(),
                          Outcome=as.character(),
                          Covariables=as.character(),
                          NCov=as.numeric(),
                          HR=as.numeric(),
                          LCL=as.numeric(),
                          UCL=as.numeric(),
                          Beta=as.numeric(),
                          rSE=as.numeric(),
                          P=as.numeric()),
                     class="data.frame")
  
  mm_coef_temp<-structure(list(Model=as.character(),
                               Nbmds=as.numeric(),
                               Outcome=as.character(),
                               Covariables=as.character(),
                               NCov=as.numeric(),
                               HR=as.numeric(),
                               LCL=as.numeric(),
                               UCL=as.numeric(),
                               Beta=as.numeric(),
                               rSE=as.numeric(),
                               P=as.numeric()),
                          class="data.frame")
  
  #for each metabolite: get metabolite coefficients from all possible models and write to table:
  
  for(j in metabolite){
    mm_coef_temp<-data.frame(getExp.coef.permetabolite(object=object,metabolite=j))
    mm_coef<-bind_rows(mm_coef,mm_coef_temp)
  }
  
  return(mm_coef)
  
}


#####################Summary statistics and multiple testing adjustment for net.coupler.out with survival object############################################################################


mult.stat.surv<-function(sum_netout,adjust_method,rule_1 = 12,rule_1_cut = 0.1,rule_2 = 5,rule_2_cut = 0.05,rule_3 = 15,rule_4 = 14, ass_rule1=16, round_number){
  
  #sum_netout: output of getExp.coef.out after adding original metabolite names
  #adjust_method: select adjustment method(s), e.g. "fdr" or c("fdr","bonferroni"), see p.adjust.methods for available options
  #rule= "set DE=1 or =2 if ..."
  #rule_1: determines first argument of rule (e.g. adjusted p-value of full model < rule_1_cut)
  #rule_2: determines second argument of rule (e.g. upperP < rule_2_cut)
  #rule_3: determines first variable of third argument of rule (e.g. highEst)
  #rule_4: determines second variable of third argument of rule (e.g. lowEst)
  #possible arguments for rule_...: adjusted p-values of full model=12=rule_1
  #                                 largest p-value across all adjacency models=upperP=5=rule_2
  #                                 highest estimated beta across all adjacency models=highEst=15=rule_3
  #                                 lowest estimated beta across all adjacency models=lowEst=14=rule_4
  #default rule: if(adjusted p< 0.1 & upperP<0.05 & (sign(highEst)==sign(lowEst)) & lowEst>0){DE=1}
  #              if(adjusted p< 0.1 & upperP<0.05 & (sign(highEst)==sign(lowEst)) & lowEst<0){DE=2}else{DE=0}
  #ass_rule1: determines first variable of second argument of association rule (e.g. bestGuess=16)
  #default association rule: if(adjusted p<0.1 & bestGuess>0){Assoc=1}
  #                          if(adjusted p<0.1 & bestGuess<0){Assoc=2}else{Assoc=0}
  #round_number: actual round number in ambiguous metabolites loop
  
  #return: for each metabolite: maximum number of models, mean hazard ratio, maximum HR, minimum HR, maximum p-value, minimum p-value, as well as number of covariates, HR and p-value for model with full adjustment:
  sum_netout_sum<-dplyr::left_join(sum_netout%>%dplyr::group_by(Outcome)%>%dplyr::summarise(avgHR=mean(HR),minHR=min(HR),maxHR=max(HR),upperP=max(P),lowerP=min(P),Nbmds=max(Nbmds)),
                                   sum_netout%>%dplyr::group_by(Outcome)%>%dplyr::filter(nchar(Covariables)==max(nchar(Covariables)))%>%dplyr::select(NCov=NCov,Outcome=Outcome,Metabolite=Metabolite,ind_HR=HR,ind_P=P),
                                   by="Outcome")
  
  #multiple testing correction:
  p_adjust_M<-p.adjust.methods[p.adjust.methods %in% adjust_method]                            #select multiple testing correction method(s)
  p_adj<-sapply(p_adjust_M,function(meth){p.adjust(sum_netout_sum$ind_P,meth)})            #calculate adjusted p-values for ind_P
  
  sum_netout_sum_adj<-bind_cols(sum_netout_sum,data.frame(p_adj))
  colnames(sum_netout_sum_adj)[12]<-adjust_method
  
  #add beta-coefficients:
  sum_netout_sum_adj<-sum_netout_sum_adj%>%dplyr::mutate(avgEst=log(avgHR),lowEst=log(minHR),highEst=log(maxHR),bestGuess=log(ind_HR))
  
  #is there an effect of specified metabolite on time-to-event?: determine effect-indicator (DE!=0: direct effect of metabolite on time-to-event, DE=0: ambiguous if adjusted ind_p<0.1):
  sum_netout_sum_adj<-mutate(sum_netout_sum_adj,
                             DE=derivedFactor(
                               "1"=(sum_netout_sum_adj[,rule_1]<rule_1_cut & sum_netout_sum_adj[,rule_2]<rule_2_cut & (sign(sum_netout_sum_adj[,rule_3])==sign(sum_netout_sum_adj[,rule_4])) & sum_netout_sum_adj[,rule_4]>0),
                               "2"=(sum_netout_sum_adj[,rule_1]<rule_1_cut & sum_netout_sum_adj[,rule_2]<rule_2_cut & (sign(sum_netout_sum_adj[,rule_3])==sign(sum_netout_sum_adj[,rule_4])) & sum_netout_sum_adj[,rule_4]<0),
                               .method="first",.default=0)
  )
  
  
  #is there any effect of specified metabolite on time-to-event?: determine association-indicator (Assoc!=0: there is an effect (no differenciation between direct or indirect), Assoc=0: there is no effect)
  sum_netout_sum_adj<-mutate(sum_netout_sum_adj,
                             Assoc=derivedFactor(
                               "1"=(sum_netout_sum_adj[,rule_1]<rule_1_cut & sum_netout_sum_adj[,ass_rule1]>0),
                               "2"=(sum_netout_sum_adj[,rule_1]<rule_1_cut & sum_netout_sum_adj[,ass_rule1]<0),
                               .method="first",.default=0)
  )
  
  #collect ambiguous effects:
  amb<-sum_netout_sum_adj%>%filter(Assoc!=0 & DE==0)
  
  #collect direct effects:
  direct<-sum_netout_sum_adj%>%filter(Assoc!=0 & DE!=0)
  
  #collect summary statistics including round information
  sum_netout_sum_adj_FV<-sum_netout_sum_adj%>%dplyr::mutate(round=as.numeric(round_number),Assoc_FV=as.character(Assoc),DE_FV=as.character(DE))
  
  return(list(sum_netout_sum_adj=sum_netout_sum_adj,amb=amb,direct=direct,sum_netout_sum_adj_FV=sum_netout_sum_adj_FV))
  
}



################################################ambiguous metabolites loop for net.coupler.out. survival analysis###############################################################################################################



amb.met.loop.out.surv<-function(exp_dat,graph_skel,dat,dat_compl,DE,survival_obj,always_set,met_map,adjust_method,round_number){
  
  #exp_dat: exposure/phenotype data
  #graph_skel: estimated DAG skeleton of samples x metabolites data matrix
  #dat: renamed samples x metabolites data matrix
  #dat_compl: renamed samples x metabolites data matrix (typically same as dat)
  #DE: indicator if direct effects were already identified, here set to NA
  #survival_obj: "survival" object
  #always_set: fixed set of covariates always included in model
  #met_map: mapping information between old and new metabolite names
  #adjust_method: select adjustment method(s), e.g. "fdr" or c("fdr","bonferroni"), see p.adjust.methods for available options
  #round_number: actual round number in ambiguous metabolites loop, initiate with round_number=1
  
  netout_amb<-list()
  netout_direct<-list()
  netout_sum<-list()
  
  repeat{
    
    print(paste(round_number,". iteration",sep=""))
    
    #estimate direct effects of metabolites on time-to-diabetes-incident: models are adjusted for all possible combinations of direct neighbors (==variables in adjacency set) -> Output is multiset of possible effects:
    net_coupler_out<-net.coupler.out(graph_skel=graph_skel,dat=dat,dat_compl = dat_compl,exp_dat = exp_dat,DE=DE,survival_obj = survival_obj,always_set = always_set,name_surv_obj=deparse(substitute(survival_obj)))       
    
    #return results:
    sum_netout<-getExp.coef.out(object=net_coupler_out,metabolite=colnames(dat))
    
    #get original metabolite names back:
    sum_netout<-merge(sum_netout,as.data.frame(met_map),by="Outcome")
    
    sum_stat_netout<-mult.stat.surv(sum_netout=sum_netout,adjust_method = adjust_method,round_number = round_number)       #calculate summary statistics and determine direct and ambiguous effects
    
    netout_amb[[length(netout_amb)+1]]<-sum_stat_netout$amb                    #summary statistics for metabolites which have ambiguous effect on time-to-diabetes-incident
    netout_direct[[length(netout_direct)+1]]<-sum_stat_netout$direct              #summary statistics for metabolites which have direct effect on time-to-diabetes-incident
    netout_sum[[length(netout_sum)+1]]<-sum_stat_netout$sum_netout_sum_adj_FV   #summary statistics for meatbolites including round-number information
    names(netout_amb)[length(netout_amb)]<-paste0(round_number,". iteration")
    names(netout_direct)[length(netout_direct)]<-paste0(round_number,". iteration")
    names(netout_sum)[length(netout_sum)]<-paste0(round_number,". iteration")
    
    numb_DE<-dim(netout_direct[[length(netout_direct)]])[1]
    numb_AMB<-dim(netout_amb[[length(netout_amb)]])[1]
    
    #check if direct effects!=0 and/or ambiguous effects!=0, otherwise quit:
    if(numb_DE==0 | numb_AMB==0){cat(paste("\n ***no direct and/or ambiguous effects identified -> stop algorithm after ",round_number,". iteration*** \n",sep="")); break}
    
    #if there are direct and ambiguous effects of metabolites on time-to-diabetes-incident -> continue:
    
    #names of direct effect metabolites:
    DE<-c(DE,netout_direct[[length(netout_direct)]]$Outcome)
    DE_1<-netout_direct[[length(netout_direct)]]$Outcome
    
    print("here1")
    
    #define new metabolite data matrix only including ambiguous effect metabolites:
    dat<-as.matrix(dat[,c(which(colnames(dat) %in% netout_amb[[length(netout_amb)]]$Outcome)),drop=FALSE])
    
    print("here2")
    print(colnames(exp_dat))
    
    #add direct effect metabolite data to exp_dat:
    exp_dat<-cbind(dat_compl[,c(which(colnames(dat_compl) %in% DE_1))],exp_dat)
    colnames(exp_dat)[1:length(DE_1)]<-DE_1
    
    print(colnames(exp_dat))
    
    #add direct effect metabolites to always set:
    always_set<-paste0(paste0(DE_1,collapse = " + "),sep=" + ",always_set)
    
    print(always_set)
    
    round_number<-round_number+1
    
  }
  
  return(list(netout_direct=netout_direct,netout_amb=netout_amb,netout_sum=netout_sum))
  
}





####################################################Data import and analysis###########################################################################


#load metabolite data (e.g. aaPCs) for complete cohort:
gpl_data <- readRDS("C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/R_data/STR_GPL.rds")

met_data<-dplyr::select(gpl_data,contains("_aa_"))

dim(met_data)                             #case-cohort: subcohort + all incident diabetes cases in complete cohort
#[1] 2731   34

#metabolite data (e.g. aaPCs) for subcohort:
gpl_data_SC <- readRDS("C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/R_data/STR_GPL_SC.rds")

met_data_SC<-dplyr::select(gpl_data_SC,contains("_aa_"))

dim(met_data_SC)       #metabolite data with rows=patients and columns=metabolites
#[1] 2092   34

is.numeric(as.matrix(met_data_SC))
#[1] TRUE

#load type-2 diabetes incident times:
load("C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/R_data/T2D_data")

dim(T2D_data)
#[1] 2731   26                                    #phenotype data including time information for 2731 patients, including diet information, lifestyle, etc.

colnames(T2D_data)
# [1] "subcohort"  "SEX"        "fall"       "Med_Hypert" "Med_HLipid" "sport"      "bike"       "alk_1"      "alk_2"      "alk_3"      "alk_4"      "alk_5"      "alk_6"     
# [14] "sta_time"   "sto_time"   "WGBperMJ"   "TMperMJ"    "age_years"  "smk1"       "smk2"       "smk3"       "educ1"      "educ2"      "educ3"      "CofCup"     "ID"

#load phenotype data for complete cohort:
Exp_data<-readRDS("C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/R_data/EXP_DATA.rds")

dim(Exp_data)
#[1] 2731   56                                   #phenotype data for 2731 patients, including diet information, lifestyle, etc.


met_data_SC_rename<-rename.met(dat=met_data_SC)$data_renamed                               #rename metabolites with short names
met_mapping_SC<-rename.met(dat=met_data_SC)$names_mapping                                     #mapping information between old and new metabolite names

met_data_rename<-rename.met(dat=met_data)$data_renamed                               #rename metabolites with short names
met_mapping<-rename.met(dat=met_data)$names_mapping                                     #mapping information between old and new metabolite names

#skeleton estimation for metabolite matrix subcohort:
met_SC_estimates_norename<-est.pcor.skel.DAG.adj(dat=met_data_SC)                       
met_SC_skel_norename<-met_SC_estimates_norename$skel_est                                                   #estimate DAG skeleton for non-renamed matrix
met_SC_DAG_norename<-met_SC_estimates_norename$DAG_est                                                     #estimate DAG for non-renamed matrix
met_SC_adj_mat_norename<-met_SC_estimates_norename$adj_matrix                                              #estimate adjacency matrix for non-renamed matrix

met_SC_estimates_rename<-est.pcor.skel.DAG.adj(dat=met_data_SC_rename)
met_SC_skel_rename<-met_SC_estimates_rename$skel_est                                                       #estimate DAG skeleton for renamed matrix
met_SC_DAG_rename<-met_SC_estimates_rename$DAG_est                                                         #estimate DAG for renamed matrix
met_SC_adj_mat_rename<-met_SC_estimates_rename$adj_matrix                                                  #estimate adjacency matrix for renamed matrix

#skeleton estimation for complete cohort metabolite matrix:
met_estimates_norename<-est.pcor.skel.DAG.adj(dat=met_data)                       
met_skel_norename<-met_estimates_norename$skel_est                                                   #estimate DAG skeleton for non-renamed matrix
met_DAG_norename<-met_estimates_norename$DAG_est                                                     #estimate DAG for non-renamed matrix
met_adj_mat_norename<-met_estimates_norename$adj_matrix                                              #estimate adjacency matrix for non-renamed matrix

met_estimates_rename<-est.pcor.skel.DAG.adj(dat=met_data_rename)
met_skel_rename<-met_estimates_rename$skel_est                                                       #estimate DAG skeleton for renamed matrix
met_DAG_rename<-met_estimates_rename$DAG_est                                                         #estimate DAG for renamed matrix
met_adj_mat_rename<-met_estimates_rename$adj_matrix                                                  #estimate adjacency matrix for renamed matrix

#create "survival" object:
t2d_surv<-Surv(T2D_data$sta_time,T2D_data$sto_time,T2D_data$fall)

#define "always"-set:
always_set<-paste0(colnames(Exp_data)[-(which(colnames(Exp_data) %in% c("SEX","subcohort","ID","age","age_years")))],collapse = " + ")
always_set<-paste0(always_set," + cluster(ID) + strata(age_years)")                                   #cluster(ID) specific for case-cohort design, strata(age_years) stratification of baseline risk according to age in years



#estimate direct effects of metabolites on time-to-diabetes-incident:
tic()
amb_met_loop_out<-amb.met.loop.out.surv(exp_dat=Exp_data,graph_skel=met_SC_skel_rename,dat=met_data_rename,dat_compl=met_data_rename,DE=NULL,survival_obj=t2d_surv,always_set=always_set,met_map=met_mapping,adjust_method="fdr",round_number=1)
toc()

#time: 4.2min

#final results:

#complete results:
net_coupler_out_iteration1<-amb_met_loop_out$netout_sum$`1. iteration`

net_coupler_out_iteration2<-amb_met_loop_out$netout_sum$`2. iteration`

net_coupler_out_iteration3<-amb_met_loop_out$netout_sum$`3. iteration`

netout_sum_1_2_3_final<-Reduce(function(x,y,z)merge(x,y,all=TRUE,by="Outcome"),list(net_coupler_out_iteration1,net_coupler_out_iteration2,net_coupler_out_iteration3))



write.xlsx(netout_sum_1_2_3_final,file="C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/HZ_netout_aaPC_ambloop_20_7_2017.xls")

#all direct effects:
net_coupler_out_direct1<-amb_met_loop_out$netout_direct$`1. iteration`
net_coupler_out_direct2<-amb_met_loop_out$netout_direct$`2. iteration`

net_coupler_out_direct_final<-rbind(net_coupler_out_direct1,net_coupler_out_direct2)

net_coupler_out_direct_final
# A tibble: 5 x 18
# Outcome     avgHR     minHR     maxHR       upperP       lowerP Nbmds  NCov   Metabolite    ind_HR        ind_P          fdr     avgEst     lowEst    highEst  bestGuess     DE  Assoc
# <chr>     <dbl>     <dbl>     <dbl>        <dbl>        <dbl> <dbl> <dbl>       <fctr>     <dbl>        <dbl>        <dbl>      <dbl>      <dbl>      <dbl>      <dbl> <fctr> <fctr>
#   1    NM20 1.4614121 1.3318502 1.5393148 0.0014633373 2.543242e-07    16    55 rPC_aa_C38_3 1.5393148 5.648558e-04 4.801274e-03  0.3794032  0.2865691  0.4313374  0.4313374      1      1
# 2    NM33 1.3058133 1.1086547 1.4776870 0.0458789204 1.160349e-05     8    54 rPC_aa_C42_5 1.4776870 1.603423e-04 1.817212e-03  0.2668261  0.1031473  0.3904780  0.3904780      1      1
# 3    NM22 0.6029201 0.4688673 0.6966110 0.0002019488 6.909858e-06     4    55 rPC_aa_C38_5 0.4688673 6.909858e-06 6.909858e-05 -0.5059705 -0.7574356 -0.3615282 -0.7574356      2      2
# 4    NM30 0.6854735 0.6126198 0.7669771 0.0150326314 3.380629e-13    16    57 rPC_aa_C42_1 0.7289033 5.043250e-03 1.850699e-02 -0.3776455 -0.4900107 -0.2652984 -0.3162142      2      2
# 5     NM6 0.7999596 0.7485887 0.8481255 0.0247436411 7.576171e-05     8    56 rPC_aa_C32_3 0.8218709 2.474364e-02 3.643585e-02 -0.2231941 -0.2895656 -0.1647267 -0.1961719      2      2

#all final ambiguous effects:
amb_met_loop_out$netout_amb$`3. iteration`
# A tibble: 5 x 18
# Outcome     avgHR     minHR     maxHR    upperP      lowerP Nbmds  NCov   Metabolite    ind_HR       ind_P        fdr     avgEst      lowEst     highEst  bestGuess     DE  Assoc
# <chr>     <dbl>     <dbl>     <dbl>     <dbl>       <dbl> <dbl> <dbl>       <fctr>     <dbl>       <dbl>      <dbl>      <dbl>       <dbl>       <dbl>      <dbl> <fctr> <fctr>
#   1    NM17 0.8065255 0.6110597 0.9657696 0.7690416 0.005109276    16    60 rPC_aa_C36_6 0.6294707 0.009917459 0.02975238 -0.2150198 -0.49256057 -0.03483002 -0.4628759      0      2
# 2    NM24 1.1859421 1.0342496 1.3427333 0.6812973 0.008297387     8    59 rPC_aa_C40_2 1.3060686 0.018227520 0.03470137  0.1705375  0.03367615  0.29470734  0.2670216      0      1
# 3    NM25 0.8159863 0.7386868 0.9051769 0.2450431 0.005109764     8    59 rPC_aa_C40_3 0.7405413 0.005569742 0.02975238 -0.2033578 -0.30288131 -0.09962484 -0.3003738      0      2
# 4    NM29 0.8289516 0.8009908 0.8703262 0.1695697 0.036505825     8    59 rPC_aa_C42_0 0.8152448 0.056753682 0.06810442 -0.1875935 -0.22190579 -0.13888723 -0.2042669      0      2
# 5    NM31 0.8471528 0.7909966 0.9189214 0.3306723 0.014083141    32    61 rPC_aa_C42_2 0.8036937 0.023134249 0.03470137 -0.1658742 -0.23446156 -0.08455471 -0.2185371      0      2



