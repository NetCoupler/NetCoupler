

#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
##############################################################shorter analysis version with loop for ambiguous metabolites#######################################################################################################
#################################################################################################################################################################################################################################

rm(list=ls())


####################################################Data import and analysis###########################################################################

#load metabolite data (e.g. SM) for complete cohort:
gpl_data <- readRDS("H:/Metabolomics/DISS I/R_obj/T2D/STR_GPL.rds")

met_data<-dplyr::select(gpl_data,contains("SM"))

dim(met_data)                             #case-cohort: subcohort + all incident diabetes cases in complete cohort
#[1] 2731   14

#metabolite data (e.g. SM) for subcohort:
gpl_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/GPL/STR_GPL_SC.rds")

met_data_SC<-dplyr::select(gpl_data_SC,contains("SM"))

dim(met_data_SC)       #metabolite data with rows=patients and columns=metabolites
#[1] 2092   14

is.numeric(as.matrix(met_data_SC))
#[1] TRUE

#load type-2 diabetes incident times:
load("H:/Metabolomics/DISS I/R_obj/T2D/T2D_data")

dim(T2D_data)
#[1] 2731   26                                    #phenotype data including time information for 2731 patients, including diet information, lifestyle, etc.

colnames(T2D_data)
# [1] "subcohort"  "SEX"        "fall"       "Med_Hypert" "Med_HLipid" "sport"      "bike"       "alk_1"      "alk_2"      "alk_3"      "alk_4"      "alk_5"      "alk_6"
# [14] "sta_time"   "sto_time"   "WGBperMJ"   "TMperMJ"    "age_years"  "smk1"       "smk2"       "smk3"       "educ1"      "educ2"      "educ3"      "CofCup"     "ID"

#load phenotype data for complete cohort:
Exp_data<- readRDS("H:/Metabolomics/DISS I/R_obj/T2D/EXP_DATA.rds")

dim(Exp_data)
#[1] 2731   56                                   #phenotype data for 2731 patients, including diet information, lifestyle, etc.


met_data_SC_rename<-rename.met(dat=met_data_SC)$data_renamed                               #rename metabolites with short names
met_mapping_SC<-rename.met(dat=met_data_SC)$names_mapping                                     #mapping information between old and new metabolite names

met_data_rename<-rename.met(dat=met_data)$data_renamed                               #rename metabolites with short names
met_mapping<-rename.met(dat=met_data)$names_mapping                                     #mapping information between old and new metabolite names

#skeleton estimation for metabolite matrix subcohort:
met_SC_estimates_norename_005<-est.pcor.skel.DAG.adj(dat=met_data_SC,alpha_val = 0.05)
met_SC_skel_norename_005<-met_SC_estimates_norename_005$skel_est                                                   #estimate DAG skeleton for non-renamed matrix
met_SC_DAG_norename_005<-met_SC_estimates_norename_005$DAG_est                                                     #estimate DAG for non-renamed matrix
met_SC_adj_mat_norename_005<-met_SC_estimates_norename_005$adj_matrix                                              #estimate adjacency matrix for non-renamed matrix

met_SC_estimates_rename_005<-est.pcor.skel.DAG.adj(dat=met_data_SC_rename,alpha_val = 0.05)
met_SC_skel_rename_005<-met_SC_estimates_rename_005$skel_est                                                       #estimate DAG skeleton for renamed matrix
met_SC_DAG_rename_005<-met_SC_estimates_rename_005$DAG_est                                                         #estimate DAG for renamed matrix
met_SC_adj_mat_rename_005<-met_SC_estimates_rename_005$adj_matrix                                                  #estimate adjacency matrix for renamed matrix

#skeleton estimation for complete cohort metabolite matrix:
met_estimates_norename_005<-est.pcor.skel.DAG.adj(dat=met_data,alpha_val = 0.05)
met_skel_norename_005<-met_estimates_norename_005$skel_est                                                   #estimate DAG skeleton for non-renamed matrix
met_DAG_norename_005<-met_estimates_norename_005$DAG_est                                                     #estimate DAG for non-renamed matrix
met_adj_mat_norename_005<-met_estimates_norename_005$adj_matrix                                              #estimate adjacency matrix for non-renamed matrix

met_estimates_rename_005<-est.pcor.skel.DAG.adj(dat=met_data_rename,alpha_val = 0.05)
met_skel_rename_005<-met_estimates_rename_005$skel_est                                                       #estimate DAG skeleton for renamed matrix
met_DAG_rename_005<-met_estimates_rename_005$DAG_est                                                         #estimate DAG for renamed matrix
met_adj_mat_rename_005<-met_estimates_rename_005$adj_matrix                                                  #estimate adjacency matrix for renamed matrix

#create "survival" object:
t2d_surv<-Surv(T2D_data$sta_time,T2D_data$sto_time,T2D_data$fall)

#define "always"-set:
always_set<-paste0(colnames(Exp_data)[-(which(colnames(Exp_data) %in% c("SEX","subcohort","ID","age","age_years")))],collapse = " + ")
always_set<-paste0(always_set," + cluster(ID) + strata(age_years)")                                   #cluster(ID) specific for case-cohort design, strata(age_years) stratification of baseline risk according to age in years



#estimate direct effects of metabolites on time-to-diabetes-incident:
tic()
amb_met_loop_out_005<-amb.met.loop.out.surv(exp_dat=Exp_data,graph_skel=met_SC_skel_rename_005,dat=met_data_rename,dat_compl=met_data_rename,DE=NULL,survival_obj=t2d_surv,always_set=always_set,met_map=met_mapping,adjust_method="fdr",round_number=1)
toc()

#final results:

#complete results:
net_coupler_out_iteration1_005<-amb_met_loop_out_005$netout_sum$`1. iteration`

net_coupler_out_iteration2_005<-amb_met_loop_out_005$netout_sum$`2. iteration`

net_coupler_out_iteration3_005<-amb_met_loop_out_005$netout_sum$`3. iteration`

netout_sum_1_2_3_final_005<-Reduce(function(x,y,z)merge(x,y,all=TRUE,by="Outcome"),list(net_coupler_out_iteration1_005,net_coupler_out_iteration2_005,net_coupler_out_iteration3_005))



#write.xlsx(netout_sum_1_2_3_final_005,file="C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/HZ_netout_SM_ambloop_005_18_9_2017.xls")

#all direct effects:
net_coupler_out_direct1_005<-amb_met_loop_out_005$netout_direct$`1. iteration`
net_coupler_out_direct2_005<-amb_met_loop_out_005$netout_direct$`2. iteration`
net_coupler_out_direct3_005<-amb_met_loop_out_005$netout_direct$`3. iteration`

net_coupler_out_direct_final_005<-rbind(net_coupler_out_direct1_005,net_coupler_out_direct2_005,net_coupler_out_direct3_005)

net_coupler_out_direct_final_005
# A tibble: 7 x 18
# Outcome     avgHR     minHR     maxHR       upperP       lowerP Nbmds  NCov     Metabolite    ind_HR        ind_P          fdr     avgEst     lowEst    highEst  bestGuess     DE  Assoc
# <chr>     <dbl>     <dbl>     <dbl>        <dbl>        <dbl> <dbl> <dbl>         <fctr>     <dbl>        <dbl>        <dbl>      <dbl>      <dbl>      <dbl>      <dbl> <fctr> <fctr>
#   1    NM10 0.8529664 0.8217661 0.8742225 0.0298472387 7.165044e-04     4    53      rSM_C20_2 0.8684133 2.379062e-02 4.758125e-02 -0.1590352 -0.1962995 -0.1344204 -0.1410875      2      2
# 2    NM12 0.7253865 0.6335816 0.8110995 0.0162965089 1.190799e-07     8    54      rSM_C24_1 0.7096051 1.879896e-03 6.579636e-03 -0.3210506 -0.4563665 -0.2093646 -0.3430467      2      2
# 3     NM7 0.7487158 0.6746416 0.8370972 0.0222113346 6.472914e-05    16    55      rSM_C16_1 0.7002194 4.003563e-03 9.341647e-03 -0.2893959 -0.3935737 -0.1778151 -0.3563615      2      2
# 4     NM3 2.1384723 1.3372969 3.0875399 0.0060942211 2.365907e-11    16    58 rSM__OH__C22_1 2.9794045 2.605694e-07 5.211388e-07  0.7600917  0.2906504  1.1273746  1.0917234      1      1
# 5     NM8 1.9840487 1.4258317 2.6258792 0.0001075808 9.144132e-10    16    58      rSM_C18_0 2.5521299 9.017145e-08 3.606858e-07  0.6851395  0.3547553  0.9654158  0.9369283      1      1
# 6     NM4 0.4939919 0.4409137 0.5593462 0.0001239414 9.196026e-06     4    58 rSM__OH__C22_2 0.4409137 3.534616e-05 7.069232e-05 -0.7052362 -0.8189062 -0.5809866 -0.8189062      2      2
# 7     NM9 0.6246849 0.6246849 0.6246849 0.0021135068 2.113507e-03     1    56      rSM_C18_1 0.6246849 2.113507e-03 2.113507e-03 -0.4705080 -0.4705080 -0.4705080 -0.4705080      2      2


#all final ambiguous effects:
amb_met_loop_out_005$netout_amb$`3. iteration`
# A tibble: 0 x 18
# ... with 18 variables: Outcome <chr>, avgHR <dbl>, minHR <dbl>, maxHR <dbl>, upperP <dbl>, lowerP <dbl>, Nbmds <dbl>, NCov <dbl>, Metabolite <fctr>, ind_HR <dbl>, ind_P <dbl>, fdr <dbl>, avgEst <dbl>, lowEst <dbl>,
#   highEst <dbl>, bestGuess <dbl>, DE <fctr>, Assoc <fctr>
