#
#
# shorter analysis version with loop for ambiguous metabolites#
#

# Data import and analysis#

# load metabolite data (e.g. SM) for complete cohort:
gpl_data <- readRDS("H:/Metabolomics/DISS I/R_obj/T2D/STR_GPL.rds")

met_data <- dplyr::select(gpl_data, contains("SM"))

dim(met_data) # case-cohort: subcohort + all incident diabetes cases in complete cohort


# metabolite data (e.g. SM) for subcohort:
gpl_data_SC <- readRDS("H:/Metabolomics/DISS I/R_Obj/GPL/STR_GPL_SC.rds")

met_data_SC <- dplyr::select(gpl_data_SC, contains("SM"))

dim(met_data_SC) # metabolite data with rows=patients and columns=metabolites


is.numeric(as.matrix(met_data_SC))


# load type-2 diabetes incident times:
load("H:/Metabolomics/DISS I/R_obj/T2D/T2D_data")

dim(T2D_data)


colnames(T2D_data)



# load phenotype data for complete cohort:
Exp_data <- readRDS("H:/Metabolomics/DISS I/R_obj/T2D/EXP_DATA.rds")

dim(Exp_data)


met_data_SC_rename <- rename.met(dat = met_data_SC)$data_renamed # rename metabolites with short names
met_mapping_SC <- rename.met(dat = met_data_SC)$names_mapping # mapping information between old and new metabolite names

met_data_rename <- rename.met(dat = met_data)$data_renamed # rename metabolites with short names
met_mapping <- rename.met(dat = met_data)$names_mapping # mapping information between old and new metabolite names

# skeleton estimation for metabolite matrix subcohort:
met_SC_estimates_norename_005 <- est.pcor.skel.DAG.adj(dat = met_data_SC, alpha_val = 0.05)
met_SC_skel_norename_005 <- met_SC_estimates_norename_005$skel_est # estimate DAG skeleton for non-renamed matrix
met_SC_DAG_norename_005 <- met_SC_estimates_norename_005$DAG_est # estimate DAG for non-renamed matrix
met_SC_adj_mat_norename_005 <- met_SC_estimates_norename_005$adj_matrix # estimate adjacency matrix for non-renamed matrix

met_SC_estimates_rename_005 <- est.pcor.skel.DAG.adj(dat = met_data_SC_rename, alpha_val = 0.05)
met_SC_skel_rename_005 <- met_SC_estimates_rename_005$skel_est # estimate DAG skeleton for renamed matrix
met_SC_DAG_rename_005 <- met_SC_estimates_rename_005$DAG_est # estimate DAG for renamed matrix
met_SC_adj_mat_rename_005 <- met_SC_estimates_rename_005$adj_matrix # estimate adjacency matrix for renamed matrix

# skeleton estimation for complete cohort metabolite matrix:
met_estimates_norename_005 <- est.pcor.skel.DAG.adj(dat = met_data, alpha_val = 0.05)
met_skel_norename_005 <- met_estimates_norename_005$skel_est # estimate DAG skeleton for non-renamed matrix
met_DAG_norename_005 <- met_estimates_norename_005$DAG_est # estimate DAG for non-renamed matrix
met_adj_mat_norename_005 <- met_estimates_norename_005$adj_matrix # estimate adjacency matrix for non-renamed matrix

met_estimates_rename_005 <- est.pcor.skel.DAG.adj(dat = met_data_rename, alpha_val = 0.05)
met_skel_rename_005 <- met_estimates_rename_005$skel_est # estimate DAG skeleton for renamed matrix
met_DAG_rename_005 <- met_estimates_rename_005$DAG_est # estimate DAG for renamed matrix
met_adj_mat_rename_005 <- met_estimates_rename_005$adj_matrix # estimate adjacency matrix for renamed matrix

# create "survival" object:
t2d_surv <- Surv(T2D_data$sta_time, T2D_data$sto_time, T2D_data$fall)

# define "always"-set:
always_set <- paste0(colnames(Exp_data)[-(which(colnames(Exp_data) %in% c("SEX", "subcohort", "ID", "age", "age_years")))], collapse = " + ")
always_set <- paste0(always_set, " + cluster(ID) + strata(age_years)") # cluster(ID) specific for case-cohort design, strata(age_years) stratification of baseline risk according to age in years

# estimate direct effects of metabolites on time-to-diabetes-incident:
tic()
amb_met_loop_out_005 <- amb.met.loop.out.surv(exp_dat = Exp_data, graph_skel = met_SC_skel_rename_005, dat = met_data_rename, dat_compl = met_data_rename, DE = NULL, survival_obj = t2d_surv, always_set = always_set, met_map = met_mapping, adjust_method = "fdr", round_number = 1)
toc()

# final results:

# complete results:
net_coupler_out_iteration1_005 <- amb_met_loop_out_005$netout_sum$`1. iteration`

net_coupler_out_iteration2_005 <- amb_met_loop_out_005$netout_sum$`2. iteration`

net_coupler_out_iteration3_005 <- amb_met_loop_out_005$netout_sum$`3. iteration`

netout_sum_1_2_3_final_005 <- Reduce(function(x, y, z) merge(x, y, all = TRUE, by = "Outcome"), list(net_coupler_out_iteration1_005, net_coupler_out_iteration2_005, net_coupler_out_iteration3_005))

# write.xlsx(netout_sum_1_2_3_final_005,file="C:/Users/helena.zacharias/Documents/Helmholtz/KORA_stress/data_analysis/Clemens_netcoupler/HZ_netout_SM_ambloop_005_18_9_2017.xls")

# all direct effects:
net_coupler_out_direct1_005 <- amb_met_loop_out_005$netout_direct$`1. iteration`
net_coupler_out_direct2_005 <- amb_met_loop_out_005$netout_direct$`2. iteration`
net_coupler_out_direct3_005 <- amb_met_loop_out_005$netout_direct$`3. iteration`

net_coupler_out_direct_final_005 <- rbind(net_coupler_out_direct1_005, net_coupler_out_direct2_005, net_coupler_out_direct3_005)

net_coupler_out_direct_final_005

# all final ambiguous effects:
amb_met_loop_out_005$netout_amb$`3. iteration`
