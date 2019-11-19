#######################
# TEST nc_make_network
# Clemens Wittenbecher
library(dagitty)
library(dplyr)
library(survival)
library(igraph)

# simulate data 
# functions generated in the simulating-DAG.R script
simulated_data <-
  bind_cols(
    simulate_survival_time(
      .data = simulated_dag_data,
      IV1 = "metabolite_3",
      IV2 = "metabolite_9",
      IV3 = "metabolite_12"
    ),
    simulated_dag_data
  )

# get a vector of variable names defining the data to generate the network model
node_names<-colnames(dplyr::tbl_df(simulated_data)%>%dplyr::select(contains("metabolite")))

# glmulti: number of strings that can be processed in the iterative multimodel procudre is restricted
# rename network variables with short names - this is necessary for later steps
renaming<-rename.met(simulated_data,node_names)
renamed_simulated_data<-renaming[[1]]
matching_table_names_newnames<-renaming[[2]]
nodes_short_names<-colnames(dplyr::tbl_df(renamed_simulated_data)%>%dplyr::select(contains("NM")))



# Make partial independence network from metabolite data
nc_make_network1(renamed_simulated_data,.05,nodes_short_names)
nc_make_network1(simulated_data,.05,node_names)
