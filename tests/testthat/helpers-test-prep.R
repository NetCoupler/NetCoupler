suppressPackageStartupMessages(library(dplyr, quietly = TRUE))

metabolite_network <- simulated_data %>%
    nc_estimate_network(starts_with("metabolite"))

