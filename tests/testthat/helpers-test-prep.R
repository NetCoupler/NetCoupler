suppressPackageStartupMessages(library(dplyr, quietly = TRUE))

metabolite_network <- simulated_data %>%
    select(matches("metabolite")) %>%
    nc_estimate_network()

