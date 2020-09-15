suppressPackageStartupMessages(library(dplyr, quietly = TRUE))

# Just in case
set.seed(21451)

metabolite_network <- simulated_data %>%
    nc_standardize(starts_with("metabolite")) %>%
    nc_estimate_network(starts_with("metabolite"))

expect_correct_model_results <- function(.actual, .expected_metabolites) {
    expect_type(.actual, "list")
    expect_identical(class(.actual)[1], "tbl_df")
    expect_identical(sort(.actual$index_node),
                     sort(.expected_metabolites)[-3])
}
