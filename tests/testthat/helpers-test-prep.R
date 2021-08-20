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

insert_random_missingness <- function(data) {
    purrr::map_df(data, ~ .[sample(
        c(TRUE, NA),
        prob = c(0.98, 0.02),
        size = length(.),
        replace = TRUE
    )])
}
