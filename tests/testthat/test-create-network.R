context("Create metabolic variable network")

# Making partial independence network from metabolite data
metabolite_data <- simulated_data %>%
    select(starts_with("metabolite"))

metabolite_network <- metabolite_data %>%
    nc_standardize() %>%
    nc_estimate_network()

test_that("network is created", {
    # For number of edges
    # Not sure how to test number of edges properly, since the set number in
    # the simulated data doesn't mean the network estimation will have the same
    # amount.
    # expect_equal(nrow(as_edge_tbl(metabolite_network)), 13)
    # For columns
    expect_equal(ncol(as_edge_tbl(metabolite_network)), 3)

    # For data frame
    expect_error(
        metabolite_data$metabolite_1 %>%
            nc_estimate_network()
    )

    # For alpha number
    expect_error(
        metabolite_data %>%
            nc_estimate_network("0.05")
    )
})

# test_that("edge table generates correct output", {
#     # TODO: Fill this out
# })
