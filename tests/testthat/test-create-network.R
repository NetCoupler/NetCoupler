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

test_that("network is constructed even with missingness", {
    metabolite_network_na <- simulated_data %>%
        insert_random_missingness() %>%
        nc_estimate_network(starts_with("metabolite"))

    expect_identical(class(metabolite_network_na)[1], "tbl_graph")

    # For number of neighbours, etc
    edges <- as_edge_tbl(metabolite_network_na)
    edges <- unique(c(edges$source_node, edges$target_node))
    expect_true(all(edges %in% names(metabolite_data)))
})

# test_that("edge table generates correct output", {
#     # TODO: Fill this out
# })
