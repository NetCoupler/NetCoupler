context("Create metabolic variable network")

test_that("network is created", {
    # Make partial independence network from metabolite data
    metabolite_network <- simulated_data %>%
        dplyr::select(dplyr::matches("metabolite")) %>%
        nc_create_network()

    expect_type(metabolite_network, "S4")
})
