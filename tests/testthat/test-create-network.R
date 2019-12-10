context("Create metabolic variable network")

test_that("network is created", {
    # Make partial independence network from metabolite data
    metabolite_data <- simulated_data %>%
        dplyr::select(dplyr::matches("metabolite"))

    metabolite_network <- metabolite_data %>%
        nc_create_network()

    expect_type(metabolite_network, "S4")

    # For data frame
    expect_error(
        metabolite_data$metabolite_1 %>%
            nc_create_network()
    )

    # For alpha number
    expect_error(
        metabolite_data %>%
            nc_create_network("0.05")
    )
})
