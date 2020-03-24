context("Create metabolic variable network")

metabolite_data <- simulated_data %>%
    dplyr::select(dplyr::matches("metabolite"))

metabolite_network <- metabolite_data %>%
    nc_create_network()

test_that("network is created", {
    # Making partial independence network from metabolite data

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

test_that("adjacency graph object is constructed", {
    # Based on adjacency matrix and partial correlation matrix
    adj_graph <- nc_adjacency_graph(metabolite_data,
                                    metabolite_network)

    expect_identical(class(adj_graph), "igraph")
})
