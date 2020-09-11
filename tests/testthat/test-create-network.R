context("Create metabolic variable network")

# Making partial independence network from metabolite data
metabolite_data <- simulated_data %>%
    dplyr::select(dplyr::matches("metabolite"))

metabolite_network <- metabolite_data %>%
    nc_estimate_network()

test_that("network is created", {
    # TODO: This might not always be this class.
    expect_s4_class(metabolite_network, "pcAlgo")

    # For number of neighbours, etc
    edges <- metabolite_network@graph@edgeL
    expect_equal(names(edges), names(metabolite_data))
    expect_true(all(purrr::map(edges, ~ length(.$edges)) > 0))

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

test_that("adjacency graph object is constructed", {
    # Based on adjacency matrix and partial correlation matrix
    adj_graph <- compute_adjacency_graph(metabolite_data,
                                    metabolite_network)

    expect_identical(class(adj_graph), "igraph")
})

test_that("network is constructed even with missingness", {
    metabolite_network_na <- simulated_data %>%
        .insert_random_missingness() %>%
        dplyr::select(dplyr::matches("metabolite")) %>%
        nc_estimate_network()

    expect_s4_class(metabolite_network_na, "pcAlgo")

    # For number of neighbours, etc
    edges <- metabolite_network_na@graph@edgeL
    expect_equal(names(edges), names(metabolite_data))
    expect_true(all(purrr::map(edges, ~ length(.$edges)) > 0))
})
