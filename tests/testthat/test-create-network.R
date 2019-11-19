context("Create metabolic variable network")

test_that("network is created", {
    # Renaming necessary for right now.
    renaming <- rename_met(simulated_data)
    renamed_simulated_data <- renaming[[1]]
    matching_table_names_newnames <- renaming[[2]]
    nodes_short_names <- renamed_simulated_data %>%
        dplyr::select(contains("NM")) %>%
        names()

    # Make partial independence network from metabolite data
    metabolite_network <-
        nc_make_network(renamed_simulated_data, .05, nodes_short_names)

    expect_type(metabolite_network, "list")
})
