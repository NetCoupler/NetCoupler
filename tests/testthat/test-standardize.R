context("Standardize metabolic variables.")

test_that("logging and scaling works", {
    std_data <- simulated_data %>%
        nc_standardize(vars(matches("metabolite")))

    metabolite_means <- std_data %>%
        dplyr::select(matches("metabolite")) %>%
        colSums() %>%
        round(0)

    expect_equal(sum(metabolite_means), 0)

    std_data <- simulated_data %>%
        nc_standardize(vars(matches("metabolite")))

    metabolite_means <- std_data %>%
        dplyr::select(matches("metabolite")) %>%
        colSums() %>%
        round(0)

    expect_equal(sum(metabolite_means), 0)
})
