context("Standardize metabolic variables.")

metabolite_data <- simulated_data %>%
    dplyr::select(dplyr::matches("metabolite"))

test_that("logging and scaling works", {
    std_data <- nc_standardize(metabolite_data)

    # expect_


})
