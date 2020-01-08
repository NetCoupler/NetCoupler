context("Standardize metabolic variables.")

standardized <- simulated_data %>%
    nc_standardize(vars(matches("metabolite")))

standardized_with_residuals <- simulated_data %>%
    nc_standardize(vars(matches("metabolite")), "age")

test_that("logging and scaling works", {
    metabolite_means <- standardized %>%
        dplyr::select(matches("metabolite")) %>%
        colSums() %>%
        round(0)

    metabolite_sd <- standardized %>%
        dplyr::select(matches("metabolite")) %>%
        purrr::map_dbl(sd) %>%
        round(0)

    expect_equal(sum(metabolite_means), 0)
    expect_true(all(metabolite_sd == 1))
    expect_false(identical(simulated_data, standardized))
})

test_that("standardization with residuals works", {
    metabolite_means_resid <- standardized_with_residuals %>%
        dplyr::select(matches("metabolite")) %>%
        colSums() %>%
        round(0)

    metabolite_sd_resid <- standardized_with_residuals %>%
        dplyr::select(matches("metabolite")) %>%
        purrr::map_dbl(sd) %>%
        round(0)

    expect_equal(sum(metabolite_means_resid), 0)
    expect_true(all(metabolite_sd_resid == 1))
    expect_false(identical(simulated_data, standardized_with_residuals))
    expect_false(identical(standardized, standardized_with_residuals))
})