context("Standardize metabolic variables.")

metabolite_data <- simulated_data %>%
    select(starts_with("metabolite"))

standardized <- nc_standardize(metabolite_data)

standardized_with_residuals <- simulated_data %>%
    nc_standardize(starts_with("metabolite"), regressed_on = "age") %>%
    select(starts_with("metabolite"))

test_that("logging and scaling works", {
    metabolite_means <- standardized %>%
        colSums(na.rm = TRUE) %>%
        round(0)

    metabolite_sd <- standardized %>%
        purrr::map_dbl(sd, na.rm = TRUE) %>%
        round(0)

    expect_equal(sum(metabolite_means), 0)
    expect_true(all(metabolite_sd == 1))
    expect_false(identical(simulated_data, standardized))
})

test_that("standardization with residuals works", {
    metabolite_means_resid <- standardized_with_residuals %>%
        colSums(na.rm = TRUE) %>%
        round(0)

    metabolite_sd_resid <- standardized_with_residuals %>%
        purrr::map_dbl(sd, na.rm = TRUE) %>%
        round(0)

    expect_equal(sum(metabolite_means_resid), 0)
    expect_true(all(metabolite_sd_resid == 1))
    expect_false(identical(simulated_data, standardized_with_residuals))
    expect_false(identical(standardized, standardized_with_residuals))
})

# TODO: Add test for highly skewed data?
