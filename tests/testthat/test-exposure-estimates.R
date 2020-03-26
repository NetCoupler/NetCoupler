context("Exposure-side multiple model computation.")

metabolite_network <- simulated_data %>%
    select(matches("metabolite")) %>%
    nc_create_network()

test_that("model estimation works and results are output", {
    exposure_estimates <- simulated_data %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = "exposure",
            .adjustment_vars = "age",
            .model_function = lm
        )

    expect_type(exposure_estimates, "list")
    expect_identical(class(exposure_estimates)[1], "tbl_df")
})

test_that("adjustment variable can be added", {
    exposure_estimates <- simulated_data %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = "exposure",
            .model_function = lm
        )

    expect_type(exposure_estimates, "list")
    expect_identical(class(exposure_estimates)[1], "tbl_df")
})

test_that("assertion checks pass", {
    # for model_function
    expect_error(simulated_data %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = "exposure",
            .model_function = "lm"
        )
    )

    # for exposure
    expect_error(simulated_data %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = c("exposure", "SomethingElse"),
            .model_function = lm
        )
    )

    # for graph
    expect_error(simulated_data %>%
        nc_exposure_estimates(
            .graph = swiss,
            .exposure = "exposure",
            .model_function = lm
        )
    )

    # for adjustment
    expect_error(simulated_data %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = "exposure",
            .adjustment_vars = 1,
            .model_function = lm
        )
    )

    expect_error(simulated_data %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = "exposure",
            .adjustment_vars = swiss,
            .model_function = lm
        )
    )
})

test_that("missingness in data still provides results", {
    missingness_data <- .insert_random_missingness(simulated_data)
    metabolite_network <- missingness_data %>%
        select(matches("metabolite")) %>%
        nc_create_network()
    exposure_estimates <- missingness_data %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = "exposure",
            .model_function = lm
        )

    expect_type(exposure_estimates, "list")
    expect_identical(class(exposure_estimates)[1], "tbl_df")
})

test_that("adding more than one adjustment variable includes it in the model", {
    two_adj_vars <- simulated_data %>%
        mutate(Random = rnorm(nrow(.))) %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = "exposure",
            .adjustment_vars = c("age", "Random"),
            .model_function = lm
        )

    expect_type(two_adj_vars, "list")
    expect_identical(class(two_adj_vars)[1], "tbl_df")

    adjusted_variables <- unique(two_adj_vars$adjusted_vars)
    adjusted_variables <- strsplit(adjusted_variables, ", ")[[1]]
    expect_identical(adjusted_variables, c("age", "Random"))
})
