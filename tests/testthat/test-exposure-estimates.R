context("Exposure-side multiple model computation.")

test_that("model estimation works and results are output", {
    metabolite_network <- simulated_data %>%
        select(matches("metabolite")) %>%
        nc_create_network()
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
