context("Exposure and outcome multiple model estimates.")

test_that("exposure side estimation outputs correctly", {
    exposure_estimates <- simulated_data %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = "exposure",
            .model_function = lm
        )

    expect_type(exposure_estimates, "list")
    expect_identical(class(exposure_estimates)[1], "tbl_df")
})

test_that("outcome side estimation outputs correctly", {
    outcome_estimates <- simulated_data %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "case_status",
            .model_function = glm,
            .model_arg_list = list(family = binomial(link = "logit")),
            .exponentiate = TRUE
        )

    expect_type(outcome_estimates, "list")
    expect_identical(class(outcome_estimates)[1], "tbl_df")
})

test_that("one or more adjustment variable can be added", {
    one_adj <- simulated_data %>%
        compute_model_estimates(
            .graph = metabolite_network,
            .external_var = "exposure",
            .external_side = "exposure",
            .adjustment_vars = "age",
            .model_function = lm
        )

    expect_type(one_adj, "list")
    expect_identical(class(one_adj)[1], "tbl_df")

    two_adj <- simulated_data %>%
        mutate(Random = rnorm(nrow(.))) %>%
        compute_model_estimates(
            .graph = metabolite_network,
            .external_var = "exposure",
            .external_side = "exposure",
            .adjustment_vars = c("age", "Random"),
            .model_function = lm
        )

    expect_type(two_adj, "list")
    expect_identical(class(two_adj)[1], "tbl_df")

    adjusted_variables <- unique(two_adj$adjusted_vars)
    adjusted_variables <- strsplit(adjusted_variables, ", ")[[1]]
    expect_identical(adjusted_variables, c("age", "Random"))
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
        nc_estimate_network()
    exposure_estimates <- missingness_data %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = "exposure",
            .model_function = lm
        )

    expect_type(exposure_estimates, "list")
    expect_identical(class(exposure_estimates)[1], "tbl_df")
})

test_that("computes when using survival::Surv and coxph", {
    skip_on_ci()
    skip_if_not_installed("survival")
    outcome_estimates <- simulated_data %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "survival::Surv(survival_time, case_status)",
            .model_function = survival::coxph
        )

    expect_type(outcome_estimates, "list")
    expect_identical(class(outcome_estimates)[1], "tbl_df")
})
