context("Outcome-side multiple model computation.")

test_that("model estimation works and results are output", {
    metabolite_network <- simulated_data %>%
        select(matches("metabolite")) %>%
        nc_create_network()
    outcome_estimates <- simulated_data %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "survival::Surv(survival_time, case_status)",
            .adjustment_vars = "age",
            .model_function = survival::coxph,
            .exponentiate = TRUE
        )

    expect_equal(ncol(outcome_estimates), 10)
    expect_type(outcome_estimates, "list")
    expect_identical(class(outcome_estimates)[1], "tbl_df")
})

test_that("assertion checks pass", {
    # for model_function
    expect_error(simulated_data %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "survival::Surv(survival_time, case_status)",
            .model_function = "survival::coxph"
        )
    )

    # for outcome
    expect_error(simulated_data %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = c("survival::Surv(survival_time, case_status)", "SomethingElse"),
            .model_function = survival::coxph
        )
    )

    # for graph
    expect_error(simulated_data %>%
        nc_outcome_estimates(
            .graph = swiss,
            .outcome = "survival::Surv(survival_time, case_status)",
            .model_function = survival::coxph
        )
    )

    # for adjustment
    expect_error(simulated_data %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "survival::Surv(survival_time, case_status)",
            .adjustment_vars = 1,
            .model_function = survival::coxph
        )
    )

    expect_error(simulated_data %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "survival::Surv(survival_time, case_status)",
            .adjustment_vars = swiss,
            .model_function = survival::coxph
        )
    )
})
