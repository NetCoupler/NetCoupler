context("Run multiple models and generate model results.")

test_that("model estimation works and results are output", {
    metabolite_network <- simulated_data %>%
        select(matches("metabolite")) %>%
        nc_create_network()
    outcome_estimates <- simulated_data %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "survival::Surv(survival_time, case_status)",
            .adjustment_vars = "Age",
            .model_function = survival::coxph
        )

    expect_equal(ncol(outcome_estimates), 9)
    expect_type(outcome_estimates, "list")
    expect_identical(class(outcome_estimates)[1], "tbl_df")
})
