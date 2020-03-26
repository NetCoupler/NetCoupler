context("Outcome-side multiple model computation.")

metabolite_network <- simulated_data %>%
    select(matches("metabolite")) %>%
    nc_create_network()

test_that("model estimation works and results are output", {
    outcome_estimates <- simulated_data %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "survival::Surv(survival_time, case_status)",
            .adjustment_vars = "age",
            .model_function = survival::coxph,
            .exponentiate = TRUE
        )

    expect_type(outcome_estimates, "list")
    expect_identical(class(outcome_estimates)[1], "tbl_df")
})

test_that("adjustment variable can be excluded", {
    outcome_estimates <- simulated_data %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "survival::Surv(survival_time, case_status)",
            .model_function = survival::coxph
        )

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

test_that("missingness in data still provides results", {
    missingness_data <- .insert_random_missingness(simulated_data)
    metabolite_network <- missingness_data %>%
        select(matches("metabolite")) %>%
        nc_create_network()
    outcome_estimates <- missingness_data %>%
        # Since netcoupler can't handle omitting na from Surv object.
        filter_at(vars(survival_time, case_status), ~ !is.na(.)) %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "survival::Surv(survival_time, case_status)",
            .model_function = survival::coxph,
            .exponentiate = TRUE
        )

    expect_type(outcome_estimates, "list")
    expect_identical(class(outcome_estimates)[1], "tbl_df")
})

test_that("adding more than one adjustment variable includes it in the model", {
    two_adj_vars <- simulated_data %>%
        mutate(Random = rnorm(nrow(.))) %>%
        nc_outcome_estimates(
            .graph = metabolite_network,
            .outcome = "survival::Surv(survival_time, case_status)",
            .adjustment_vars = c("age", "Random"),
            .model_function = survival::coxph
        )

    expect_type(two_adj_vars, "list")
    expect_identical(class(two_adj_vars)[1], "tbl_df")

    adjusted_variables <- unique(two_adj_vars$adjusted_vars)
    adjusted_variables <- strsplit(adjusted_variables, ", ")[[1]]
    expect_identical(adjusted_variables, c("age", "Random"))
})
