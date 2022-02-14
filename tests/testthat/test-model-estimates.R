context("Exposure and outcome multiple model estimates.")

std_sim_data <- simulated_data %>%
    nc_standardize(starts_with("metabolite"))

metabolite_names <- simulated_data %>%
    select(starts_with("metabolite")) %>%
    names()

exposure_estimates <- std_sim_data %>%
    nc_estimate_exposure_links(
        edge_tbl = as_edge_tbl(metabolite_network),
        exposure = "exposure",
        model_function = lm
    )

outcome_estimates <- std_sim_data %>%
    nc_estimate_outcome_links(
        edge_tbl = as_edge_tbl(metabolite_network),
        outcome = "outcome_continuous",
        model_function = lm
    )

test_that("Correct metabolites are classified as direct", {
    exposure_direct_effect_vars <- exposure_estimates %>%
        dplyr::filter(effect == "direct") %>%
        pull(index_node)

    expected_exposure_associations <- paste0("metabolite_", c(1, 10, 8))
    # Expect at least the three associations to be identified.
    expect_equal(sum(exposure_direct_effect_vars %in% expected_exposure_associations), 3)

    outcome_direct_effect_vars <- outcome_estimates %>%
        dplyr::filter(effect == "direct") %>%
        pull(index_node)

    expected_outcome_associations <- paste0("metabolite_", c(3, 9, 12))
    # Expect at least the three associations to be identified.
    expect_equal(sum(outcome_direct_effect_vars %in% expected_outcome_associations), 3)
})

test_that("estimations outputs as correct class", {
    expect_correct_model_results(exposure_estimates, metabolite_names)
    expect_correct_model_results(outcome_estimates, metabolite_names)
})

test_that("one or more adjustment variables can be added", {
    one_adj <- std_sim_data %>%
        nc_estimate_exposure_links(
            edge_tbl = as_edge_tbl(metabolite_network),
            exposure = "exposure",
            adjustment_vars = "age",
            model_function = lm
        )

    expect_correct_model_results(one_adj, metabolite_names)

    two_adj <- std_sim_data %>%
        mutate(Random = rnorm(nrow(.))) %>%
        nc_estimate_exposure_links(
            edge_tbl = as_edge_tbl(metabolite_network),
            exposure = "exposure",
            adjustment_vars = c("age", "Random"),
            model_function = lm
        )

    expect_correct_model_results(two_adj, metabolite_names)

    two_adj_outcome <- std_sim_data %>%
        mutate(Random = rnorm(nrow(.))) %>%
        nc_estimate_outcome_links(
            edge_tbl = as_edge_tbl(metabolite_network),
            outcome = "outcome_continuous",
            adjustment_vars = c("age", "Random"),
            model_function = lm
        )

    expect_correct_model_results(two_adj, metabolite_names)
})

test_that("assertion checks pass", {
    # for model_function
    expect_error(std_sim_data %>%
        nc_estimate_exposure_links(
            edge_tbl = as_edge_tbl(metabolite_network),
            exposure = "exposure",
            model_function = "lm"
        )
    )

    # for exposure
    expect_error(std_sim_data %>%
        nc_estimate_exposure_links(
            edge_tbl = as_edge_tbl(metabolite_network),
            exposure = c("exposure", "SomethingElse"),
            model_function = lm
        )
    )

    # for graph
    expect_error(std_sim_data %>%
        nc_estimate_exposure_links(
            edge_tbl = swiss,
            exposure = "exposure",
            model_function = lm
        )
    )

    # for adjustment
    expect_error(std_sim_data %>%
        nc_estimate_exposure_links(
            edge_tbl = as_edge_tbl(metabolite_network),
            exposure = "exposure",
            adjustment_vars = 1,
            model_function = lm
        )
    )

    expect_error(std_sim_data %>%
        nc_estimate_exposure_links(
            edge_tbl = as_edge_tbl(metabolite_network),
            exposure = "exposure",
            adjustment_vars = swiss,
            model_function = lm
        )
    )
})

test_that("missingness in data still provides results", {
    set.seed(21451)
    missingness_data <- insert_random_missingness(std_sim_data)
    metabolite_network <- missingness_data %>%
        select(starts_with("metabolite")) %>%
        nc_estimate_network()
    exposure_estimates <- missingness_data %>%
        nc_estimate_exposure_links(
            edge_tbl = as_edge_tbl(metabolite_network),
            exposure = "exposure",
            model_function = lm
        )

    expect_correct_model_results(exposure_estimates, metabolite_names)
})

# test_that("computes when using survival::Surv and coxph", {
#     skip_on_ci()
#     skip_if_not_installed("survival")
#     outcome_estimates <- std_sim_data %>%
#         nc_estimate_outcome_links(
#             edge_tbl = as_edge_tbl(metabolite_network),
#             outcome = "survival::Surv(survival_time, case_status)",
#             model_function = survival::coxph
#         )
#
#     expect_correct_model_results(outcome_estimates, metabolite_names)
# })

test_that("Factor confounders are extracted properly", {
    multimodel_exposure <- std_sim_data %>%
        mutate(Sex = sample(rep(c("F", "M"), times = nrow(.) / 2))) %>%
        nc_estimate_exposure_links(
            edge_tbl = as_edge_tbl(metabolite_network),
            exposure = "exposure",
            adjustment_vars = c("age", "Sex"),
            model_function = lm
        )

    expect_correct_model_results(multimodel_exposure, metabolite_names)

    # Same expected metabolites
    exposure_direct_effect_vars <- multimodel_exposure %>%
        dplyr::filter(effect == "direct") %>%
        pull(index_node)

    expected_exposure_associations <- paste0("metabolite_", c(1, 10, 8))
    expect_identical(exposure_direct_effect_vars, expected_exposure_associations)
})
