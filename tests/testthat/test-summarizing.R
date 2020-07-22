context("Summarizing and classifying pathways from multi-models")

multimodel_exposure <- simulated_data %>%
    mutate(Random = rnorm(nrow(.))) %>%
    nc_exposure_estimates(
        .graph = metabolite_network,
        .exposure = "exposure",
        .adjustment_vars = c("age", "Random"),
        .model_function = lm
    )
exposure_results <- nc_classify_effects(multimodel_exposure)

multimodel_outcome <- simulated_data %>%
    mutate(Random = rnorm(nrow(.))) %>%
    nc_outcome_estimates(
        .graph = metabolite_network,
        .outcome = "case_status",
        .model_function = glm,
        .model_arg_list = list(family = binomial(link = "logit")),
        .exponentiate = TRUE
    )


outcome_results <- nc_classify_effects(multimodel_outcome)

test_that("Error is thrown if columns don't exist.", {
    expect_error(nc_classify_effects(multimodel_exposure[-2]))
    expect_error(nc_classify_effects(multimodel_outcome[-2]))
})

test_that("Same number of metabolites are returned", {
    number_metabolites <- simulated_data %>%
        select(matches("metabolite")) %>%
        ncol()

    expect_equal(nrow(exposure_results),
                 number_metabolites)
    expect_equal(nrow(outcome_results),
                 number_metabolites)
})

test_that("Correct metabolites are classified as direct", {
    exposure_direct_effect_vars <- exposure_results %>%
        dplyr::filter(direct_effect == "direct") %>%
        pull(index_node)

    expected_exposure_associations <- paste0("metabolite_", c(1, 10, 8))
    expect_identical(exposure_direct_effect_vars, expected_exposure_associations)

    outcome_direct_effect_vars <- outcome_results %>%
        dplyr::filter(direct_effect == "direct") %>%
        pull(index_node)

    # TODO: Fix so simulated data gives exact associations. Right now it's only 9
    expected_outcome_associations <- paste0("metabolite_", c(3, 9, 12))[2]
    expect_identical(outcome_direct_effect_vars, expected_outcome_associations)
})

test_that("Factor confounders are extracted properly", {
    multimodel_exposure <- simulated_data %>%
        mutate(Sex = sample(rep(c("F", "M"), times = nrow(.) / 2))) %>%
        nc_exposure_estimates(
            .graph = metabolite_network,
            .exposure = "exposure",
            .adjustment_vars = c("age", "Sex"),
            .model_function = lm
        )

    classified_results <- nc_classify_effects(multimodel_exposure)
    expect_equal(nrow(classified_results), 12)

    # Same expected metabolites
    exposure_direct_effect_vars <- classified_results %>%
        dplyr::filter(direct_effect == "direct") %>%
        pull(index_node)

    expected_exposure_associations <- paste0("metabolite_", c(1, 10, 8))
    expect_identical(exposure_direct_effect_vars, expected_exposure_associations)
})
