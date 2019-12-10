context("Extract estimates from model")

# TODO: Deprecate this in 0.0.4
test_that("estimates are correctly extracted from outcome side", {
    skip_if_not_installed("glmulti")
    library(glmulti)

    metabolite_network <- simulated_data %>%
        dplyr::select(dplyr::matches("metabolite")) %>%
        nc_create_network()

    survival_object <<-
        survival::Surv(simulated_data$survival_time, simulated_data$case_status)
    net_coupler_case <- suppressMessages(suppressWarnings(NetCoupler::net_coupler_out(
        graph_skel = metabolite_network,
        dat = simulated_data %>%
            dplyr::select(dplyr::contains("metabolite"), case_status),
        DE = NULL,
        adjustment_data = simulated_data %>%
            dplyr::select(Age),
        survival_obj = "survival_object"
    )))

    extracted_estimates <- suppressWarnings(getExp.coef.out(object = net_coupler_case,
                    metabolite = simulated_data %>%
                        dplyr::select(dplyr::contains("metabolite")) %>%
                        names())) %>%
        tibble::as_tibble()

    expect_identical(class(extracted_estimates)[1], "tbl_df")
})
