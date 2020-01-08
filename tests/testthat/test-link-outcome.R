context("Networks links with the outcome")

# TODO: Deprecate this in 0.0.4
test_that("DAG network estimates with the outcome", {
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
            dplyr::select(age),
        survival_obj = "survival_object"
    )))

    expect_type(net_coupler_case, "list")
})
