context("Networks links with the outcome")

test_that("DAG network estimates with the outcome", {
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

    expect_type(net_coupler_case, "list")
})
