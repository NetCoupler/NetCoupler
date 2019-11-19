context("Extract estimates from model")

test_that("estimates are correctly extracted from outcome side", {
    renaming <- rename_met(simulated_data)
    renamed_simulated_data <- renaming[[1]]
    matching_table_names_newnames <- renaming[[2]]
    nodes_short_names <- renamed_simulated_data %>%
        dplyr::select(contains("NM")) %>%
        names()

    metabolite_network <-
        nc_make_network(renamed_simulated_data, .05, nodes_short_names)
    survival_object <<-
        survival::Surv(simulated_data$survival_time, simulated_data$case_status)
    net_coupler_case <- suppressMessages(suppressWarnings(NetCoupler::net_coupler_out(
        graph_skel = metabolite_network$skel_est,
        dat = simulated_data %>%
            dplyr::select(contains("metabolite"), case_status),
        DE = NULL,
        adjustment_data = simulated_data %>%
            dplyr::select(Age),
        survival_obj = "survival_object"
    )))

    extracted_estimates <- suppressWarnings(getExp.coef.out(object = net_coupler_case,
                    metabolite = simulated_data %>%
                        dplyr::select(contains("metabolite")) %>% names())) %>%
        tibble::as_tibble()

    expect_identical(class(extracted_estimates)[1], "tbl_df")
})
