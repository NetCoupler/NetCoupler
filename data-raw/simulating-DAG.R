library(dagitty)
library(dplyr)
library(survival)

# Create DAG to base simulation -------------------------------------------

# including exposure and network-variables (a-p)

# TODO: Update this so it is using lavaan instead?
dag_graph <- dagitty('dag {
    exposure -> metabolite_1
    exposure -> metabolite_8
    exposure -> metabolite_10
    metabolite_1 -> metabolite_2
    metabolite_2 -> metabolite_3
    metabolite_3 -> metabolite_4
    metabolite_4 -> metabolite_5
    metabolite_5 -> metabolite_6
    metabolite_2 -> metabolite_7
    metabolite_7 -> metabolite_8
    metabolite_8 -> metabolite_9
    metabolite_9 -> metabolite_10
    metabolite_10 -> metabolite_11
    metabolite_1 -> metabolite_12
    metabolite_3 -> metabolite_12
    metabolite_4 -> metabolite_12
    metabolite_5 -> metabolite_9
    metabolite_6 -> metabolite_10
    metabolite_10 -> metabolite_12
    metabolite_3 -> outcome_continuous
    metabolite_9 -> outcome_continuous
    metabolite_12 -> outcome_continuous
    exposure [exposure]
    outcome_continuous [outcome]
}')

simulated_dag_data <- simulateSEM(dag_graph, N = 2000) %>%
    mutate(across(matches("metabolite_"), ~ . + 6)) %>%
    as_tibble()

#' Survival time simulation.
#'
#' Based on the Gompertz distribution. Bender et al. (2003)  http://epub.ub.uni-muenchen.de/
#'
#' @param .data The simulated dataset based on the DAG.
#' @param IV1 Influential variable 1. Used in calculation of survival.
#' @param IV2 Influential variable 2.
#' @param IV3 Influential variable 3.
#'
#' @return
simulate_survival_time <- function(.data, IV1, IV2, IV3) {

    # Variation of these parameters will modify the distribution of survival times
    # (and thus the prevalence of the event at a given censoring date)

    gompertz_distribution_params <- data.frame(
        # TODO: Confirm it is width
        # Between 0 and 1
        width = round(runif(1, 0, 1), 2),
        alpha = 0.7,
        # TODO: Confirm lambda event meaning.
        # Steepness of slope
        lambda_event = 0.7 * 9 ** (-8)
    )

    X1 <- .data[[IV1]]
    X2 <- .data[[IV2]]
    X3 <- .data[[IV3]]

    number_observations <- nrow(.data)
    randomness <- rnorm(number_observations, mean = 10, sd = 1)

    exposure_effect_estimates <- exp(0.5 * X1 + 0.5 * X2 + (-0.5) * X3 + randomness)

    survival_time <-
        with(gompertz_distribution_params,
             (1 / alpha) * log(1 - ((alpha * log(width)) /
                                        (lambda_event * exposure_effect_estimates)))
        )

    sorted_survival_time <- sort(survival_time)

    # Simulate case status at specific censoring time
    incidence_cases <- number_observations / 10
    tenth_percentile <- (sorted_survival_time[incidence_cases] +
                             sorted_survival_time[incidence_cases + 1]) / 2

    case_status <- rep.int(0, times = number_observations)

    # TODO: Redo to not use for loop
    for (obs in 1:number_observations) {
        if (survival_time[obs] < tenth_percentile) {
            case_status[obs] <- 1
        }
        else {
            survival_time[obs] <- tenth_percentile
        }
    }

    data.frame(survival_time, case_status)
}

simulated_data <-
    bind_cols(
        simulate_survival_time(
            .data = simulated_dag_data,
            IV1 = "metabolite_3" ,
            IV2 = "metabolite_9",
            IV3 = "metabolite_12"
        ),
        simulated_dag_data
    ) %>%
    mutate(age = rnorm(n(), mean = 50, sd = 10)) %>%
    as_tibble()

usethis::use_data(simulated_data, overwrite = TRUE)
