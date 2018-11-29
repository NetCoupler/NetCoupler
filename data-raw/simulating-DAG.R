library(dagitty)
library(dplyr)
library(survival)

# Create DAG to base simulation -------------------------------------------

# including exposure and network-variables (a-p)

dag_graph <- dagitty('dag{
    exposure -> b [beta=.15]
    exposure -> l [beta=.15]
    exposure -> o [beta=.15]
    a -> b [beta=.3]
    b -> c [beta=.2]
    c -> d [beta=.5]
    d -> e [beta=.3]
    e -> f [beta=.6]
    b -> k [beta=.3]
    k -> l [beta=.3]
    l -> m [beta=.5]
    m -> o [beta=.3]
    o -> p [beta=.5]
    a -> v [beta=.2]
    c -> v [beta=.5]
    d -> v [beta=.3]
    e -> m [beta=.5]
    k -> o [beta=.3]
    o -> p [beta=.65]
}')

simulated_dag_data <- simulateSEM(dag_graph, N = 20000) %>%
    setNames(c("exposure", paste0("metabolite_", 1:(length(.) - 1))))

#' Survival time simulation.
#'
#' Based on the Gompertz distribution.
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
        alpha = 0.5,
        # TODO: Confirm lambda event meaning.
        # Steepness of slope
        lambda_event = 0.7 * 10 ** (-11)
    )

    X1 <- .data[, IV1]
    X2 <- .data[, IV2]
    X3 <- .data[, IV3]

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
    )

usethis::use_data(simulated_data, overwrite = TRUE)
