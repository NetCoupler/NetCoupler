library(dplyr)
library(lavaan)
library(lavaanPlot)
library(simsurv)

# Create DAG to base simulation -------------------------------------------

dag_model <- "
# Network connections
metabolite_1 ~ metabolite_2 + metabolite_12
metabolite_2 ~ metabolite_3 + metabolite_7
metabolite_3 ~ metabolite_4 + metabolite_12
metabolite_4 ~ metabolite_5 + metabolite_12
metabolite_5 ~ metabolite_6 + metabolite_9
metabolite_6 ~ metabolite_10
metabolite_7 ~ metabolite_8
metabolite_8 ~ metabolite_9
metabolite_9 ~ metabolite_10
metabolite_10 ~ metabolite_11 + metabolite_12

# Exposure connections
metabolite_1 + metabolite_8 + metabolite_10 ~ exposure

# Outcome connections
outcome_continuous ~ metabolite_3 + metabolite_9 + metabolite_12

# Covariance (confounders)
exposure + outcome_continuous ~~ age
"

sim_data <- as_tibble(simulateData(dag_model, sample.nobs = 2000)) %>%
    mutate(id = 1:n())

# Check if DAG was created correctly
# sim_model_fit <- sem(model = dag_model, data = sim_data)
# lavaanPlot(model = sim_model_fit)


# Simulate survival times -------------------------------------------------

sim_event_time_data <-
    simsurv(
        dist = "gompertz",
        # "Scale" of distribution. Meaning how much it distributes from the left or right.
        # Combine with maxt to create right skewed distributions.
        # - Higher number gives a narrower numerical distribution.
        # - Lower number gives a wider numerical distribution.
        lambdas = 0.01,
        # "Shape" of distribution. Meaning where the peak occurs and how much it spreads.
        # - Higher number gives a more uniform distribution on left and right side,
        # with narrower spread.
        # - Lower number gives a more left skewed distribution with wider spread.
        gammas = 0.5,
        x = sim_data,
        # Where to "censor" time, aka when the study ended.
        maxt = 5,
        betas = c(
            metabolite_3 = 0.5,
            metabolite_9 = -0.5,
            metabolite_12 = 0.75,
            age = 0.3
        )
    ) %>%
    as_tibble() %>%
    rename(outcome_binary = status, outcome_event_time = eventtime)

simulated_data <-
    full_join(sim_data,
              sim_event_time_data) %>%
    mutate(age = age + 50)

# Save dataset ------------------------------------------------------------

usethis::use_data(simulated_data, overwrite = TRUE)
