# Simulated dataset with an underlying Directed Graph structure for the metabolites.

Simulated dataset with an underlying Directed Graph structure for the
metabolites.

## Usage

``` r
simulated_data
```

## Format

The simulated dataset is a
[tibble](https://tibble.tidyverse.org/reference/tibble-package.html)
with the following variables:

- Two outcome variables (`outcome_continuous` and `outcome_binary`)
  along with survival time (`outcome_event_time`) that is used for the
  `outcome_binary` variable

- A generic `exposure` variable as continuous

- 12 `metabolite_*` variables

- An `age` variable used as a confounder
