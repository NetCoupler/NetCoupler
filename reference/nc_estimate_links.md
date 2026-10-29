# Compute model estimates between an external (exposure or outcome) variable and a network.

**\[experimental\]**

This is the main function that identifies potential links between
external factors and the network. There are two functions to estimate
and classify links:

- `nc_estimate_exposure_links()`: Computes the model estimates for the
  exposure side.

- `nc_estimate_outcome_links()`: Computes the model estimates for the
  exposure side.

## Usage

``` r
nc_estimate_exposure_links(
  data,
  edge_tbl,
  exposure,
  adjustment_vars = NA,
  model_function,
  model_arg_list = NULL,
  exponentiate = FALSE,
  classify_option_list = classify_options()
)

nc_estimate_outcome_links(
  data,
  edge_tbl,
  outcome,
  adjustment_vars = NA,
  model_function,
  model_arg_list = NULL,
  exponentiate = FALSE,
  classify_option_list = classify_options()
)
```

## Arguments

- data:

  The data.frame or tibble that contains the variables of interest,
  including the variables used to make the network.

- edge_tbl:

  Output graph object from
  [`nc_estimate_network()`](https://netcoupler.github.io/NetCoupler/reference/nc_estimate_network.md),
  converted to an edge table using
  [`as_edge_tbl()`](https://netcoupler.github.io/NetCoupler/reference/as_edge_tbl.md).

- exposure, outcome:

  Character. The exposure or outcome variable of interest.

- adjustment_vars:

  Optional. Variables to adjust for in the models.

- model_function:

  A function for the model to use (e.g.
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html),
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html),
  survival::coxph()). Can be any model as long as the function has the
  arguments `formula` and `data`. Type in the model function as a bare
  object (without `()`, for instance as `lm`).

- model_arg_list:

  Optional. A list containing the named arguments that will be passed to
  the model function. A simple example would be
  `list(family = binomial(link = "logit"))` to specify that the `glm`
  model is a logistic model and not a linear one. See the examples for
  more on the usage.

- exponentiate:

  Logical. Whether to exponentiate the log estimates, as computed with
  e.g. logistic regression models.

- classify_option_list:

  A list with classification options for direct, ambigious, or no
  effects. Used with the
  [`classify_options()`](https://netcoupler.github.io/NetCoupler/reference/classify_options.md)
  function with the arguments:

  - `single_metabolite_threshold`: Default of 0.05. P-values from models
    with only the index metabolite (no neighbour adjustment) are
    classified as effects if below this threshold. For larger sample
    sizes and networks, we recommend lowering the threshold to reduce
    risk of false positives.

  - `network_threshold`: Default of 0.1. P-values from any models that
    have direct neighbour adjustments are classified as effects if below
    this threshold. This is assumed as a one-sided p-value threshold.
    Like the threshold above, a lower value should be used for larger
    sample sizes and networks.

  - `direct_effect_adjustment`: Default is NA. After running the
    algorithm once, sometimes it's useful to adjust for the direct
    effects identified to confirm whether other links exist.

## Value

Outputs a
[tibble](https://tibble.tidyverse.org/reference/tibble-package.html)
that contains the model estimates from either the exposure or outcome
side of the network as well as the effect classification. Each row
represents the "no neighbour node adjusted" model and has the results
for the outcome/exposure to index node pathway. Columns for the outcome
are:

- `outcome` or `exposure`: The name of the variable used as the external
  variable.

- `index_node`: The name of the metabolite used as the index node from
  the network. In combination with the outcome/exposure variable, they
  represent the individual model used for the classification.

- `estimate`: The estimate from the outcome/exposure and index node
  model.

- `std_error`: The standard error from the outcome/exposure and index
  node model.

- `fdr_p_value`: The False Discovery Rate-adjusted p-value from the
  outcome/exposure and index node model.

- `effect`: The NetCoupler classified effect between the index node and
  the outcome/exposure. Effects are classified as "direct" (there is a
  probable link based on the given thresholds), "ambigious" (there is a
  potential link but not all thresholds were passed), and "none" (no
  potential link seen).

The tibble output also has an attribute that contains all the models
generated *before* classification. Access it with
`attr(output, "all_models_df")`.

## See also

`vignette("examples")` article has more details on how to use NetCoupler
with different models.

## Examples

``` r

standardized_data <- simulated_data %>%
    nc_standardize(starts_with("metabolite"))

metabolite_network <- simulated_data %>%
    nc_standardize(starts_with("metabolite"),
                   regressed_on = "age") %>%
    nc_estimate_network(starts_with("metabolite"))
edge_table <- as_edge_tbl(metabolite_network)

results <- standardized_data %>%
  nc_estimate_exposure_links(
    edge_tbl = edge_table,
    exposure = "exposure",
    model_function = lm
   )
results
#> # A tibble: 12 × 6
#>    exposure index_node    estimate std_error fdr_p_value effect   
#>  * <chr>    <chr>            <dbl>     <dbl>       <dbl> <chr>    
#>  1 exposure metabolite_1   0.173      0.0228      0      direct   
#>  2 exposure metabolite_10  0.318      0.0219      0      direct   
#>  3 exposure metabolite_11  0.0543     0.0232      0.0409 ambiguous
#>  4 exposure metabolite_12  0.0242     0.0231      0.380  none     
#>  5 exposure metabolite_2  -0.0430     0.0231      0.106  ambiguous
#>  6 exposure metabolite_3   0.0411     0.0231      0.123  ambiguous
#>  7 exposure metabolite_4   0.00344    0.0232      0.920  none     
#>  8 exposure metabolite_5   0.0479     0.0232      0.0717 ambiguous
#>  9 exposure metabolite_6  -0.0189     0.0230      0.506  none     
#> 10 exposure metabolite_7  -0.162      0.0229      0      direct   
#> 11 exposure metabolite_8  -0.355      0.0216      0      direct   
#> 12 exposure metabolite_9   0.0571     0.0230      0.0292 ambiguous

# Get results of all models used prior to classification

```
