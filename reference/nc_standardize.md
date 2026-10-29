# Standardize the metabolic variables.

**\[experimental\]**

Can standardize by either 1)
[`log()`](https://rdrr.io/r/base/Log.html)-transforming and then
applying [`scale()`](https://rdrr.io/r/base/scale.html) (mean-center and
scaled by standard deviation), or 2) if `regressed_on` variables are
given, then log-transforming, running a linear regression to obtain the
[`stats::residuals()`](https://rdrr.io/r/stats/residuals.html), and
finally scaled. Use `regressed_on` to try to remove influence of
potential confounding.

## Usage

``` r
nc_standardize(data, cols = everything(), regressed_on = NULL)
```

## Arguments

- data:

  Data frame.

- cols:

  Metabolic variables that will make up the network.

- regressed_on:

  Optional. A character vector of variables to regress the metabolic
  variables on. Use if you want to standardize the metabolic variables
  on variables that are known to influence them, e.g. sex or age.
  Calculates the residuals from a linear regression model.

## Value

Outputs a
[tibble](https://tibble.tidyverse.org/reference/tibble-package.html)
object, with the original metabolic variables now standardized.

## See also

[nc_estimate_links](https://netcoupler.github.io/NetCoupler/reference/nc_estimate_links.md)
for more detailed examples or the
[`vignette("NetCoupler")`](https://netcoupler.github.io/NetCoupler/articles/NetCoupler.md).

## Examples

``` r

# Don't regress on any variable
simulated_data %>%
  nc_standardize(starts_with("metabolite_"))
#> # A tibble: 2,000 × 18
#>    metabolite_1 metabolite_2 metabolite_3 metabolite_4 metabolite_5 metabolite_6
#>           <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>
#>  1        0.318        0.105       1.30         NA           -0.181      -0.441 
#>  2       -1.06         0.622      -1.04         -0.113       -2.14       -0.774 
#>  3       -0.467        1.04        1.13          0.888        0.910      -0.280 
#>  4        1.16         0.162       0.115        -0.664        0.730      -0.202 
#>  5        0.336       -2.80        1.19          1.66         1.02        1.25  
#>  6        0.114        0.319      -0.0494       -1.34        -0.192       0.442 
#>  7        0.586        0.232       0.533         1.84         0.822      -1.58  
#>  8        0.959       -0.763      -0.375         0.523        1.29        1.24  
#>  9       -0.589       -1.01       -0.241        -1.28        -0.789      -0.0464
#> 10       -0.516       -0.575       0.484        -0.369       -0.432       1.48  
#> # ℹ 1,990 more rows
#> # ℹ 12 more variables: metabolite_7 <dbl>, metabolite_8 <dbl>,
#> #   metabolite_9 <dbl>, metabolite_10 <dbl>, outcome_continuous <dbl>,
#> #   metabolite_12 <dbl>, metabolite_11 <dbl>, exposure <dbl>, age <dbl>,
#> #   id <int>, outcome_event_time <dbl>, outcome_binary <int>

# Extract residuals by regressing on a variable
simulated_data %>%
  nc_standardize(starts_with("metabolite_"), "age")
#> # A tibble: 2,000 × 18
#>    metabolite_1 metabolite_2 metabolite_3 metabolite_4 metabolite_5 metabolite_6
#>           <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>
#>  1        0.296        0.121       1.28        NA            -0.198      -0.457 
#>  2       -1.07         0.626      -1.04        -0.0980       -2.16       -0.781 
#>  3       -0.483        1.05        1.11         0.920         0.902      -0.289 
#>  4        1.18         0.147       0.137       -0.714         0.745      -0.191 
#>  5        0.339       -2.80        1.19         1.66          1.02        1.25  
#>  6        0.108        0.322      -0.0530      -1.33         -0.197       0.440 
#>  7        0.584        0.234       0.532        1.84          0.822      -1.59  
#>  8        0.971       -0.770      -0.362        0.497         1.30        1.25  
#>  9       -0.616       -0.994      -0.264       -1.22         -0.809      -0.0623
#> 10       -0.528       -0.568       0.476       -0.346        -0.440       1.48  
#> # ℹ 1,990 more rows
#> # ℹ 12 more variables: metabolite_7 <dbl>, metabolite_8 <dbl>,
#> #   metabolite_9 <dbl>, metabolite_10 <dbl>, outcome_continuous <dbl>,
#> #   metabolite_12 <dbl>, metabolite_11 <dbl>, exposure <dbl>, age <dbl>,
#> #   id <int>, outcome_event_time <dbl>, outcome_binary <int>

# Works with factors too
simulated_data %>%
  dplyr::mutate(Sex = as.factor(sample(rep(c("F", "M"), times = nrow(.) / 2)))) %>%
  nc_standardize(starts_with("metabolite_"), c("age", "Sex"))
#> # A tibble: 2,000 × 19
#>    metabolite_1 metabolite_2 metabolite_3 metabolite_4 metabolite_5 metabolite_6
#>           <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>
#>  1        0.231       0.0720       1.25         NA           -0.225      -0.475 
#>  2       -1.13        0.579       -1.07         -0.104       -2.18       -0.799 
#>  3       -0.422       1.10         1.14          0.926        0.928      -0.272 
#>  4        1.12        0.102        0.112        -0.719        0.720      -0.207 
#>  5        0.276      -2.85         1.17          1.65         0.997       1.24  
#>  6        0.171       0.369       -0.0268       -1.32        -0.171       0.457 
#>  7        0.522       0.187        0.506         1.84         0.796      -1.61  
#>  8        0.911      -0.816       -0.387         0.491        1.28        1.23  
#>  9       -0.683      -1.04        -0.292        -1.23        -0.837      -0.0806
#> 10       -0.466      -0.522        0.502        -0.340       -0.415       1.50  
#> # ℹ 1,990 more rows
#> # ℹ 13 more variables: metabolite_7 <dbl>, metabolite_8 <dbl>,
#> #   metabolite_9 <dbl>, metabolite_10 <dbl>, outcome_continuous <dbl>,
#> #   metabolite_12 <dbl>, metabolite_11 <dbl>, exposure <dbl>, age <dbl>,
#> #   id <int>, outcome_event_time <dbl>, outcome_binary <int>, Sex <fct>
```
