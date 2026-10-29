# Classification options for direct, ambigious, and no effect.

Classification options for direct, ambigious, and no effect.

## Usage

``` r
classify_options(
  single_metabolite_threshold = 0.05,
  network_threshold = 0.1,
  direct_effect_adjustment = NA
)
```

## Arguments

- single_metabolite_threshold, network_threshold,
  direct_effect_adjustment:

  See the `classify_option_list` argument in
  [nc_estimate_links](https://netcoupler.github.io/NetCoupler/reference/nc_estimate_links.md)
  for details.

## Value

List with options for the classification.
