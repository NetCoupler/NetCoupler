# Create an estimate of the metabolic network as an undirected graph.

**\[experimental\]**

The main NetCoupler network creator. Uses the input data to estimate the
underlying undirected graph. The default uses the PC algorithm,
implemented within NetCoupler with
[`pc_estimate_undirected_graph()`](https://netcoupler.github.io/NetCoupler/reference/pc_estimate_undirected_graph.md)
Defaults to using the PC algorithm to calculate possible edges. Any
missing values in the input data are removed by this function, since
some computations can't handle missingness.

## Usage

``` r
nc_estimate_network(data, cols = everything(), alpha = 0.01)
```

## Arguments

- data:

  Data that would form the underlying network.

- cols:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Variables to include by using
  [`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
  style selection.

- alpha:

  The alpha level to use to test whether an edge exists or not. Default
  is 0.01.

## Value

Outputs a
[`tidygraph::tbl_graph()`](https://tidygraph.data-imaginist.com/reference/tbl_graph.html)
with the start and end nodes, as well as the edge weights.

## See also

See
[nc_estimate_links](https://netcoupler.github.io/NetCoupler/reference/nc_estimate_links.md)
for examples on using NetCoupler and
[pc_estimate_undirected_graph](https://netcoupler.github.io/NetCoupler/reference/pc_estimate_undirected_graph.md)
for more details on the PC-algorithm network estimation method.
