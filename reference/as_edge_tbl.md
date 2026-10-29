# Convert network graphs to edge tables as tibbles/data.frames.

**\[experimental\]**

## Usage

``` r
as_edge_tbl(network_object)
```

## Arguments

- network_object:

  Network graph from
  [`nc_estimate_network()`](https://netcoupler.github.io/NetCoupler/reference/nc_estimate_network.md).

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble-package.html),
with at least two columns:

- `source_node`: The starting node (variable).

- `target_node`: The ending node (variable) that links with the source
  node.

- `adjacency_weight`: (Optional) The "weight" given to the edge, which
  represents the strength of the link between two nodes.

## See also

See
[nc_estimate_links](https://netcoupler.github.io/NetCoupler/reference/nc_estimate_links.md)
for examples on using NetCoupler.
