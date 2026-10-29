# Estimate the undirected graph of the metabolic data.

Uses the PC-algorithm and is mostly a wrapper around
[`pcalg::skeleton()`](https://rdrr.io/pkg/pcalg/man/skeleton.html).

## Usage

``` r
pc_estimate_undirected_graph(data, alpha = 0.01)
```

## Arguments

- data:

  Input numeric data that forms the basis of the underlying graph.

- alpha:

  Significance level threshold applied to each test to determine
  conditional dependence for if an edge exists.

## Value

A `pcAlgo` object that contains the DAG skeleton, aka undirected graph.

## Details

This function estimates the "skeleton of a DAG", meaning a graph without
arrowheads, aka an undirected graph. The default estimation method used
is the "PC-stable" method, which estimates the *order-independent*
skeleton of the DAG, meaning the order of the variables given does not
impact the results (older versions of the algorithm were
order-dependent). The method also assumes no latent variables.

An edge is determined by testing for conditional dependence between two
nodes based on the
[`pcalg::gaussCItest()`](https://rdrr.io/pkg/pcalg/man/condIndFisherZ.html).
Conditional *independence* exists when the nodes have zero partial
correlation determined from a p-value based hypothesis test against the
correlation matrix of the data from the nodes. The estimated edges
exists between the *start* and *end* nodes when the *start* and *end*
variables are conditionally dependent given the subset of remaining
variables.

## See also

The help documentation of
[`pcalg::skeleton()`](https://rdrr.io/pkg/pcalg/man/skeleton.html) has
more details.
