% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{pc_dag_estimates}
\alias{pc_dag_estimates}
\title{Estimate equivalence class of DAG from the PC algorithm.}
\usage{
pc_dag_estimates(data, alpha = 0.01)
}
\arguments{
\item{data}{Input data, samples by metabolite matrix or as data.frame.}

\item{alpha}{Significance level threshold applied to each test.}
}
\value{
Outputs a \code{pcAlgo} object.
}
\description{
Is mostly a wrapper around \code{\link[pcalg:pc]{pcalg::pc()}}. Estimates an order-independent
skeleton.
}
\keyword{internal}
