% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize.R
\name{nc_standardize}
\alias{nc_standardize}
\title{Standardize the metabolic variables.}
\usage{
nc_standardize(data, cols = everything(), regressed_on = NULL)
}
\arguments{
\item{data}{Data frame.}

\item{cols}{Metabolic variables that will make up the network.}

\item{regressed_on}{Optional. A character vector of variables to regress the
metabolic variables on. Use if you want to standardize the metabolic variables
on variables that are known to influence them, e.g. sex or age. Calculates
the residuals from a linear regression model.}
}
\value{
Outputs a \link[tibble:tibble-package]{tibble} object, with the original metabolic
variables now standardized.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

Can standardize by either 1) \code{\link[=log]{log()}}-transforming and then applying \code{\link[=scale]{scale()}}
(mean-center and scaled by standard deviation), or 2) if \code{regressed_on}
variables are given, then log-transforming, running a linear regression to obtain
the \code{\link[stats:residuals]{stats::residuals()}}, and finally scaled. Use \code{regressed_on} to try to
remove influence of potential confounding.
}
\examples{

# Don't regress on any variable
simulated_data \%>\%
  nc_standardize(starts_with("metabolite_"))

# Extract residuals by regressing on a variable
simulated_data \%>\%
  nc_standardize(starts_with("metabolite_"), "age")

# Works with factors too
simulated_data \%>\%
  dplyr::mutate(Sex = as.factor(sample(rep(c("F", "M"), times = nrow(.) / 2)))) \%>\%
  nc_standardize(starts_with("metabolite_"), c("age", "Sex"))

}
\seealso{
\link{nc_estimate_links} for more detailed examples or the \code{vignette("NetCoupler")}.
}
