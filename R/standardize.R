
#' Standardize the metabolic variables.
#'
#' Can standardize by either 1) [log()]-transforming and then applying [scale()]
#' (mean-center and scaled by standard deviation), or 2) if `.regressed_on`
#' variables are given, then log-transforming, running a linear regression to obtain
#' the [stats::residuals()], and finally scaled. Use `.regressed_on` to try to
#' remove influence of potential confounding.
#'
#' @param .tbl Data frame.
#' @param .vars Metabolic variables that will make up the network.
#' @param .regressed_on Optional. A character vector of variables to regress the
#'   metabolic variables on. Use if you want to standardize the metabolic variables
#'   on variables that are known to influence them, e.g. sex or age. Calculates
#'   the residuals from a linear regression model.
#'
#' @return Outputs a [tibble][tibble::tibble-package] object, with the original metabolic
#'   variables now standardized.
#' @export
#'
#' @examples
#'
#' # Don't regress on any variable.
#' simulated_data %>%
#'   nc_standardize(vars(matches("metabolite_"))) %>%
#'   tibble::as_tibble()
#'
#' # Don't regress on any variable.
#' simulated_data %>%
#'   nc_standardize(vars(matches("metabolite_")), "age") %>%
#'   tibble::as_tibble()
nc_standardize <- function(.tbl, .vars, .regressed_on = NULL) {
    if (!is.null(.regressed_on)) {
        assertive.types::assert_is_character(.regressed_on)
        .tbl %>%
            dplyr::mutate_at(.vars, .funs = .log_regress_standardize,
                             regressed_on = .tbl[.regressed_on])
    } else {
        .tbl %>%
            dplyr::mutate_at(.vars, .funs = .log_standardize)
    }
}

.log_standardize <- function(x) {
    as.numeric(scale(log(x)))
}

.log_regress_standardize <- function(x, regressed_on) {
    # TODO: Decide which method to regress by. lm only?
    logged_x <- log(x)
    residual_x <- stats::residuals(stats::glm.fit(y = logged_x, x = regressed_on))
    as.numeric(scale(residual_x))
}
