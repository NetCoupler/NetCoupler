
# Metabolic variables should be log-transformed and standardized on potential
# confounders {{this the best way?}} by regressing on age, sex, BMI, and prevalence
# of hypertension and using the residual variance in further models. Metabolites were
# then scaled (mean of zero and standard deviation of one). {{both regressed and
# scaled?}}

nc_standardize <- function(.tbl, .regressed_on = NULL) {
    if (!is.null(.regressed_on)) {
        assertive.types::assert_is_character(.regressed_on)
        .tbl %>%
            dplyr::mutate_all(.log_regress_standardize, .regressed_on)
    } else {
        .tbl %>%
            dplyr::mutate_all(.log_standardize)
    }
}

.log_standardize <- function(x) {
    scale(log(x))
}

.log_regress_standardize <- function(x, regressed_on) {
    # TODO: Decide which method to regress by. lm only?
    logged_x <- log(x)
    residual_x <- stats::residuals(stats::glm.fit(logged_x = regressed_on, y = x))
    scale(residual_x)
}
