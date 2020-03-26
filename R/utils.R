
.insert_random_missingness <- function(.tbl) {
    purrr::map_df(.tbl, ~ .[sample(
        c(TRUE, NA),
        prob = c(0.98, 0.02),
        size = length(.),
        replace = TRUE
    )])
}
