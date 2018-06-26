
# 1.4.1 Round numeric columns in dataframes
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame
  # digits: number of digits to round
  numeric_columns <- sapply(x, class) == "numeric"
  x[numeric_columns] <- round(x[numeric_columns], digits)
  x
}

# 1.4.2 get object-names as string-variables
name.as.string <- function(v1) {
  deparse(substitute(v1))
}

# Get rownames as variable in dplyr#

draw.rownames.out <- function(dat) {
  dat %>% do(mutate(., Metabolite = as.factor(rownames(.))))
}

draw.rownames <- draw.rownames.out

draw_rownames_Exp <- function(.data) .data %>% do(mutate(., Exp = rownames(.)))
