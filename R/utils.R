
# Taken from forcats
lvls_reorder <- function (f, idx, ordered = NA) {
    if (!is.numeric(idx)) {
        stop("`idx` must be numeric", call. = FALSE)
    }
    if (!setequal(idx, lvls_seq(f)) || length(idx) != nlevels(f)) {
        stop("`idx` must contain one integer for each level of `f`",
             call. = FALSE)
    }
    refactor(f, levels(f)[idx], ordered = ordered)
}

# Taken from forcats
refactor <- function (f, new_levels, ordered = NA) {
    if (is.na(ordered)) {
        ordered <- is.ordered(f)
    }
    new_f <- factor(f, levels = new_levels, exclude = NULL, ordered = ordered)
    attributes(new_f) <- utils::modifyList(attributes(f), attributes(new_f))
    new_f
}

# Taken from forcats
fct_rev <- function(f) {
    f <- factor(f)
    lvls_reorder(f, rev(lvls_seq(f)))
}

# Taken from forcats
lvls_seq <- function(f) {
    seq_along(levels(f))
}
