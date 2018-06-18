
#' Rename feature names in order to avoid clash with glm#
#'
#' @param dat Dataset of metabolites (columns) by samples (rows)
#'
#' @return Outputs a list of renamed column variables.
#' @export
#'
rename.met <- function(dat) {

  # dat: samples x metabolites data matrix

  Ll <- paste("NM", c(1:dim(dat)[2]), sep = "") # generate shorter metabolite names

  names_mapping <- cbind(colnames(dat), Ll) # mapping of old and new metabolite names
  colnames(names_mapping) <- c("Metabolite", "Outcome")

  data_renamed <- dat
  colnames(data_renamed) <- Ll # is character!

  return(list(data_renamed = data_renamed, names_mapping = names_mapping))
}
