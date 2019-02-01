
#' Rename feature names in order to avoid clash with glm#
#'
#' @param dat Dataset of metabolites (columns) by samples (rows)
#'
#' @return Outputs a list of renamed column variables.
#' @export
#'
rename.met <- function(.data,.nodes) {

  # dat: samples x metabolites data matrix

  Ll <- paste("NM",  c(1:length(.nodes)), sep = "") # generate shorter metabolite names

  names_mapping <- cbind(.nodes, Ll) # mapping of old and new metabolite names
  colnames(names_mapping) <- c("Metabolite", "Outcome")

  data_renamed <- .data
  colnames(data_renamed) <- Ll # is character!
  data_renamed<-dplyr::bind_cols(.data,data_renamed)
  return(list(data_renamed = data_renamed, names_mapping = names_mapping))
}


