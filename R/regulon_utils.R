#' Filter a regulon object by TF list, likelihood and correlation thresholds
#'
#' @param regulon A regulon object returned by corto()
#' @param tf_list Character vector with names of transcription factors to retain
#' @param likelihood_threshold Minimum likelihood (0 to 1) to keep an interaction
#' @param correlation_threshold Minimum absolute correlation to keep an interaction
#'
#' @return A filtered regulon object or NULL if empty
#' @export
filter_regulon <- function(regulon, 
                           tf_list = names(regulon), 
                           likelihood_threshold = 0, 
                           correlation_threshold = 0) {
  
  regulon_filtered <- list()
  
  for (tf in tf_list) {
    if (is.null(regulon[[tf]])) next
    
    targets <- names(regulon[[tf]]$tfmode)
    correlations <- regulon[[tf]]$tfmode
    likelihoods <- regulon[[tf]]$likelihood

    keep <- (likelihoods >= likelihood_threshold) & (abs(correlations) >= correlation_threshold)

    if (sum(keep) == 0) next

    regulon_filtered[[tf]] <- list(
      tfmode = correlations[keep],
      likelihood = likelihoods[keep]
    )
  }

  if (length(regulon_filtered) == 0) {
    message("Regulon is empty after applying filters.")
    return(NULL)
  }
  
  return(regulon_filtered)
}

#' Convert a regulon object into a tibble
#'
#' @param regulon A regulon object (filtered or full)
#' @param tf_list Character vector of TFs to include (default: all)
#' @param likelihood_threshold Minimum likelihood required for inclusion
#' @param export Logical. If TRUE, saves the result as a TSV file
#' @param filename Name of the output file if export = TRUE
#'
#' @return A tibble with columns: tf, target, correlation, and likelihood
#' @export
getregulon <- function(regulon, 
                       tf_list = names(regulon), 
                       likelihood_threshold = 0, 
                       export = FALSE, 
                       filename = "regulon.tsv") {
  
  tsv_list <- lapply(tf_list, function(name) {
    if (is.null(regulon[[name]])) {
      return(NULL)
    } else {
      tf <- rep(name, length(regulon[[name]]$tfmode))
      likelihoods <- regulon[[name]]$likelihood
      targets <- names(regulon[[name]]$tfmode)
      correlations <- regulon[[name]]$tfmode

      keep <- likelihoods >= likelihood_threshold
      if (sum(keep) == 0) return(NULL)

      tsv_temp <- tibble::tibble(
        tf = tf[keep],
        target = targets[keep],
        correlation = correlations[keep],
        likelihood = likelihoods[keep]
      )
      tsv_temp <- dplyr::arrange(tsv_temp, tf, dplyr::desc(likelihood))
      return(tsv_temp)
    }
  })
  
  tsv <- dplyr::bind_rows(tsv_list)
  
  if (export) {
    readr::write_tsv(tsv, filename)
    message("TSV file exported to: ", filename)
  }
  
  return(tsv)
}
