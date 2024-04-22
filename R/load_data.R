#' Load RNA-seq and clinical data
#'
#' This function reads RNA-seq and clinical data from CSV files and assigns them to the global environment.
#'
#' @importFrom utils read.csv
#' @param RNAseq_path Path to the CSV file containing RNA-seq data.
#' @param Clinical_path Path to the CSV file containing clinical data.
#'
#' @return None
#'
#' @export
load_data <- function(RNAseq_path, Clinical_path) {
  rna_df <- read.csv(RNAseq_path, header = TRUE, check.names = FALSE, row.names = 1)
  clinical_df <- read.csv(Clinical_path, header = TRUE, row.names = 1)

  data <- list(rna_df = rna_df, clinical_df = clinical_df)
  return(data)
}
