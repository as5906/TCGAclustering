#' Remove rows not shared between two data frames
#'
#' This function removes rows that are not shared between two data frames and returns a merged data frame.
#'
#' @param df1 First data frame.
#' @param df2 Second data frame.
#'
#' @return Merged data frame with shared rows.
#'
#' @export
remove_unshared_rows <- function(df1, df2) {
  merged_df <- merge(df1, df2, by = "row.names", all = TRUE)
  merged_df <- merged_df[complete.cases(merged_df), ]
  rownames(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[, -1]
  return(merged_df)
}
