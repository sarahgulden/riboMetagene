#' Load saved 5' end count tables from a folder
#'
#' @param dir Folder with the <sample_name>.rds files written by extract_5prime()
#' @return A named list of count tables, named by sample_name
#' @export
read_5prime_tables <- function(dir) {
  files  <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  counts <- lapply(files, readRDS)
  names(counts) <- tools::file_path_sans_ext(basename(files))   # removes only the final ".rds"
  counts
}
