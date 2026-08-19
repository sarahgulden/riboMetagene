#' List BAM files and assign sample names
#'
#' Scans a folder for BAM files and extracts a sample name from each
#' filename (filename without the .bam extension). Any experiment-specific
#' metadata (tissue, timepoint, etc.) should be parsed separately from
#' the sample_name.
#'
#' @param base_dir Path to the folder containing BAM files
#' @param pattern File pattern to match (default: ".bam$")
#' @return A data.frame with columns: filepath, filename, sample_name
#' @export
list_bam_samples <- function(base_dir, pattern = "\\.bam$") {

  all_files <- list.files(base_dir, pattern = pattern, recursive = TRUE)

  filename <- basename(all_files)

  sample_name <- tools::file_path_sans_ext(filename)

  data.frame(
    filepath    = file.path(base_dir, all_files),
    filename    = filename,
    sample_name = sample_name,
    stringsAsFactors = FALSE
  )
}
