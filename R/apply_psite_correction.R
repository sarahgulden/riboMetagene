#' Apply P-site correction to a signal dataframe
#'
#' Shifts each readlength's positions so that its P-site offset
#' (from find_psite_offsets()) becomes the new position 0.
#'
#' @param signal_df Long dataframe from build_signal_df() (needs readlength, position)
#' @param offsets Dataframe from find_psite_offsets() (readlength, peak_position)
#' @return signal_df with an added position_corrected column
#' @export
apply_psite_correction <- function(signal_df, offsets) {

  corrected <- merge(signal_df, offsets, by = "readlength", all.x = TRUE)

  corrected$position_corrected <- corrected$position - corrected$peak_position

  return(corrected)
}
