# TRUE if the BAM ends with the standard BGZF end-of-file marker (= file is complete)
bam_is_complete <- function(bam_file) {

  eof <- as.raw(c(0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
                  0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
                  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00))

  size <- file.size(bam_file)
  if (is.na(size) || size < 28) return(FALSE)

  con <- file(bam_file, "rb")
  on.exit(close(con))
  seek(con, where = size - 28)          # jump to the last 28 bytes
  identical(readBin(con, what = "raw", n = 28), eof)
}
