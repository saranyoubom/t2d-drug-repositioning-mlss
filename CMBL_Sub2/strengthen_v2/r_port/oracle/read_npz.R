## read_npz.R - a minimal reader for NumPy .npz / .npy files (little-endian float64, int64, int32, bool; C or Fortran order), to compare the R results with the arrays the Python scripts saved (positive_control/inj_<tag>.npz).
## read_npy(raw) -> an R array (0-d arrays become a plain number); read_npz(path) -> a named list of arrays. No numpy needed.
read_npy <- function(b) {
  stopifnot(identical(b[1:6], as.raw(c(0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59)))); major <- as.integer(b[7])
  hl <- if (major == 1L) readBin(b[9:10], "integer", n = 1, size = 2, signed = FALSE, endian = "little") else readBin(b[9:12], "integer", n = 1, size = 4, endian = "little"); off <- if (major == 1L) 10L else 12L
  hdr <- rawToChar(b[(off + 1L):(off + hl)]); descr <- sub(".*'descr': *'([^']*)'.*", "\\1", hdr); fortran <- grepl("'fortran_order': *True", hdr)
  shp <- sub(".*'shape': *\\(([^)]*)\\).*", "\\1", hdr); shape <- if (nzchar(trimws(shp))) as.integer(strsplit(gsub(" ", "", shp), ",", fixed = TRUE)[[1]]) else integer(0); shape <- shape[!is.na(shape)]
  n <- if (length(shape)) prod(shape) else 1L; data <- b[(off + hl + 1L):length(b)]
  x <- switch(descr, "<f8" = readBin(data, "double", n = n, size = 8, endian = "little"), "<i8" = readBin(data, "double", n = n, size = 8, endian = "little"), "<i4" = readBin(data, "integer", n = n, size = 4, endian = "little"),
              "|b1" = as.logical(as.integer(data[seq_len(n)])), stop("dtype not supported: ", descr))
  if (descr == "<i8") x <- as.numeric(readBin(data, "integer", n = 2L * n, size = 4, endian = "little")[c(TRUE, FALSE)]) + 4294967296 * readBin(data, "integer", n = 2L * n, size = 4, endian = "little")[c(FALSE, TRUE)]
  if (length(shape) <= 1L) return(if (length(shape) == 0L) x[1] else x)
  if (fortran) array(x, dim = shape) else aperm(array(x, dim = rev(shape)), rev(seq_along(shape)))
}
read_npz <- function(path) {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE)); fl <- unzip(path, exdir = tmp); out <- list()
  for (f in fl) out[[sub("[.]npy$", "", basename(f))]] <- read_npy(readBin(f, "raw", file.size(f)))
  out
}
