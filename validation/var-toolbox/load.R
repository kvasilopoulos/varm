dir0 <- "/Users/VASILK04/Documents/research/varm/validation/var-toolbox"
rd <- function(n) { d <- read.csv(file.path(dir0, "data", paste0(n, ".csv")), header = FALSE, skip = 2, stringsAsFactors = FALSE)
  nm <- unlist(read.csv(file.path(dir0, "data", paste0(n, ".csv")), header = FALSE, nrows = 1, skip = 1, stringsAsFactors = FALSE))[-1]
  m <- as.matrix(d[, -1]); storage.mode(m) <- "double"; colnames(m) <- nm; m }
