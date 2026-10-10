# Minimal MATLAB v5 .mat reader (numeric matrices only; others skipped). Returns named list.
read_mat5 <- function(f) {
  r <- readBin(f, "raw", file.size(f)); pos <- 129L; out <- list()
  szof <- c(`1`=1,`2`=1,`3`=2,`4`=2,`5`=4,`6`=4,`7`=4,`9`=8,`12`=8,`13`=8,`14`=8)
  rd <- function(b, p) { # returns list(type, data(raw), next)
    hi <- readBin(b[(p+2):(p+3)], "integer", 1, 2, signed=FALSE, endian="little")
    lo <- readBin(b[p:(p+1)], "integer", 1, 2, signed=FALSE, endian="little")
    if (hi != 0) { # small element
      ty <- lo; n <- hi
      return(list(ty=ty, d=b[(p+4):(p+3+n)], nx=p+8L))
    } else {
      w <- lo; n <- readBin(b[(p+4):(p+7)], "integer", 1, 4, endian="little")
      pad <- if (w == 15L) 0L else (8 - n %% 8) %% 8
      list(ty=w, d=if (n>0) b[(p+8):(p+7+n)] else raw(0), nx=p+8L+n+pad)
    }
  }
  num <- function(ty, d) switch(as.character(ty),
    `1`=as.integer(readBin(d,"integer",length(d),1,signed=TRUE)), `2`=as.integer(d),
    `3`=readBin(d,"integer",length(d)%/%2,2,endian="little"), `4`=readBin(d,"integer",length(d)%/%2,2,signed=FALSE,endian="little"),
    `5`=readBin(d,"integer",length(d)%/%4,4,endian="little"), `6`=readBin(d,"integer",length(d)%/%4,4,signed=FALSE,endian="little"),
    `7`=readBin(d,"numeric",length(d)%/%4,4,endian="little"), `9`=readBin(d,"numeric",length(d)%/%8,8,endian="little"),
    `12`=readBin(d,"integer",length(d)%/%8,8,endian="little"), as.numeric(NA))
  parse_matrix <- function(d) {
    e1 <- rd(d, 1L); fl <- readBin(e1$d[1:4], "integer", 1, 4, endian="little"); cls <- bitwAnd(fl, 255)
    e2 <- rd(d, e1$nx); dims <- num(e2$ty, e2$d)
    e3 <- rd(d, e2$nx); nm <- rawToChar(e3$d)
    if (cls %in% 6:13) {
      e4 <- rd(d, e3$nx); v <- num(e4$ty, e4$d)
      list(name=nm, val=array(as.numeric(v), dims))
    } else if (cls == 4) { # char
      e4 <- rd(d, e3$nx); list(name=nm, val=intToUtf8(num(e4$ty, e4$d)))
    } else list(name=nm, val=NULL)
  }
  while (pos < length(r)) {
    e <- rd(r, pos); pos <- e$nx
    d <- if (e$ty == 15L) memDecompress(e$d, "gzip") else e$d
    if (e$ty == 15L) { e2 <- rd(d, 1L); d <- e2$d; ty <- e2$ty } else ty <- e$ty
    if (ty == 14L) { m <- parse_matrix(d); if (!is.null(m$val)) out[[m$name]] <- m$val }
  }
  out
}
