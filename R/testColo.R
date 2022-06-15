#' @name testColo
#' @rdname testColo
#' 
#' @title co-localization
#'
#' @description ...
#'
#' @template man
#' 
#' @param ... further arguments; passed to \code{\link[ks]{kde}}.
#'
#' @details ...
#'
#' @return a \code{\link[S4Vectors]{DataFrame-class}} object.
#'
#' @examples
#' data(foo)
#' testColo(foo, by = "subset", gridsize = c(5, 5))
#' 
#' @author Helena L. Crowell
#' 
#' @importFrom ks kde
#' @importFrom S4Vectors DataFrame
#' @importFrom matrixStats colRanges
#' @importFrom BiocParallel bpmapply SerialParam
#' @export
NULL

testColo <- \(x, 
  by, xy = c("x", "y"),
  BPPARAM = SerialParam(),
  ...)
{
  # validity checks
  stopifnot(
    is(BPPARAM, "BiocParallelParam"))
  .check_xy(x, xy)
  
  df <- .df(x)
  
  lys <- split(df[, xy], df[[by]])
  idx <- combn(length(lys), 2)
  pcc <- bpmapply(
    SIMPLIFY = FALSE,
    BPPARAM = BPPARAM,
    i = idx[1, ], j = idx[2, ],
    \(i, j) {
      a <- lys[[i]]; b <- lys[[j]]
      if (nrow(a) < 2 | 
          nrow(b) < 2) 
        return(NULL)
      m <- rbind(a[, xy], b[, xy])
      r <- colRanges(as.matrix(m))
      l <- r[, 1]; r <- r[, 2]
      d1 <- c(kde(a, xmin = l, xmax = r)$estimate)
      d2 <- c(kde(b, xmin = l, xmax = r)$estimate)
      corr <- cor(d1, d2, method = "pearson")
      data.frame(from = names(lys)[i], to = names(lys)[j], corr)
    })
  rmv <- vapply(pcc, is.null, logical(1))
  pcc <- do.call(rbind, pcc[!rmv])
  return(pcc)
}
