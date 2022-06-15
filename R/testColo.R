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
#' @param by character string; 
#'   specifies the variable to test by for co-localization.
#' @param group character string; 
#'   an option additional variable to group by, e.g., 
#'   in order to perform a series of independent tests.
#' @param n scalar integer; 
#'   threshold on the number of cells above 
#'   which to consider a group for testing.
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

testColo <- \(x, by, 
  group = NULL, 
  xy = c("x", "y"), n = 20, 
  BPPARAM = SerialParam(), ...)
{
  # validity checks
  .check_xy(x, xy)
  stopifnot(
    is(BPPARAM, "BiocParallelParam"),
    is.numeric(n), length(n) == 1, round(n) == n)

  df <- .df(x)
  lys <- if (is.null(group)) {
    list(foo = df)
  } else {
    split(df, df[[group]])
  }
  lys <- lapply(lys, \(.) split(.[, xy], .[[by]]))
  idx <- lapply(lys, \(.) combn(names(.), 2))
  
  names(gs) <- gs <- names(lys)
  pcc <- lapply(gs, \(g) {
    pcc <- bpmapply(
      SIMPLIFY = FALSE,
      BPPARAM = BPPARAM,
      i = idx[[g]][1, ], 
      j = idx[[g]][2, ],
      \(i, j) {
        a <- lys[[g]][[i]]
        b <- lys[[g]][[j]]
        if (nrow(a) < n | nrow(b) < n) 
          return(NULL)
        m <- rbind(a[, xy], b[, xy])
        r <- colRanges(as.matrix(m))
        l <- r[, 1]; r <- r[, 2]
        d1 <- c(kde(a, xmin = l, xmax = r, ...)$estimate)
        d2 <- c(kde(b, xmin = l, xmax = r, ...)$estimate)
        corr <- cor(d1, d2, method = "pearson")
        data.frame(group = g, from = i, to = j, corr)
      })
    rmv <- vapply(pcc, is.null, logical(1))
    pcc <- do.call(rbind, pcc[!rmv])
  })
  rmv <- vapply(pcc, is.null, logical(1))
  pcc <- do.call(rbind, pcc[!rmv])
  rownames(pcc) <- NULL
  if (is.null(group)) {
    pcc$group <- NULL
  } else {
    rnm <- grep("group", names(pcc))
    names(pcc)[rnm] <- group
  }
  head(pcc)
}
