#' @name testQuadrat
#' @rdname testQuadrat
#' 
#' @title quadrat test
#'
#' @description 
#' Test for complete spatial randomness (CSR) of points based on 
#' quadrat counts (see \code{\link[spatstat.core]{quadrat.test}}).
#'
#' @template man
#' 
#' @param by character vector; 
#'   specifies a variable or set of variables to group points by. 
#'   Each grouping will be tested separately.
#' @param ... further arguments;
#'   passed to \code{\link[spatstat.core]{quadrat.test}}.
#'
#' @details 
#' In brief, the window of observation is divided into \eqn{n} tiles, 
#' and the observed number of data points \eqn{O} in each tile \eqn{i} 
#' is counted. Then, the expected number of points in each quadrat \eqn{E} 
#' is calculated, and the Pearson statistic computed according to: 
#' \deqn{\chi^2=\sum_{i=1}^n((O_i-E)^2/E)}
#'
#' @return a \code{\link[S4Vectors]{DataFrame-class}} object.
#'
#' @examples
#' data(foo)
#' testQuadrat(foo, by = "subset")
#' 
#' @author Helena L. Crowell
#' 
#' @importFrom methods is
#' @importFrom spatstat.geom owin ppp
#' @importFrom spatstat.core quadrat.test
#' @export
NULL

testQuadrat <- \(x, 
  by = NULL, 
  xy = c("x", "y"), 
  n = 100, ...)
{
  # validity checks
  .check_xy(x, xy)
  df <- .df(x)
    
  if (is.null(by)) {
    idx <- gl(1, nrow(df)) 
  } else {
    idx <- df[, by]
    if (length(by) > 1)
      idx <- as.list(idx)
  }
  res <- by(df, idx, simplify = FALSE, \(.) {
    if ((.n <- nrow(.)) < n) return(NULL)
    x <- .[[xy[1]]]
    y <- .[[xy[2]]]
    w <- owin(range(x), range(y))
    pp <- ppp(x, y, w)
    qt <- quadrat.test(pp)
    DataFrame(
      row.names = NULL, 
      .[1, by, drop = FALSE],
      n = .n, 
      p.value = qt$p.value,
      statistic = qt$statistic,
      expected = qt$expected[1],
      observed = I(list(qt$observed)),
      residuals = I(list(qt$residuals)))
  })
  rmv <- vapply(res, is.null, logical(1))
  res <- do.call(rbind, res[!rmv])
  return(res)
}
