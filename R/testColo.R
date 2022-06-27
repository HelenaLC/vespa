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
#' @param group character vector; 
#'   optional additional variable(s) to group by, e.g., 
#'   in order to perform a series of independent tests.
#' @param nc scalar integer; 
#'   threshold on the number of cells above 
#'   which to consider a group for testing.
#'
#' @details ...
#'
#' @return a \code{\link[S4Vectors]{DataFrame-class}} object.
#'
#' @examples
#' data(foo)
#' res <- testColo(foo, 
#'   by = "subset", 
#'   gridsize = c(5, 5))
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
  xy = c("x", "y"), nc = 20, 
  BPPARAM = SerialParam(), ...)
{
  # validity checks
  .check_xy(x, xy)
  stopifnot(
    is(BPPARAM, "BiocParallelParam"),
    is.numeric(nc), length(nc) == 1, round(nc) == nc)
  
  df <- .df(x)
  if (is.factor(df[[by]]))
    df[[by]] <- droplevels(df[[by]])
  
  lys <- if (is.null(group)) {
    list(foo = df)
  } else {
    args <- c(as.list(df[group]), sep = ";")
    df$.group <- do.call(paste, args)
    split(df, df$.group)
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
        n <- nrow(a <- lys[[g]][[i]])
        m <- nrow(b <- lys[[g]][[j]])
        if (n < nc | m < nc) return(NULL)
        c <- rbind(a[, xy], b[, xy])
        r <- colRanges(as.matrix(c))
        l <- r[, 1]; r <- r[, 2]
        d1 <- c(kde(a, xmin = l, xmax = r, ...)$estimate)
        d2 <- c(kde(b, xmin = l, xmax = r, ...)$estimate)
        corr <- cor(d1, d2, method = "pearson")
        data.frame(.group = g, i, j, corr, n, m)
      })
    rmv <- vapply(pcc, is.null, logical(1))
    pcc <- do.call(rbind, pcc[!rmv])
  })
  rmv <- vapply(pcc, is.null, logical(1))
  pcc <- do.call(rbind, pcc[!rmv])
  
  if (!is.null(group)) {
    # add uniquely mappable metadata
    ids <- unique(pcc$.group)
    tmp <- vapply(ids, \(.) {
      i <- !is.na(match(df$.group, .))
      vapply(df[i, ], \(.) length(unique(.)) == 1, logical(1))
    }, logical(ncol(df)))
    j <- names(df)[rowSums(tmp) == length(ids)]
    j <- setdiff(j, c(".group", group, by))
    fd <- df[match(pcc$.group, df$.group), ]
    pcc <- cbind(fd[group], pcc, fd[j])
    # drop grouping column 
    pcc$.group <- NULL
  }
  rownames(pcc) <- NULL
  return(pcc)
}
