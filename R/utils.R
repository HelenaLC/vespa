#' @importFrom methods is
#' @importFrom SummarizedExperiment colData
.check_xy <- \(x, xy) {
  stopifnot(
    !is.na(xy),
    length(xy) == 2,
    is.character(xy))
  if (is(x, "SingleCellExperiment")) x <- colData(x)
  stopifnot(sum(!is.na(match(xy, names(x)))) == 2)
}

#' @importFrom methods is
#' @importFrom SummarizedExperiment assayNames
.check_assay <- \(x, assay) {
  if (is.null(assay) ||
      !is(x, "SingleCellExperiment"))
    return()
  stopifnot(!is.na(assay),
    length(assay) == 1,
    is.character(assay),
    sum(grepl(assay, assayNames(x))) == 1)
}

#' @importFrom methods is
#' @importFrom SummarizedExperiment assay colData
.df <- \(x, assay = NULL, ...) {
  if (is(x, "SingleCellExperiment")) {
    a <- if (!is.null(assay)) {
      a <- assay(x, assay)
      if (!is.matrix(a))
        a <- as.matrix(a)
      i <- intersect(c(...), rownames(a))
      a <- if (length(i) > 0) 
        t(a[i, , drop = FALSE])
    }
    if (is.null(a))
      a <- matrix(0, ncol(x), 0)
    x <- data.frame(
      colData(x), a,
      row.names = NULL,
      check.names = FALSE)
  } else if (is(x, "DFrame")) {
    x <- as.data.frame(x)
  }
  if (!is.data.frame(x))
    stop(
      "'x' must be of class 'data.frame', ",
      "'DFrame' or 'SingleCellExperiment'.")
  return(x)
}
