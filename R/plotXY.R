#' @name plotXY
#' @rdname plotXY
#'
#' @title spatial plots
#'
#' @description 
#' Flexible \code{ggplot}-based spatial plots (by xy-coordinates)
#' with support for continuous/discrete/logical coloring, highlighting 
#' (via differences in point size), and faceting (by variables of interest).
#'
#' @template man
#' 
#' @param color_by,facet_by character vector; 
#'   specifies the variable(s) to color/facet by.
#'   Note than, if \code{color/facet_by} corresponds to more than 
#'   one variable, you can only \code{facet/color_by} a single variable.
#'   Also, the variables specified by \code{color_by}
#'   must be of a unique type (e.g., all numeric).
#' @param highlight logical, integer or character vector;
#'   specifies which points to highlight (i.e., enlarge).
#'   Alternatively, a character string corresponding 
#'   to a logical variable already stored in \code{x}.
#' @param assay character string; specifies which \code{assay} to use
#'   when coloring by features. Valid values are \code{assayNames(x)}. 
#'   Ignored unless \code{x} is a \code{SingleCellExperiment}.
#' @param scale function; applied to variables specified by \code{color_by}.
#'   Should accept numeric values and return an equal-length numeric vector.
#' @param thm character string; specifies the plot background.
#'   When \code{"white"}, a border will be added as well.
#' @param size,alpha scalar numeric; aesthetics passed to \code{geom_point}.
#' @param nrow,ncol scalar integer; number of rows/columns for faceting.
#'   Ignored when faceting by exactly two variables, in which case 
#'   rows/columns will correspond to the first/second variable.
#' @param pal_continuous,pal_discrete,pal_logical character string; 
#'   color palettes to use when values to \code{color_by}
#'   are continuous, discrete and logical, respectively.
#'
#' @return ab object of class \code{\link[ggplot2]{ggplot}}.
#'
#' @examples
#' data(foo)
#' 
#' # color cells by compartment
#' plotXY(foo, color_by = "subset")
#' 
#' # color cells by whether or not they
#' # belong to the epithelial compartment
#' foo$epi <- foo$subset == "epi"
#' plotXY(foo, color_by = "epi")
#' 
#' # scale expression to values b/w 0 and 1 using 
#' # low (1%) and high (99%) quantiles as boundaries
#' .scale01 <- \(x, q = 0.01) {
#'   qs <- quantile(x, c(q, 1-q))
#'   x <- (x-qs[1]) / diff(qs)
#'   x[x < 0 | is.na(x)] <- 0
#'   x[x > 1] <- 1; return(x)
#' } 
#' # color by scaled expression values
#' # highlighting epithelial cells
#' plotXY(foo, 
#'   color_by = c("AQP8", "EPCAM"), 
#'   scale = .scale01, highlight = "epi")
#'
#' @author Helena L. Crowell
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @export
NULL

.mar <- margin(t = 0.25, b = 0.25, unit = "line")
.thm_black <- theme(
  panel.background = element_rect(fill = "black"),
  strip.text = element_text(color = "white", margin = .mar),
  strip.background = element_rect(fill = "black", color = "black"))
.thm_white <- theme(
  panel.background = element_rect(fill = "white"),
  strip.text = element_text(color = "black", margin = .mar),
  strip.background = element_rect(fill = "lightgrey", color = "white"))

plotXY <- \(x,
  xy = c("x", "y"),
  color_by = NULL,
  facet_by = NULL,
  highlight = NULL,
  assay = "logcounts", 
  scale = NULL,
  size = 1, alpha = 1, 
  nrow = NULL, ncol = NULL,
  thm = c("black", "white"),
  pal_continuous = brewer.pal(9, "OrRd"),
  pal_discrete = brewer.pal(12, "Paired"),
  pal_logical = c("pink", "deeppink"))
{
  # argument validation
  thm <- match.arg(thm)
  .check_assay(x, assay)
  .check_xy(x, xy)
  
  if (length(color_by) > 1 && length(facet_by) > 1)
    stop("Can only 'facet_by' one when provided", 
      " with more than one value to 'color_by'.")
  
  val <- if (is(x, "SingleCellExperiment")) c(rownames(x), names(colData(x))) else names(x)
  if (!is.null(color_by)) stopifnot(is.character(color_by), color_by %in% val)
  if (!is.null(facet_by)) stopifnot(is.character(facet_by), facet_by %in% val)
  if (!is.null(scale)) {
    y <- runif(n <- 100)
    z <- tryCatch(scale(y), error = function(e) e)
    . <- c(is.function(scale), !inherits(z, "error"), length(z) == n, is.numeric(z))
    if (any(!.)) stop(
      "'scale' should be a function applicable to numeric values and ",
      "that returns a numeric vector of the same length as its input.")
  }

  df <- .df(x, assay, color_by)
  aes <- list()
  
  if (!is.null(highlight)) {
    if (is.character(highlight)) {
      df$.size_by <- df[[highlight]]
      stopifnot(is.logical(df$.size_by))
    } else {
      df$.size_by <- seq(nrow(df))
      fd <- df[highlight, ]$.size_by
      df$.size_by <- df$.size_by %in% fd
    }
    aes <- c(aes, scale_size_manual(values = c("FALSE" = 0.5, "TRUE" = 1.5)))
    size_by <- ".size_by"
  } else {
    size_by <- NULL
  }
  
  typ <- unique(vapply(df[, color_by], class, character(1)))
  
  if (!is.null(scale)) {
    if (!all(typ == "numeric"))
      stop("Cannot apply 'scale' function on ", 
        "non-numeric variables to 'color_by'.")
    for (. in color_by) df[[.]] <- scale(df[[.]])
  }
  
  if ((n <- length(color_by)) > 1) {
    if (length(typ) != 1)
      stop("Values to 'color_by' must be of",
        " one type; found columns of type: ", 
        paste(sQuote(typ), collapse = ", "))
    i <- setdiff(names(df), color_by)
    fd <- df[, i, drop = FALSE]
    fd <- do.call(rbind, replicate(n, fd, FALSE))
    fd$.name <- rep(color_by, each = nrow(df))
    fd$.value <- unlist(df[, color_by])
    color_by <- ".value"
    facet_by <- ".name"
    df <- fd; rm(fd)
  }
  
  aes <- c(aes, if (!is.null(color_by)) {
    if (scale_type(df[[color_by]]) == "continuous") {
      list(
        theme(
          legend.position = "bottom",
          legend.key.width = unit(1.5, "lines"),
          legend.key.height = unit(0.5, "lines"),
          legend.title = element_text(vjust = 1)),
        scale_colour_gradientn(
          if (color_by != ".value") color_by,
          colors = pal_continuous))
    } else {
      if (is.logical(df[[color_by]])) {
        list(
          guides(col = "none"),
          scale_color_manual(NULL, values = pal_logical))
      } else {
        list(
          theme(legend.key.size = unit(0.5, "lines")),
          scale_color_manual(color_by, values = pal_discrete),
          guides(col = guide_legend(override.aes = list(alpha = 1, size = 2))))
      }
    }
  } else {
    df$.color_by <- switch(thm, black = "white", "black")
    color_by <- ".color_by"; scale_color_identity(NULL)
  })
  
  ggplot(df, aes_string(xy[1], xy[2], 
    col = color_by, size = size_by)) +
    geom_point(
      shape = 16, alpha = alpha, 
      if (!is.null(size_by) && size_by != ".size_by") size = size) + 
    guides(size = "none") +
    ( if (!is.null(facet_by)) {
      if (length(facet_by) == 2) {
        facet_grid(facet_by, switch = "y")
      } else facet_wrap(facet_by, nrow = nrow, ncol = ncol) 
    } ) +
    coord_fixed() + theme_void() + 
    aes + get(paste0(".thm_", thm)) +
    theme(strip.text = element_text(face = "bold"))
}
