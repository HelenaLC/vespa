% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotXY.R
\name{plotXY}
\alias{plotXY}
\title{spatial plots}
\arguments{
\item{x}{\code{data.frame} or \code{\link{SingleCellExperiment}}.}

\item{xy}{character vector;
specifies the variables containing 
spatial (x- and y-)coordinates.}

\item{BPPARAM}{a \code{\link[BiocParallel]{BiocParallelParam}}
instance defining the parallel back-end to be used.}

\item{color_by, facet_by}{character vector; 
specifies the variable(s) to color/facet by.
Note than, if \code{color/facet_by} corresponds to more than 
one variable, you can only \code{facet/color_by} a single variable.
Also, the variables specified by \code{color_by}
must be of a unique type (e.g., all numeric).}

\item{highlight}{logical, integer or character vector;
specifies which points to highlight (i.e., enlarge).
Alternatively, a character string corresponding 
to a logical variable already stored in \code{x}.}

\item{assay}{character string; specifies which \code{assay} to use
when coloring by features. Valid values are \code{assayNames(x)}. 
Ignored unless \code{x} is a \code{SingleCellExperiment}.}

\item{scale}{function; applied to variables specified by \code{color_by}.
Should accept numeric values and return an equal-length numeric vector.}

\item{scale_xy}{character string; 
specifies whether to rescale coordinates, e.g., in order to make panels 
comparable when variables to \code{facet_by} occupy different spatial 
regions. If \code{TRUE} (the default), each group's coordinates are scaled 
separately (retaining their original aspect ratio!); ignored if \code{facet_by = NULL}.}

\item{thm}{character string; specifies the plot background.
When \code{"white"}, a border will be added as well.}

\item{size, alpha}{scalar numeric; aesthetics passed to \code{geom_point}.}

\item{nrow, ncol}{scalar integer; number of rows/columns for faceting.
Ignored when faceting by exactly two variables, in which case 
rows/columns will correspond to the first/second variable.}

\item{pal_continuous, pal_discrete, pal_logical}{character string; 
color palettes to use when values to \code{color_by}
are continuous, discrete and logical, respectively.}
}
\value{
ab object of class \code{\link[ggplot2]{ggplot}}.
}
\description{
Flexible \code{ggplot}-based spatial plots (by xy-coordinates)
with support for continuous/discrete/logical coloring, highlighting 
(via differences in point size), and faceting (by variables of interest).
}
\examples{
data(foo)

# color cells by compartment
plotXY(foo, color_by = "subset")

# color cells by whether or not they
# belong to the epithelial compartment
foo$epi <- foo$subset == "epi"
plotXY(foo, color_by = "epi")

# scale expression to values b/w 0 and 1 using 
# low (1\%) and high (99\%) quantiles as boundaries
.scale01 <- \(x, q = 0.01) {
  qs <- quantile(x, c(q, 1-q))
  x <- (x-qs[1]) / diff(qs)
  x[x < 0 | is.na(x)] <- 0
  x[x > 1] <- 1; return(x)
} 
# color by scaled expression values
# highlighting epithelial cells
plotXY(foo, 
  color_by = c("AQP8", "EPCAM"), 
  scale = .scale01, highlight = "epi")
  
# demo of not/scaling xy-coordinates
foo$foo <- foo$x > median(foo$x)
plotXY(foo, facet_by = "foo", scale_xy = TRUE)
plotXY(foo, facet_by = "foo", scale_xy = FALSE)

}
\author{
Helena L. Crowell
}
