#' Print method for 'corr' objects
#'
#' @param x 'corr' object
#' @param ... additional arguments
#'
#' @exportS3Method
#' @method print corr
#' @return No value
print.corr <- function(x, ...) {
  cat("Correlation table:\n")
  cat("Table dimensions:", nrow(x$r), "rows,", ncol(x$r), "columns\n")
}


#' Print method for 'coors' objects
#'
#' @param x 'coors' object
#' @param ... additional arguments
#' @return No value
#' @exportS3Method
#' @method print coors
print.coors <- function(x, ...) {
  cat("Coordinates:\n")
  cat("Table dimensions:", nrow(x$coors), "rows,", ncol(x$coors), "columns\n")
}


#' Print method for 'metanet' objects
#'
#' @param x 'metanet' object
#' @param ... Additional arguments
#' @exportS3Method
#' @method print metanet
#' @return No value
print.metanet <- function(x, ...) {
  pcutils::dabiao(" metanet ", print = TRUE)
  print.igraph(x)
}

#' Print method for 'ggig' objects
#'
#' @param x 'ggig' object
#' @param ... Additional arguments
#'
#' @return No value
#' @exportS3Method
#' @method print ggig
print.ggig <- function(x, ...) {
  pcutils::dabiao(" ggig ", print = TRUE)
  pcutils::dabiao("use `plot()` to visualize ggig object.", char = " ", print = TRUE)
}

#' Print method for 'robust' objects
#' @param x 'robust' object
#' @param ... Additional arguments
#' @return No value
#' @exportS3Method
#' @method print robust
print.robust <- function(x, ...) {
  pcutils::dabiao(" robust ", print = TRUE)
  pcutils::dabiao("use `plot()` to visualize robust object.", char = " ", print = TRUE)
}

#' Print method for 'robustness' objects
#' @param x 'robustness' object
#' @param ... Additional arguments
#' @return No value
#' @exportS3Method
#' @method print robustness
print.robustness <- function(x, ...) {
  pcutils::dabiao(" robustness ", print = TRUE)
  pcutils::dabiao("use `plot()` to visualize robustness object.", char = " ", print = TRUE)
}

#' Print method for 'vulnerability' objects
#' @param x 'vulnerability' object
#' @param ... Additional arguments
#' @return No value
#' @exportS3Method
#' @method print vulnerability
print.vulnerability <- function(x, ...) {
  pcutils::dabiao(" vulnerability ", print = TRUE)
  pcutils::dabiao("use `plot()` to visualize vulnerability object.", char = " ", print = TRUE)
}

#' Print method for 'cohesion' objects
#' @param x 'cohesion' object
#' @param ... Additional arguments
#' @return No value
#' @exportS3Method
#' @method print cohesion
print.cohesion <- function(x, ...) {
  pcutils::dabiao(" cohesion ", print = TRUE)
  pcutils::dabiao("use `plot()` to visualize cohesion object.", char = " ", print = TRUE)
}
