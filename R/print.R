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
