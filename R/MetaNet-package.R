#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import ggplot2
#' @import igraph
#' @importFrom dplyr select filter left_join arrange count mutate group_by summarise pull rename
#' @importFrom igraph V E
#' @importFrom utils data combn  head tail write.csv
#' @importFrom stats aggregate median var sd setNames runif relevel coef fitted cor coefficients formula lm pt
#' @importFrom stats time na.omit kmeans p.adjust density approx ks.test smooth.spline
#' @importFrom pcutils lib_ps mmscale get_cols trans guolv hebing update_param tidai
#' @importFrom graphics legend hist par text
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom reshape2 acast melt
## usethis namespace: end
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))

#' Show MetaNet logo
#'
#' @return picture
#' @export
show_MetaNet_logo <- function() {
    pcutils::read.file(system.file("figures/hexSticker1.png", package = "MetaNet"), all_yes = TRUE)
}
