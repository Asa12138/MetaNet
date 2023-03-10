% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{c_net_lay}
\alias{c_net_lay}
\title{Layout coordinates}
\usage{
c_net_lay(
  go,
  method = igraph::nicely(),
  order_by = NULL,
  order_ls = NULL,
  seed = 1234,
  line_curved = 1
)
}
\arguments{
\item{go}{igraph}

\item{method}{as_star(), as_tree(), in_circle(), nicely(), on_grid(), on_sphere(),randomly(), with_dh(), with_fr(), with_gem(), with_graphopt(), with_kk(),with_lgl(), with_mds(),as_line(), as_arc(), as_polygon(), as_polyarc(). see \code{\link[igarph]{layout_}}.}

\item{order_by}{order nodes according to a node attribute}

\item{order_ls}{manual the discrete variable with a vector, or continuous variable with "desc" to decreasing}

\item{seed}{random seed}

\item{line_curved}{consider line curved, only for some layout methods like as_line(), as_polygon().default:0}
}
\value{
coordinates for nodes,columns: name, X, Y
}
\description{
Layout coordinates
}
\examples{
c_net_lay(co_net)->coors
c_net_plot(co_net,coors)
c_net_plot(co_net,c_net_lay(co_net,in_circle()),vertex.size=2)
c_net_plot(co_net,c_net_lay(co_net,in_circle(),order_by="v_class"),vertex.size=2)
c_net_plot(co_net,c_net_lay(co_net,in_circle(),order_by="size",order_ls="desc"))
c_net_plot(co_net,c_net_lay(co_net,as_polygon(3)))
}
