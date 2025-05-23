% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3-1.layout.R
\name{g_layout_nice}
\alias{g_layout_nice}
\alias{get_big_lay_nice}
\alias{g_layout_circlepack}
\alias{g_layout_treemap}
\alias{g_layout_backbone}
\alias{g_layout_stress}
\title{Layout with group nicely}
\usage{
g_layout_nice(go, group = "module", mode = "circlepack", group_zoom = 1, ...)

get_big_lay_nice(go, group = "module", mode = "circlepack", ...)

g_layout_circlepack(go, group = "module", ...)

g_layout_treemap(go, group = "module", ...)

g_layout_backbone(go, group = "module", ...)

g_layout_stress(go, group = "module", ...)
}
\arguments{
\item{go}{igraph or metanet}

\item{group}{group name (default: module)}

\item{mode}{circlepack, treemap, backbone, stress}

\item{group_zoom}{zoom for each group}

\item{...}{add}
}
\value{
coors
}
\description{
Layout with group nicely
}
\examples{
\donttest{
data("c_net")
module_detect(co_net, method = "cluster_fast_greedy") -> co_net_modu
if (requireNamespace("ggraph")) {
  plot(co_net_modu, coors = g_layout_nice(co_net_modu, group = "module"))
  plot(co_net_modu, coors = g_layout_nice(co_net_modu, group = "module", mode = "treemap"))
}
}
}
\seealso{
Other g_layout: 
\code{\link{g_layout_poly_sector}()},
\code{\link{g_layout_polygon}()},
\code{\link{g_layout}()}
}
\concept{g_layout}
