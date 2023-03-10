% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{g_lay_polygon}
\alias{g_lay_polygon}
\alias{g_lay_polyarc}
\title{Layout with group as a polygon}
\usage{
g_lay_polygon(
  go,
  group = "v_group",
  group2 = NULL,
  group2_order = NULL,
  line_curved = 0.5
)

g_lay_polyarc(
  go,
  group = "v_group",
  group2 = NULL,
  group2_order = NULL,
  space = pi/4
)
}
\arguments{
\item{go}{igraph}

\item{group}{group name (default:v_group)}

\item{group2}{group2 name, will order nodes in each group according to group2_order (default:v_class)}

\item{group2_order}{group2_order}
}
\value{
coors
}
\description{
Layout with group as a polygon

Layout with group as a polyarc
}
\examples{
g_lay_polygon(multi1)->oridata
c_net_plot(multi1,oridata)
g_lay_polyarc(multi1,group2_order = c(LETTERS[4:1]))->oridata
c_net_plot(multi1,oridata)
}
