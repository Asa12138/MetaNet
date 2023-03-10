% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{c_net_set}
\alias{c_net_set}
\title{Set basic attributes from totu table}
\usage{
c_net_set(
  go,
  ...,
  vertex_group = "v_group",
  vertex_class = "v_class",
  vertex_size = "size",
  edge_type = "e_type",
  edge_class = "e_class",
  edge_width = "width"
)
}
\arguments{
\item{go}{a igraph object}

\item{...}{some dataframes to annotate go}

\item{vertex_class}{choose which column to be vertex_class (map to vertex_color)}

\item{vertex_size}{choose which column to be vertex_size (map to vertex_size)}
}
\value{
a igraph object
}
\description{
Set basic attributes from totu table
}
\examples{
data("otutab")
t(otutab) -> totu
metadata[,3:10] -> env

data("c_net")
co_net<-c_net_set(co_net,taxonomy,data.frame("Abundance"=colSums(totu)),vertex_class="Phylum",vertex_size="Abundance")
co_net_rmt<-c_net_set(co_net_rmt,taxonomy,data.frame("Abundance"=colSums(totu)),vertex_class="Phylum",vertex_size="Abundance")
co_net2<-c_net_set(co_net2,taxonomy,data.frame(name=colnames(env),env=colnames(env)),vertex_class=c("Phylum","env"))
co_net2<-c_net_set(co_net2,data.frame("Abundance"=colSums(totu)),vertex_size="Abundance")
}
