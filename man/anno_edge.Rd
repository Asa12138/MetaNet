% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2-2.manipulate.R
\name{anno_edge}
\alias{anno_edge}
\title{Use dataframe to annotate edges of an igraph}
\usage{
anno_edge(go, anno_tab, verbose = TRUE)
}
\arguments{
\item{go}{metanet an igraph object}

\item{anno_tab}{a dataframe using to annotate (with rowname or a name column)}

\item{verbose}{logical}
}
\value{
a annotated igraph object
}
\description{
Use dataframe to annotate edges of an igraph
}
\examples{
data("c_net")
anno <- data.frame("from" = "s__Pelomonas_puraquae", "to" = "s__un_g__Rhizobium", new_atr = "new")
anno_edge(co_net, anno) -> anno_net
}
\seealso{
Other manipulate: 
\code{\link{anno_vertex}()},
\code{\link{c_net_annotate}()},
\code{\link{c_net_filter}()},
\code{\link{c_net_load}()},
\code{\link{c_net_save}()},
\code{\link{c_net_union}()},
\code{\link{get_e}()},
\code{\link{get_n}()},
\code{\link{get_v}()},
\code{\link{is_metanet}()}
}
\concept{manipulate}
