% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6.modules.R
\name{module_detect}
\alias{module_detect}
\alias{c_net_module}
\title{Detect the modules}
\usage{
module_detect(
  go,
  method = "cluster_fast_greedy",
  n_node_in_module = 0,
  delete = FALSE
)
}
\arguments{
\item{go}{an igraph object}

\item{method}{cluster_method: "cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass"}

\item{n_node_in_module}{transfer the modules less than n_node_in_module to "others"}

\item{delete}{logical, delete others modules? default:FALSE, the others module will be "others".}
}
\value{
an igraph object
}
\description{
Detect the modules
}
\examples{
data("c_net")
module_detect(co_net) -> co_net_modu
}
\seealso{
Other module: 
\code{\link{filter_n_module}()},
\code{\link{get_community}()},
\code{\link{get_module_eigen}()},
\code{\link{get_module}()},
\code{\link{module_eigen}()},
\code{\link{module_net}()},
\code{\link{summary_module}()},
\code{\link{to_module_net}()},
\code{\link{zp_analyse}()}
}
\concept{module}
