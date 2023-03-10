% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{c_net_build}
\alias{c_net_build}
\title{Construct a network from correlation table}
\usage{
c_net_build(
  cor_res,
  r_thres = 0.6,
  p_thres = 0.05,
  use_p_adj = T,
  del_single = T
)
}
\arguments{
\item{cor_res}{result of c_net_cal()}

\item{r_thres}{r_threshold (default:>0.6)}

\item{p_thres}{p_threshold (default:<0.05)}

\item{use_p_adj}{use the p.adjust instead of p-value (default: T)}

\item{del_single}{should delete single vertexes?}
}
\value{
an igraph object
}
\description{
Construct a network from correlation table
}
\examples{
data("otutab")
t(otutab) -> totu
metadata[,3:10] -> env
c_net_cal(totu) -> corr
c_net_build(corr,r_thres=0.65) -> co_net
c_net_cal(totu,env) -> corr2
c_net_build(corr2) -> co_net2
}
