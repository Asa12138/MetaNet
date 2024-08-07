% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1.calculate.R
\name{c_net_calculate}
\alias{c_net_calculate}
\alias{c_net_cal}
\title{Calculate correlation for one or two t(otutab), or distance for one t(otutab).}
\usage{
c_net_calculate(
  totu,
  totu2 = NULL,
  method = "spearman",
  filename = FALSE,
  p.adjust.method = NULL,
  p.adjust.mode = "all",
  threads = 1,
  verbose = TRUE
)
}
\arguments{
\item{totu}{t(otutab), row are samples, column are features.}

\item{totu2}{t(otutab2) or NULL, row are samples, column are features.}

\item{method}{"spearman" (default), "pearson", "sparcc", or distance index from \code{\link[vegan]{vegdist}}.}

\item{filename}{the prefix of saved .corr file or FALSE.}

\item{p.adjust.method}{see \code{\link[stats]{p.adjust}}}

\item{p.adjust.mode}{see \code{\link{p.adjust.table}}}

\item{threads}{threads, default: 1.}

\item{verbose}{verbose, default: TRUE.}
}
\value{
a corr object with 3 elements:
\item{r}{default: spearman correlation}
\item{p.value}{default: p-value of spearman correlation}
\item{p.adjust}{default p.adjust.method = NULL}
}
\description{
Calculate correlation for one or two t(otutab), or distance for one t(otutab).
}
\examples{
data("otutab", package = "pcutils")
t(otutab) -> totu
c_net_calculate(totu) -> corr
metadata[, 3:10] -> env
c_net_calculate(totu, env) -> corr2
}
\seealso{
Other calculate: 
\code{\link{cal_sim}()},
\code{\link{fast_cor}()},
\code{\link{p.adjust.table}()},
\code{\link{read_corr}()}
}
\concept{calculate}
