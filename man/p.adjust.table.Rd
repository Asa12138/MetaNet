% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate.R
\name{p.adjust.table}
\alias{p.adjust.table}
\title{p.adjust apply on a correlation table (matrix or data.frame)}
\usage{
p.adjust.table(pp, method = "BH")
}
\arguments{
\item{pp}{table}

\item{method}{see \code{\link[stats]{p.adjust}}}
}
\value{
matrix
}
\description{
p.adjust apply on a correlation table (matrix or data.frame)
}
\examples{
matrix(abs(rnorm(100,0.01)),10,10)->pp
p.adjust.table(pp)
}
