% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2-1.build.R
\name{multi_net_build}
\alias{multi_net_build}
\title{Multi-omics network build}
\usage{
multi_net_build(
  ...,
  mode = "full",
  method = "spearman",
  filename = FALSE,
  p.adjust.method = NULL,
  r_threshold = 0.6,
  p_threshold = 0.05,
  use_p_adj = TRUE,
  delete_single = TRUE
)
}
\arguments{
\item{...}{some omics abundance tables}

\item{mode}{"full"}

\item{method}{"spearman" or "pearson"}

\item{filename}{the prefix of saved .corr file or FALSE}

\item{p.adjust.method}{see \code{\link[stats]{p.adjust}}}

\item{r_threshold}{r_threshold (default: >0.6)}

\item{p_threshold}{p_threshold (default: <0.05)}

\item{use_p_adj}{use the p.adjust instead of p-value (default: TRUE)}

\item{delete_single}{should delete single vertexes?}
}
\value{
metanet
}
\description{
Multi-omics network build
}
\examples{
data("multi_test")
multi1 <- multi_net_build(list(Microbiome = micro, Metabolome = metab, Transcriptome = transc))
multi1 <- c_net_set(multi1, micro_g, metab_g, transc_g,
  vertex_class = c("Phylum", "kingdom", "type")
)
multi1 <- c_net_set(multi1, data.frame("Abundance1" = colSums(micro)),
  data.frame("Abundance2" = colSums(metab)), data.frame("Abundance3" = colSums(transc)),
  vertex_size = paste0("Abundance", 1:3)
)
c_net_plot(multi1)
}
\seealso{
Other build: 
\code{\link{c_net_build}()},
\code{\link{c_net_from_edgelist}()},
\code{\link{c_net_set}()},
\code{\link{c_net_update}()}
}
\concept{build}
