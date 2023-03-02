#'@title Microbiota community data
#'@description Microbiota community data
#'
#'@docType data
#'@name otutab
#'@usage otutab
#'@format three co_nets
#'\describe{
#' \item{tbls}{contians otutable rawdata}
#' \item{metas}{contians metadata}
#' \item{otus}{contians taxomomy table}
#'}
#'
#'
#'
NULL

#'@title MetaNet networks
#'@description MetaNet co_nets
#'
#'@docType data
#'@name c_net
#'@usage c_net
#'@format three co_nets
#'\describe{
#' \item{tbls}{contians otutable rawdata}
#' \item{metas}{contians metadata}
#' \item{otus}{contians taxomomy table}
#'}
#'
#'
#'
NULL

#'@title MetaNet networks
#'@description MetaNet co_nets
#'
#'@docType data
#'@name multi_test
#'@usage c_net
#'@format three co_nets
#'\describe{
#' \item{tbls}{contians otutable rawdata}
#' \item{metas}{contians metadata}
#' \item{otus}{contians taxomomy table}
#'}
#'
#'
#'
NULL

#'@title MetaNet networks
#'@description MetaNet co_nets
#'
#'@docType data
#'@name multi_net
#'@usage c_net
#'@format three co_nets
#'\describe{
#' \item{tbls}{contians otutable rawdata}
#' \item{metas}{contians metadata}
#' \item{otus}{contians taxomomy table}
#'}
#'
#'
#'
NULL
if(F){
  data("otutab")
  t(otutab) -> totu
  metadata[,3:10] -> env
  c_net_cal(totu) -> corr
  c_net_build(corr,r_thres=0.65) -> co_net
  c_net_cal(totu,env) -> corr2
  c_net_build(corr2) -> co_net2

}
