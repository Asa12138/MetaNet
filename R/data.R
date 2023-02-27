#'@title testdata (pc_otu class) for pc_tax package.
#'@description MetaNet co_nets.
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

if(F){
  data("otutab")
  t(otutab) -> totu
  metadata[,3:10] -> env
  c_net_cal(totu) -> corr
  c_net_build(corr,r_thres=0.65) -> co_net
  c_net_cal(totu,env) -> corr2
  c_net_build(corr2) -> co_net2

}
