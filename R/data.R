#'@title MetaNet networks
#'@description MetaNet co_nets
#'
#'@docType data
#'@name co_net
NULL

#'@title MetaNet networks
#'@description MetaNet co_nets
#'
#'@docType data
#'@name co_net2
NULL

#'@title MetaNet networks
#'@description MetaNet co_nets
#'
#'@docType data
#'@name co_net_rmt
NULL

#'@title Edgelist
#'@description Edgelist for `c_net_from_edgelist()`
#'
#'@docType data
#'@name arc_count
#'
NULL

#'@title Edgelist
#'@description Edgelist for `c_net_from_edgelist()`
#'
#'@docType data
#'@name arc_taxonomy
#'
NULL

#'@title MetaNet networks metadata
#'@description MetaNet co_nets
#'
#'@docType data
#'@name metab_g
NULL

#'@title MetaNet networks metadata
#'@description MetaNet co_nets
#'
#'@docType data
#'@name transc_g
NULL

#'@title MetaNet networks metadata
#'@description MetaNet co_nets
#'
#'@docType data
#'@name micro_g
NULL

#'@title MetaNet networks abundance
#'@description MetaNet co_nets
#'
#'@docType data
#'@name metab
NULL

#'@title MetaNet networks abundance
#'@description MetaNet co_nets
#'
#'@docType data
#'@name micro
NULL

#'@title MetaNet networks abundance
#'@description MetaNet co_nets
#'
#'@docType data
#'@name transc
NULL

#'@title MetaNet networks
#'@description MetaNet co_nets
#'
#'@docType data
#'@name multi1
NULL

if(F){
  # data("otutab",package = "pcutils")
  t(otutab) -> totu
  metadata[,3:10] -> env
  c_net_cal(totu) -> corr
  c_net_build(corr,r_thres=0.65) -> co_net
  c_net_cal(totu,env) -> corr2
  c_net_build(corr2) -> co_net2
}
