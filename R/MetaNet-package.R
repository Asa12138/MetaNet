#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import ggplot2
#' @import igraph
#' @importFrom dplyr select filter left_join arrange count mutate group_by summarise pull
#' @importFrom igraph V E
#' @importFrom utils data combn  head tail
#' @importFrom stats aggregate median var sd setNames runif relevel coef fitted cor coefficients
#' @importFrom stats time na.omit kmeans p.adjust density approx ks.test smooth.spline
#' @importFrom pcutils lib_ps mmscale get_cols trans guolv hebing update_param
#' @importFrom graphics legend
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom reshape2 acast melt
## usethis namespace: end
NULL

if(F){
  # data("multi_net")
  # imgurl="inst/figures/multi.png"
  # png(imgurl,width = 1000,height = 800,res =140,bg = NA)
  # par(mar = rep(0,4) + 0.1)
  # modu_net(3,24)%>%c_net_plot(mark_module = T,labels_num = 0,edge.color=c("#5fa8d3","#ff4d6d"),
  #                         legend = F,main="",edge.width=2.5,vertex.frame.width=1.5)
  # #c_net_plot(multi1,labels_num = 0,vertex_size_range = c(3,8),legend = F,main="")
  # dev.off()
  # library(showtext)
  # font_families()
  # #font_add("MarkerFelt","/System/Library/Fonts/MarkerFelt.ttc")
  # font_add("BoB","/Users/asa/Library/Fonts/SmileySans-Oblique.ttf")
  # imgurl="inst/figures/multi.png"
  # library(hexSticker)
  # s=sticker(imgurl,dpi=600,asp = 0.9,
  #           s_x=1, s_y=.7, s_width=.7,s_height = .7,
  #           package="MetaNet", p_size=48,p_color = "#2f3e46",p_y = 1.45,
  #           p_family ="BoB",p_fontface = "plain",
  #           h_fill = "#d8e2dc",h_size = 2,h_color = "#84a98c",
  #           filename="inst/figures/hexSticker1.png")
  # s
}

#' Show MetaNet logo
#'
#' @return picture
#' @export
show_MetaNet_logo=function(){
  pcutils::read.file(system.file("figures/hexSticker1.png",package = "MetaNet"))
}
