#========3.layout========
#' Layout coordinates
#'
#' @param go igraph
#' @param niter iterations number
#'
#' @export
#'
c_net_lay<-function(go,niter = 500){
  coors<-igraph::layout_with_fr(go,niter = niter,grid = "nogrid")
  coors<-data.frame(name = V(go)$name,X=coors[,1],Y=coors[,2])
  return(coors)
}

#' Layout with group
#'
#' @param go igraph object
#' @param group group name (default:module)
#' @param zoom1 big network layout size
#' @param zoom2 average sub_network layout size
#' @param layout1 layout1 method, one of
#' {1.numeric: use adjusted circle}
#' {2.a dataframe: rowname is group, two columns are X and Y}
#' {3.function: \code{\link[igraph]{layout_}} default: in_circle()}
#' @param layout2 one of \code{\link[igraph]{layout_}}
#'
#' @return coors
#' @export
#'
#' @examples
#' data("c_net")
#' modu_dect(co_net,method="cluster_fast_greedy") -> co_net_modu
#' g_lay(co_net_modu,group ="module",zoom1=30,zoom2 = 5,layout2 =nicely())->oridata
#' #g_lay(co_net_modu,group ="module",zoom1=30,zoom2 = 5,layout1=data.frame(X=c(1,1,1,1,2,2,2,2,3,3,3,3,3.5),Y=c(1,2,3,4,1,2,3,4,1,2,3,4,2)),layout2 =nicely())->oridata
#' modu_plot(co_net_modu,coors = oridata)
#'
#' g_lay_nice(co_net_modu,group ="module")->oridata
#' modu_plot(co_net_modu,coors = oridata)
#'
g_lay<-function(go,group="module",zoom1=20,zoom2=3,layout1=in_circle(),layout2=in_circle()){
  stopifnot(is_igraph(go))
  if(!group%in%vertex_attr_names(go))stop("no group named ",group," !")
  igraph::as_data_frame(go,what = "vertices")%>%dplyr::select(name,!!group)->nodeGroup
  colnames(nodeGroup)=c("ID","group")
  nodeGroup$group<-as.factor(nodeGroup$group)

  xs = as.data.frame(table(nodeGroup$group))
  r = xs$Freq/10
  scale_f=mmscale(r,0.5,2)
  names(scale_f)=levels(nodeGroup$group)

  big_lay<-\(zoom1,layout1){
    if(is.data.frame(layout1))da=layout1
    else if(is.numeric(layout1)){
      #弧度修正
      arg=cumsum(scale_f/sum(scale_f)*360)
      x = rep(0, length(r))
      y = rep(0, length(r))
      for (i in 1:length(r)) {
        x[i] =  sin(arg[i] * 3.14/180)
        y[i] =  cos(arg[i] * 3.14/180)
      }
      da = data.frame(x = x, y = y,row.names = levels(nodeGroup$group))
    }
    else {da = data.frame(layout_(make_ring(length(r)),layout1),
                          row.names = levels(nodeGroup$group))
    }
    #group的中心点
    da[,1]=mmscale(da[,1],-zoom1,zoom1)
    da[,2]=mmscale(da[,2],-zoom1,zoom1)
    return(da)
  }

  da = big_lay(zoom1,layout1)

  #每个group内部的分布
  oridata=data.frame()
  for (i in levels(nodeGroup$group)) {
    nodeGroup%>%dplyr::filter(group==i)%>%dplyr::pull(ID)->tmpid
    induced_subgraph(go,tmpid,impl = "copy_and_delete")->tmp_net
    igraph::layout_(tmp_net,layout2)->data

    zoom2*scale_f[i]->r2
    data[,1]=mmscale(data[,1],-r2,r2)
    data[,2]=mmscale(data[,2],-r2,r2)
    data <- data.frame(name =tmpid ,X = data[,1] + da[i, 1], Y = data[,2] +da[i, 2])

    oridata = rbind(oridata, data)
  }
  print("Big layout:")
  print(da)
  return(oridata)
}

#' Layout with group nicely
#'
#' @export
#'
#' @rdname g_lay
g_lay_nice<-function(go,group="module"){
  lib_ps("ggraph")
  stopifnot(is_igraph(go))
  if(!group%in%vertex_attr_names(go))stop("no group named ",group," !")
  igraph::as_data_frame(go,what = "vertices")%>%dplyr::select(name,!!group)->nodeGroup
  colnames(nodeGroup)=c("ID","group")
  nodeGroup$group<-as.factor(nodeGroup$group)

  edge = data.frame(from = paste("model_", nodeGroup$group,sep = ""), to = nodeGroup$ID)

  vertices_t <- data.frame(name = unique(c(as.character(edge$from),
                                           as.character(edge$to))))
  vertices_t$size = sample(1:10, nrow(vertices_t), replace = TRUE)

  mygraph <- igraph::graph_from_data_frame(edge, vertices = vertices_t)
  data = ggraph::create_layout(mygraph, layout = "circlepack",weight = size)
  coor = data %>% dplyr::filter(leaf == TRUE) %>% dplyr::select(name,x,y)
  colnames(coor) = c("name", "X", "Y")
  return(coor)
}

#========4.plot========
#' Plot a co-network
#'
#' @param go a igraph object
#' @param coors the coordinates you saved
#' @param ... additional parameters for plot.igraph()
#' @return a network plot
#' @export
#'
#' @examples
#' data("c_net")
#' c_net_plot(co_net)
#' c_net_plot(co_net2)
c_net_plot <- function(go, coors = NULL, ...) {
  lib_ps("igraph")
  # mode1，普通的单网络
  if (is.null(coors)) coors <- layout.fruchterman.reingold(go)
  else {
    coors<-coors[match(V(go)$name,coors$name),2:3]%>%as.matrix()
  }

  set.seed(123)
  vertex.attributes(go) %>% as.data.frame()-> tmp_v
  #show labels
  tmp_v%>%dplyr::top_n(10,size)%>%dplyr::pull(name)->toplabel
  tmp_v$label=ifelse(tmp_v$name%in%toplabel,tmp_v$name,NA)

  plot(go, ...,
       main = "Co-occurrence network", layout = coors,
       vertex.label.font = 1, vertex.label.color = "black",
       vertex.label.cex = 0.05 * V(go)$size,
       vertex.label=tmp_v$label,
       edge.curved = TRUE, margin = c(0, 0, 0, 0)
  )

  table(E(go)$cor>0)->numofe

  nfe<-\(x) ifelse(is.na(x),0,x)
  legend(1.2, 1, cex = 0.7, legend = c(paste0("+: ",nfe(numofe["TRUE"])), paste0("- : ",nfe(numofe["FALSE"]))),
         col = c("#48A4F0","#E85D5D"), bty = "n", title = "Correlation", lty = 1)

  f1 <- length(tmp_v$group %>% unique()) > 1
  f2 <- length(tmp_v$class %>% unique()) > 1
  f3 <- length(E(go)$class %>% unique()) > 1

  if (!f1 && f2) {
    legend(-2, 1,
           cex = 0.7, legend = tmp_v$class %>% unique(),
           col = "black", pt.bg = tmp_v$color %>% unique(), bty = "n", pch = 21
    )
  }
  # mode2，两个表交互节点
  if (f1 & !f3) {
    legend(-2, 1,
           cex = 0.7, legend = (dplyr::filter(tmp_v, group == "group1"))$class %>% unique(),
           col = "black", pt.bg = (dplyr::filter(tmp_v, group == "group1"))$color %>% unique(), bty = "n", pch = 21
    )
    legend(-2, -.5,
           cex = 0.7, legend = (dplyr::filter(tmp_v, group == "group2"))$class %>% unique(),
           col = "black", pt.bg = (dplyr::filter(tmp_v, group == "group2"))$color %>% unique(), bty = "n", pch = 22
    )
  }
  # mode3，两个网络本身加上网络间交互
  if (f1 & f3) {
    legend(-2, 1,
           cex = 0.7, legend = (dplyr::filter(tmp_v, group == "group1"))$class %>% unique(),
           col = "black", pt.bg = (dplyr::filter(tmp_v, group == "group1"))$color %>% unique(), bty = "n", pch = 21
    )
    legend(-2, -0.5,
           cex = 0.7, legend = (dplyr::filter(tmp_v, group == "group2"))$class %>% unique(),
           col = "black", pt.bg = (dplyr::filter(tmp_v, group == "group2"))$color %>% unique(), bty = "n", pch = 22
    )
    table(E(go)$class=="intra")->numofec
    legend(1.2, .7, cex = 0.7,legend = c(paste0("intra: ",nfe(numofec["TRUE"])), paste0("inter: ",nfe(numofec["FALSE"]))),col = "black", bty = "n", title = "", lty = c(1, 5))
  }
}


#' Output a graphml file for Gephi
#' @param go igraph
#' @param file graphml file name
#'
#' @export
#'
#' @examples
#' output_gephi(co_net,"test")
output_gephi<-function(go,file){
  igraph::write.graph(go,paste0(file,".graphml"),format = "graphml")
}

#' Input a graphml file exported by Gephi
#'
#' @param file graphml file exported by Gephi
#'
#' @return list contains the igraph object and coordinates
#' @export
#' @examples
#' input_gephi("test.graphml")->gephi
#' c_net_plot(gephi$go,gephi$coors)
input_gephi<-function(file){
  igraph::read.graph(file,format = "graphml")->gephi
  vertex_attr(gephi)%>%as.data.frame()->tmp_v
  #extract coors
  if(!all(c("x","y","r","g","b","id")%in%colnames(tmp_v)))stop("The file is not exported by Gephi, please use igraph::read.graph()")
  coors<-tmp_v[,c("x","y")]
  coors<-data.frame(name = tmp_v$name,X=coors[,1],Y=coors[,2])
  coors%>%mutate(X=mmscale(X,-40,40),Y=mmscale(Y,-40,40))->coors
  #transfrom color
  rgb2code(tmp_v[,c("r","g","b")])%>%dplyr::pull(code)->tmp_v$color
  E(gephi)$color=ifelse(E(gephi)$cor > 0, "#48A4F0", "#E85D5D")
  #scale size
  tmp_v$size=mmscale(tmp_v$size,1,5)
  E(gephi)$width=mmscale(E(gephi)$width,0.05,0.2)
  #delete
  tmp_v%>%dplyr::select(-c("label","x","y","r","g","b","id"))%>%as.list()->vertex.attributes(gephi)
  edge.attributes(gephi)["Edge Label"]=edge.attributes(gephi)["id"]=NULL
  return(list(go=gephi,coors=coors))
}

#' Transfer a igraph object to a ggig
#'
#' @param go igraph
#' @param net_par_res net_par() result
#' @param coors coordinates for nodes,columns: name, X, Y
#'
#' @return ggig object
#' @export
#'
#' @examples
#' coors=layout_with_fr(co_net, niter=999,grid="nogrid")
#' coors<-data.frame(name = V(co_net)$name,X=coors[,1],Y=coors[,2])
#' to.ggig(co_net,coors=coors)->ggig
#' plot(ggig)
to.ggig<-function(go,net_par_res=NULL,coors=NULL){
  if(is.null(net_par_res))net_par(go)->net_par_res
  if(is.null(coors)){
    layout.fruchterman.reingold(go)->coors
    coors<-data.frame(name = V(go)$name,X=coors[,1],Y=coors[,2])
  }

  net_par_res$v_index%<>%dplyr::left_join(.,coors,by="name",suffix = c("", ".1"))
  net_par_res$e_index%<>%dplyr::left_join(.,coors,by=c("from"="name"))%>%rename(X1=X,Y1=Y)%>%
    dplyr::left_join(.,coors,by=c("to"="name"))%>%rename(X2=X,Y2=Y)

  #show labels
  net_par_res$v_index%>%dplyr::top_n(10,size)%>%dplyr::pull(name)->toplabel
  net_par_res$v_index$label=ifelse(net_par_res$v_index$name%in%toplabel,net_par_res$v_index$name,NA)
  class(net_par_res)<-c("ggig","list")
  return(net_par_res)
}

#' Plot a ggig
#'
#' @param ggig ggig object
#'
#' @return ggplot
#' @exportS3Method
#'
plot.ggig<-function(ggig){
  lib_ps("ggplot2","ggnewscale")
  cols=unique(ggig$v_index$color)
  names(cols)=unique(ggig$v_index$class)

  p <- ggplot() +
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = inter, linewidth = width),data = ggig$e_index,alpha=0.5) +
    scale_linewidth(range = c(0.2,1.5))+scale_color_manual(values =c(positive="#48A4F0",negative="#E85D5D") )+
    ggnewscale::new_scale_color()+
    geom_point(aes(X, Y,fill = class,size = size),pch = 21, data = ggig$v_index) +
    scale_fill_manual(values = ggig$v_index$color %>% unique())+
    ggnewscale::new_scale("size")+
    geom_text(aes(X, Y,col = class,size = size,label=label), data = ggig$v_index,show.legend = F)+
    scale_color_manual(values = ggig$v_index$color %>% unique())+
    scale_size(range = c(1,4))+
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    coord_fixed(ratio = 1)+
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  p
}

