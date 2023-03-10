#========3.layout========
#' Layout coordinates
#'
#' @param go igraph
#' @param method as_star(), as_tree(), in_circle(), nicely(), on_grid(), on_sphere(),randomly(), with_dh(), with_fr(), with_gem(), with_graphopt(), with_kk(),with_lgl(), with_mds(),as_line(), as_arc(), as_polygon(), as_polyarc(). see \code{\link[igarph]{layout_}}.
#' @param order_by order nodes according to a node attribute
#' @param order_ls manual the discrete variable with a vector, or continuous variable with "desc" to decreasing
#' @param seed random seed
#' @param line_curved consider line curved, only for some layout methods like as_line(), as_polygon().default:0
#'
#' @return coordinates for nodes,columns: name, X, Y
#' @export
#' @examples
#' c_net_lay(co_net)->coors
#' c_net_plot(co_net,coors)
#' c_net_plot(co_net,c_net_lay(co_net,in_circle()),vertex.size=2)
#' c_net_plot(co_net,c_net_lay(co_net,in_circle(),order_by="v_class"),vertex.size=2)
#' c_net_plot(co_net,c_net_lay(co_net,in_circle(),order_by="size",order_ls="desc"))
#' c_net_plot(co_net,c_net_lay(co_net,as_polygon(3)))
c_net_lay<-function(go,method=igraph::nicely(),order_by=NULL,order_ls=NULL,seed=1234,line_curved=1){
  set.seed(seed)
  if("igraph_layout_spec"%in%class(method))coors<-igraph::layout_(go,method)
  else if ("poly"%in%class(method)) {coors=method(go,group2 =order_by,group2_order =order_ls)}
  else if ("layout"%in%class(method)) coors=method(go)
  else stop("No valid method")

  #order
  if(is.matrix(coors)){
    if(is.null(order_by))coors<-data.frame(name =V(go)$name,X=coors[,1],Y=coors[,2])
    else {
      get_v(go)->tmp_v
      ordervec=tmp_v[,order_by]
      if(is.numeric(ordervec)){
        name=tmp_v[order(ordervec,decreasing = is.null(order_ls)),"name"]}
      else {
        ordervec=change_fac_lev(ordervec,order_ls)
        name=tmp_v[order(ordervec),"name"]
      }
      coors<-data.frame(name =name,X=coors[,1],Y=coors[,2])
    }}

  #if line type, need to consider edge.curved

  if(line_curved){
    if("line"%in%class(method)){
    tmp_e=data.frame(igraph::as_data_frame(go))[,c("from","to")]
    if(nrow(tmp_e)>0)curved=data.frame(tmp_e,curved=line_curved)
    else curved=NULL
    coors=list(coors=coors,curved=curved)
    }
  }

  return(coors)
}

#' @export
is_lay=\(x){
  any(class(x)%in%c("igraph_layout_spec","layout"))
}

#' @export
get_coors=\(coors,go,...){
  edge_curved=NULL
  if(is.null(coors)) coors <- igraph::nicely()

  if(is_lay(coors))coors=c_net_lay(go,coors,...)
  if(all(class(coors)=="list")){
    if(!is.null(coors$curved)){edge_curved=dplyr::left_join(get_e(go),coors$curved,
                                 by=c("from","to"),suffix = c(".x", ""))%>%
      dplyr::select("from","to","curved")}
    if(!is.null(coors$coors))coors<-coors$coors
    else coors=c_net_lay(go,igraph::nicely(),...)
  }
  if(is.matrix(coors))coors=data.frame(name=V(go)$name,X=coors[,1],Y=coors[,2])
  if(is.data.frame(coors)) {
    coors<-coors[match(V(go)$name,coors$name),]
    return(list(coors=coors,curved=edge_curved))
  }
  return("coors wrong")
}

#' Layout as a line
#'
#' @param angle anticlockwise rotation angle
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#'
#' @examples
#' as_line()(co_net)
as_line=\(angle=0){fun=\(go){
  nv=length(V(go))
  data.frame(x=seq(-cos(angle),cos(angle),len=nv),
             y=seq(-sin(angle),sin(angle),len=nv))%>%as.matrix()%>%round(.,4)}
class(fun)=c("line","layout","function")
fun
}

#' Layout as a arc
#'
#' @param angle anticlockwise rotation angle
#' @param arc the radian of arc
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#'
#' @examples
#' as_arc()(co_net)
as_arc=\(angle=0,arc=pi){fun=\(go){
  #(0,0) is the midpoint of circle
  nv=length(V(go))
  theta=seq(-arc/2+angle,arc/2+angle,len=nv)
  coor=data.frame(x=cos(theta),y=sin(theta))
  as.matrix(coor)%>%round(.,4)}
class(fun)=c("layout","function")
fun
}

#deprecated,(0,0) in on arc
if(F){
  as_arc1=\(angle=0,arc=pi){\(go){
    nv=length(V(go))
    theta=seq(-arc/2,arc/2,len=nv)
    coor=data.frame(x=sin(theta+angle)-sin(angle),
                    y=cos(angle)-cos(theta+angle))
    as.matrix(coor)%>%round(.,4)}
  }
}

#' Layout as a polygon
#'
#' @param n how many edges of this polygon
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#'
#' @examples
#' as_polygon()(co_net)
as_polygon=\(n=3,line_curved=0.5){fun=\(go,group2 =NULL,group2_order =NULL){
  V(go)$poly_group=rep(paste0("omic",1:n),len=length(go))
  if(n<2)stop("n should bigger than 1")
  g_lay_polygon(go,group ="poly_group",group2 =group2 ,group2_order = group2_order,line_curved=line_curved)->oridata
  oridata
}
class(fun)=c("poly","layout","function")
fun
}

#' Layout as a poly_arc
#'
#' @param n how many arcs of this poly_arc
#' @param space the space between each arc, default: pi/3
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#'
#' @examples
#' as_polyarc()(co_net)
as_polyarc=\(n=3,space = pi/3){fun=\(go,group2 =NULL,group2_order =NULL){
  V(go)$poly_group=rep(paste0("omic",1:n),len=length(go))
  if(n<2)stop("n should bigger than 1")
  g_lay_polyarc(go,"poly_group",space = space,group2 =group2 ,group2_order = group2_order)->oridata
  oridata
}
class(fun)=c("poly","layout","function")
fun
}

#' Interactive layout
#'
#' @examples
#' make_graph("Zachary")%>%manual_lay
manual_lay<-function(go,canvas.width =600,canvas.height=400){
  flag="y"
  if(length(V(go))>100){
    print("Too big network, recommend using Gephi to layout,still use tkplot?")
    flag=readline("yes/no(y/n):")}
  if(tolower(flag)%in%c("yes","y")){
    x <- igraph::tkplot(go,
                vertex.size = 10,
                canvas.width = canvas.width,
                canvas.height = canvas.height
    )
    #Here: Move nodes within the tkplot window to a layout you like!
    coors <- tkplot.getcoords(x)
  }
  else print("Quit.");return()
  coors<-data.frame(name = V(go)$name,X=coors[,1],Y=coors[,2])
  coors
}

#' Layout with group
#'
#' @param go igraph object
#' @param group group name (default:module)
#' @param zoom1 big network layout size
#' @param zoom2 average sub_network layout size, or numeric vector, or "auto"
#' @param layout1 layout1 method, one of
#' {1.a dataframe or matrix: rowname is group, two columns are X and Y}
#' {2.function: layout method for \code{\link{c_net_lay}} default: in_circle()}
#' @param layout2 one of functions: layout method for \code{\link{c_net_lay}}, or a list of functions.
#'
#' @return coors
#' @export
#'
#' @examples
#' data("c_net")
#' modu_dect(co_net,method="cluster_fast_greedy") -> co_net_modu
#' g_lay(co_net_modu,group ="module",zoom1=30,zoom2 = 1:5,layout2 =as_line())->oridata
#' #g_lay(co_net_modu,group ="module",zoom1=30,zoom2 = 5,layout1=data.frame(X=c(1,1,1,1,2,2,2,2,3,3,3,3,3.5),Y=c(1,2,3,4,1,2,3,4,1,2,3,4,2)),layout2 =nicely())->oridata
#' plot(co_net_modu,coors = oridata)
#'
#' g_lay_nice(co_net_modu,group ="module")->oridata
#' plot(co_net_modu,coors = oridata)
#'
g_lay<-function(go,group="module",layout1=in_circle(),zoom1=20,layout2=in_circle(),zoom2=3,show_big_lay=F,...){
  stopifnot(is_igraph(go))
  if(!group%in%vertex_attr_names(go))stop("no group named ",group," !")
  get_v(go)%>%dplyr::select(name,!!group)->nodeGroup
  colnames(nodeGroup)=c("ID","group")
  nodeGroup$group<-as.factor(nodeGroup$group)

  big_lay<-\(zoom1,layout1){
    c_net_update(make_ring(nlevels(nodeGroup$group)))->tmp_da_net
    da=get_coors(layout1,tmp_da_net)[["coors"]]
    da$name=NULL
    if(!all(levels(nodeGroup$group)%in%rownames(da)))rownames(da)=levels(nodeGroup$group)

    #center of each group
    scale_f=(ceiling(max(da)-min(da)))
    scale_f=ifelse(scale_f==0,1,scale_f/2)
    da=da/scale_f*zoom1
    colnames(da)=c("X","Y")
    return(da)
  }
  da = big_lay(zoom1,layout1)

  #layout of vertexes in one group
  layoutls=list()
  if(is_lay(layout2))layoutls=rep(list(layout2),nlevels(nodeGroup$group))
  else if(all(class(layout2)=="list")){
    if(is.null(names(layout2)))layoutls=rep(layout2,len=nlevels(nodeGroup$group))
    else {
      for (i in levels(nodeGroup$group)) {
        if(i%in%names(layout2))layoutls[[i]]=layout2[[i]]
        else {layoutls[[i]]=in_circle();warning("layout of ",i," not set, use in_circle()")}
      }
    }
  }
  if(is.null(names(layoutls)))names(layoutls)=levels(nodeGroup$group)

  zoom2=rep(zoom2,nlevels(nodeGroup$group))
  if(zoom2[1]=="auto")zoom2=ceiling((table(nodeGroup$group))^(1/3))
  names(zoom2)=levels(nodeGroup$group)

  #get coors
  oridata=list();oridata2=list()
  for (i in levels(nodeGroup$group)) {
    nodeGroup%>%dplyr::filter(group==i)%>%dplyr::pull(ID)->tmpid
    igraph::subgraph(go,tmpid)->tmp_net

    if(T){
      get_coors(layoutls[[i]],tmp_net,...)->coors
      data=coors$coors
      if("igraph_layout_spec"%in%class(layoutls[[i]]))data[,c("X","Y")]=igraph::norm_coords(as.matrix(data[,c("X","Y")]))
      data[,"X"]=data[,"X"]*zoom2[i]+ da[i, "X"]
      data[,"Y"]=data[,"Y"]*zoom2[i]+ da[i, "Y"]
    }
    oridata[[i]] = data
    oridata2[[i]]= coors$curved
  }
  oridata=do.call(rbind,oridata)
  oridata2=do.call(rbind,oridata2)

  if(show_big_lay){
    print("Big layout:")
    print(da)
    plot(da,pch=21,bty="n",bg=get_cols(nrow(da),"col2"),main="Big layout coordinates")
    text(da$X*0.8,da$Y*0.9,rownames(da))}

  return(list(coors=oridata,curved=oridata2))
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
  get_v(go)%>%dplyr::select(name,!!group)->nodeGroup
  colnames(nodeGroup)=c("ID","group")
  nodeGroup$group<-as.factor(nodeGroup$group)

  edge = data.frame(from = paste("group_", nodeGroup$group,sep = ""), to = nodeGroup$ID)

  vertices_t <- data.frame(name = unique(c(as.character(edge$from),
                                           as.character(edge$to))))
  vertices_t$size = sample(1:10, nrow(vertices_t), replace = TRUE)

  mygraph <- igraph::graph_from_data_frame(edge, vertices = vertices_t)
  data = ggraph::create_layout(mygraph, layout = "circlepack",weight = size)
  coor = data %>% dplyr::filter(leaf == TRUE) %>% dplyr::select(name,x,y)
  colnames(coor) = c("name", "X", "Y")
  return(coor)
}

#' @export
plot_gg_circle<-function(go=co_net,group="v_class",sep=0){
  #???????????????
  lib_ps("ggraph")
  stopifnot(is_igraph(go))
  #construct hierarchy or just give it.
  {
    if(!group%in%vertex_attr_names(go))stop("no group named ",group," !")
    get_v(go)%>%dplyr::select(name,v_group,!!group)->nodeGroup
    colnames(nodeGroup)=c("ID","v_group","group")

    nodeGroup$v_group<-as.factor(nodeGroup$v_group)
    nodeGroup$group<-as.factor(nodeGroup$group)

    hierarchy=rbind(data.frame(from ="origin",to=levels(nodeGroup$v_group)),
                    data.frame(distinct(nodeGroup,v_group,group)%>%rename(from="v_group",to="group")),
                    data.frame(from =nodeGroup$group, to = nodeGroup$ID))
    if(sep){
      #????????????
      hierarchy=rbind(data.frame(from ="origin",to=paste0(levels(nodeGroup$v_group),"_a")),
                      data.frame(distinct(nodeGroup,v_group,group)%>%rename(from="v_group",to="group")),
                      data.frame(from =nodeGroup$group, to = nodeGroup$ID))
      hierarchy=rbind(hierarchy,
                      data.frame(from =rep(paste0(levels(nodeGroup$v_group),"_a"),each=2),
                                 to=paste0(rep(levels(nodeGroup$v_group),each=2),c("","_b"))),
                      data.frame(from=rep(paste0(levels(nodeGroup$v_group),"_b"),each=sep),
                                 to=paste0(rep(paste0(levels(nodeGroup$v_group),"_bbb"),each=sep),seq_len(sep))))
    }
  }
  edge = get_e(go)
  # create a vertices data.frame. One line per object of our hierarchy
  vertices  <-  data.frame(
    name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to)))
  )
  # Let's add a column with the group of each name. It will be useful later to color points
  vertices$gg_group  <-  hierarchy$from[ match( vertices$name, hierarchy$to )]
  vertices=left_join(vertices,get_v(go))

  # Create a graph object
  mygraph <- graph_from_data_frame(hierarchy, vertices=vertices )

  # The connection object must refer to the ids of the leaves:
  edge$from  <-  match( edge$from, vertices$name)
  edge$to  <-  match( edge$to, vertices$name)

  # Basic usual argument
  p=ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
    geom_conn_bundle(data = get_con(from = edge$from, to = edge$to,e_type=edge$e_type,width=edge$width),
                      alpha=0.2, aes(width=width,colour=e_type),tension=0.5) +
    scale_edge_width(range = c(0.2,1))+
    #scale_edge_colour_discrete()+
    theme_void()

  #?????????????????????
  if(F){
    offset=1.05
    for(i in c("Phylum","Kingdom")){
      p=p+geom_node_point(aes(filter = leaf, x = x*offset, y=y*offset, colour=get(i),shape=v_group),size=3)+
        ggnewscale::new_scale_color()
      offset=offset+0.1
    }
  }
  p+geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=v_class,shape=v_group),size=3) +
    scale_colour_manual(values= pcutils::get_cols(14,"col3"),name=group,na.translate=F)+
    scale_shape(na.translate=F)+guides(shape=guide_none())
}

#' Layout with group as a polygon
#'
#' @param go igraph
#' @param group group name (default:v_group)
#' @param group2 group2 name, will order nodes in each group according to group2_order (default:v_class)
#' @param group2_order group2_order
#'
#' @return coors
#' @export
#'
#' @examples
#' g_lay_polygon(multi1)->oridata
#' c_net_plot(multi1,oridata)
#' g_lay_polyarc(multi1,group2_order = c(LETTERS[4:1]))->oridata
#' c_net_plot(multi1,oridata)
g_lay_polygon<-function(go,group="v_group",group2=NULL,group2_order=NULL,line_curved=0.5){
  n=length(unique(igraph::vertex.attributes(go)[[group]]))

  if(n<2)stop("n should bigger than 1")
  angle_ls=-pi/2+(0:(n-1))*2*pi/n
  fun_ls=sapply(angle_ls, \(i)as_line(i))

  g_lay(go,group =group,zoom1=1,zoom2 =0.9*(ifelse(n>2,tan(pi/n),2)),
        layout2 =fun_ls,order_by=group2,order_ls=group2_order)->oridata

  oridata$curved$curved=line_curved
  oridata
}

#' Layout with group as a polyarc
#'
#' @rdname g_lay_polygon
#' @export
g_lay_polyarc<-function(go,group="v_group",group2=NULL,group2_order=NULL,space=pi/4){
  get_v(go)->tmp_v
  group1=as.factor(tmp_v[,group])
  n=nlevels(group1)
  if(n<2)stop("n should bigger than 1")
  #consider each group numbers!!!
  g_num=table(group1)
  sep=space/n
  arc_r=(2*pi-sep*n)*g_num/length(group1)

  #coordinate
  coors=data.frame()
  theta1=0
  for (i in names(arc_r)) {
    tmp_t=seq(theta1,theta1+arc_r[i],len=g_num[i])
    tmp_v1=tmp_v[tmp_v[,group]==i,]
    {if(!is.null(group2)){
      ordervec=tmp_v1[,group2]
      if(is.numeric(ordervec)){
        name=tmp_v1[order(ordervec,decreasing = is.null(group2_order)),"name"]}
      else {
        ordervec=change_fac_lev(ordervec,group2_order)
        name=tmp_v1[order(ordervec),"name"]
      }
    }
    else name=tmp_v1$name}
    tmp_coor=data.frame(name=name,X=cos(tmp_t),Y=sin(tmp_t))
    coors=rbind(coors,tmp_coor)
    theta1=theta1+arc_r[i]+sep
  }
  coors
}

#========4.plot========

#' Plot a metanet
#'
#' @param go metanet object
#'
#' @exportS3Method
plot.metanet=function(go,...){
  if(is.null(get_n(go)$n_type)) c_net_plot(go,...)
  else if(get_n(go)$n_type=="skeleton")skeleton_plot(go,...)
  else if(get_n(go)$n_type=="module")c_net_plot(go,group_legend_title = "Module",...)
  else c_net_plot(go,...)
}

#' Plot a network
#'
#' @param go a igraph object
#' @param coors the coordinates you saved
#' @param ... additional parameters for \code{\link[igraph]{igraph.plotting}}
#' @param labels_num show how many labels,>1 indicates number, <1 indicates fraction ,default:5
#' @param legend_number legend with numbers
#' @param legend all legends
#' @param legend_position legend_position, default: c(left_leg_x=-1.9,left_leg_y=1,right_leg_x=1.2,right_leg_y=1)
#' @param legend_cex 	character expansion factor relative to current par("cex"), default: 1
#' @param lty_legend logical
#' @param size_legend logical
#' @param edge_legend logical
#' @param col_legend logical
#' @param width_legend logical
#' @param lty_legend_title lty_legend_title
#' @param size_legend_title size_legend_title
#' @param edge_legend_title edge_legend_title
#' @param edge_legend_order edge_legend_order vector, e.g. c("positive","negative")
#' @param width_legend_title width_legend_title
#' @param col_legend_order col_legend_order vector,
#' @param group_legend_title group_legend_title, length must same to the numbers of v_group
#' @param group_legend_order group_legend_order vector
#' @param mark_module logical, mark the modules?
#' @param mark_color mark colors
#' @param seed random seed, default:1234, make sure each plot is the same.
#' @param plot_module logical, plot module?
#'
#' @return a network plot
#' @export
#'
#' @examples
#' data("c_net")
#' c_net_plot(co_net)
#' c_net_plot(co_net_rmt)
#' c_net_plot(co_net2)
#' c_net_plot(multi1)
c_net_plot <- function(go, coors = NULL,...,labels_num=5,
                       legend_number=F,legend=T,legend_cex=1,
                       legend_position=c(left_leg_x=-1.9,left_leg_y=1,right_leg_x=1.2,right_leg_y=1),
                       lty_legend=F,lty_legend_title="Edge class",
                       size_legend=F,size_legend_title="Node Size",
                       edge_legend=T,edge_legend_title="Edge type",edge_legend_order=NULL,
                       width_legend=F,width_legend_title="Edge width",
                       col_legend=T,col_legend_order=NULL,
                       group_legend_title=NULL,group_legend_order=NULL,
                       plot_module=F,mark_module=F,mark_color=NULL,seed=1234){
  lib_ps("igraph")
  set.seed(seed)

  #modules
  if(plot_module){
    go=to_module_net(go)
    if(is.null(group_legend_title))group_legend_title="Module"
  }

  #get network type
  main="Network"
  if(!is.null(get_n(go)$n_type)){
    switch (get_n(go)$n_type,
            "single" = {main="Correlation network"},
            "bipartite" = {main="Bipartite network"},
            "multi_full" ={main="Multi-omics network"},
            "module"={main=paste0(get_n(go)$n_modules,"-modules network");plot_module=T},
            "skeleton" ={main=paste0(get_n(go)$skeleton," skeleton network")},
            default={main="Network"}
    )
  }

  get_v(go)-> tmp_v
  get_e(go)->tmp_e

  #get coordinates
  coors<-get_coors(coors,go,seed=seed)
  if(is.null(coors$curved))edge_curved=0
  else edge_curved=coors$curved$curved
  edge_curved[is.na(edge_curved)]=0

  coors=coors$coors[,c("X","Y")]%>%as.matrix()

  #scale the size and width
{
  node_size_text1=c()
  node_size_text2=c()
  for(i in unique(tmp_v$v_group)){
    node_size_text1=c(node_size_text1,min(tmp_v[tmp_v$v_group==i,"size"]))
    node_size_text2=c(node_size_text2,max(tmp_v[tmp_v$v_group==i,"size"]))
    tmp_v[tmp_v$v_group==i,"size"]%<>%pcutils::mmscale(.,3,12)
  }
  tmp_e$width=pcutils::mmscale(tmp_e$width,0.5,1)

  #some custom parameters
  params=list(...)
  params_name=names(params)
  if("vertex.size"%in%params_name)tmp_v$size=params[["vertex.size"]]
  if("vertex.color"%in%params_name)tmp_v$color=condance(data.frame(tmp_v$color,tidai(tmp_v$v_class,params[["vertex.color"]])))
  if("vertex.shape"%in%params_name)tmp_v$shape=condance(data.frame(tmp_v$shape,tidai(tmp_v$v_group,params[["vertex.shape"]])))
  if("edge.color"%in%params_name)tmp_e$color=condance(data.frame(tmp_e$color,tidai(tmp_e$e_type,params[["edge.color"]])))
  if("edge.lty"%in%params_name)tmp_e$lty=condance(data.frame(tmp_e$lty,tidai(tmp_e$e_class,params[["edge.lty"]])))
  if("edge.width"%in%params_name)tmp_e$width=params[["edge.width"]]
  if("vertex.label"%in%params_name)tmp_e$label=params[["vertex.label"]]
}
  #show labels
  {if(labels_num>=1){tmp_v%>%dplyr::top_n(labels_num,size)%>%dplyr::pull(name)%>%head(labels_num)->toplabel}
  else {tmp_v%>%dplyr::top_frac(labels_num,size)%>%dplyr::pull(name)%>%head(ceiling(labels_num*nrow(tmp_v)))->toplabel}
  tmp_v$label=ifelse(tmp_v$name%in%toplabel,tmp_v$label,NA)
  }

  #modules set
  {
  if(mark_module){
    new_modu=as_module(setNames(tmp_v$module,tmp_v$name))
    new_modu[["others"]]=NULL

    module_color=pcutils::get_cols(length(new_modu))
    if(!is.null(mark_color))module_color=condance(data.frame(module_color,tidai(names(new_modu),mark_color)))

    module_color=setNames(module_color,names(new_modu))
    module_color=module_color[names(module_color)!="others"]
  }
  else {
    new_modu=module_color=NULL
  }}

  #main plot
  {
    igraph::plot.igraph(go, layout = coors,
                vertex.size=tmp_v$size,
                vertex.color=tmp_v$color,
                vertex.shape=tmp_v$shape,
                edge.color=tmp_e$color,
                edge.lty=tmp_e$lty,
                edge.width=tmp_e$width,
                mark.groups = new_modu,
                mark.col = pcutils::add_alpha(module_color[names(new_modu)],0.3),
                mark.border =module_color[names(new_modu)],
                ...,
                main = main,
                vertex.label.font = 1, vertex.label.color = "black",
                vertex.label.cex = 0.08*tmp_v$size,
                vertex.label=tmp_v$label,
                edge.arrow.size=0.3,
                edge.arrow.width=0.5,
                edge.curved = edge_curved, margin = c(0, 0, 0, 0))
  }

  if(!legend)return(message(""))

  #produce legends
  legend_position_default=c(left_leg_x=-1.9,left_leg_y=1,right_leg_x=1.2,right_leg_y=1)
  if(is.null(legend_position))legend_position=legend_position_default
  else {
    if(is.null(names(legend_position)))legend_position=setNames(legend_position,names(legend_position_default)[seq_along(legend_position)])
    legend_position=condance(data.frame(legend_position_default,tidai(names(legend_position_default),legend_position)))
    names(legend_position)=names(legend_position_default)
  }
  left_leg_x=legend_position["left_leg_x"]
  left_leg_y=legend_position["left_leg_y"]
  right_leg_x=legend_position["right_leg_x"]
  right_leg_y=legend_position["right_leg_y"]

  if(size_legend){
    legend(right_leg_x, right_leg_y, cex = 0.7*legend_cex,title.font = 2,title = size_legend_title,title.adj = 0,
           legend = c(paste(node_size_text1,collapse = "/"),paste(node_size_text2,collapse = "/")),adj = 0,title.cex = 0.8*legend_cex,
           col="black", bty = "n",pch=21,pt.cex = c(min(tmp_v$size),max(tmp_v$size))*legend_cex/5)
    right_leg_y=right_leg_y-0.5*legend_cex
    }

  if(edge_legend){
    edges=change_fac_lev(tmp_e$e_type,edge_legend_order)
    edges=levels(edges)
    edge_cols=distinct(tmp_e,color,e_type)
    edge_cols=setNames(edge_cols$color,edge_cols$e_type)
    if(legend_number){
      eee=table(tmp_e$e_type)
      le_text= paste(edges,eee[edges],sep = ": ")
    }
    else le_text= edges
    legend(right_leg_x, right_leg_y, cex = 0.7*legend_cex,title.font = 2,title = edge_legend_title,title.adj = 0,
           legend = le_text,adj = 0,title.cex = 0.8*legend_cex,
           col=edge_cols[edges], bty = "n",  lty = 1)
    right_leg_y=right_leg_y-(length(unique(tmp_e$color))*0.12+0.2)*legend_cex
  }

  if(width_legend){
    legend(right_leg_x, right_leg_y, cex = 0.7*legend_cex,title.font = 2,title = width_legend_title,title.adj = 0,
           legend = c(min(E(go)$width),max(E(go)$width))%>%round(.,3),adj = 0,title.cex = 0.8*legend_cex,
           col="black", bty = "n",  lty = 1,lwd = c(min(tmp_e$width),max(tmp_e$width)))
    right_leg_y=right_leg_y-0.5*legend_cex
  }

  if(lty_legend){
    edges=levels(factor(tmp_e$e_class))
    edge_ltys=distinct(tmp_e,lty,e_class)
    edge_ltys=setNames(edge_ltys$lty,edge_ltys$e_class)

    if(legend_number){
      eee=table(tmp_e$e_class)
      le_text= paste(edges,eee[edges],sep = ": ")
    }
    else le_text= edges
    legend(right_leg_x,right_leg_y, cex = 0.7*legend_cex,title.font = 2,title = lty_legend_title,title.adj = 0,
           legend = le_text,adj = 0,title.cex = 0.8*legend_cex,
           col="black", bty = "n",  lty =edge_ltys[edges])
    }

  if(col_legend){
    pchls=c("circle"=21,"square"=22)
    vgroups=change_fac_lev(tmp_v$v_group,group_legend_order)
    vgroups=levels(vgroups)
    if(is.null(group_legend_title))group_legend_title=setNames(vgroups,vgroups)
    else if(is.null(names(group_legend_title)))group_legend_title=setNames(rep(group_legend_title,len=length(vgroups)),vgroups)

    for(g_i in vgroups){
      tmp_v1=tmp_v[tmp_v$v_group==g_i,c("v_class","color","shape")]

      vclass=change_fac_lev(tmp_v1$v_class,col_legend_order)
      vclass=levels(vclass)

      node_cols=distinct(tmp_v1,color,v_class)
      node_cols=setNames(node_cols$color,node_cols$v_class)
      node_shapes=distinct(tmp_v1,shape,v_class)
      node_shapes=setNames(node_shapes$shape,node_shapes$v_class)

      if(legend_number){
        eee=table(tmp_v1$v_class)
        le_text= paste(vclass,eee[vclass],sep = ": ")
      }
      else le_text= vclass
      legend(left_leg_x, left_leg_y, cex = 0.7*legend_cex, adj = 0,
             legend = le_text,title.cex = 0.8*legend_cex,
             title =group_legend_title[g_i],title.font = 2,title.adj = 0,
             col = "black", pt.bg = node_cols[vclass], bty = "n", pch = pchls[node_shapes[vclass]])

      left_leg_y=left_leg_y-(length(vclass)*0.12+0.2)*legend_cex
  }}
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
#' coors=c_net_lay(co_net)
#' to.ggig(co_net,coors=coors)->ggig
#' plot(ggig)
#' to.ggig(multi1,coors=c_net_lay(multi1))->ggig
#' plot(ggig,labels_num=0.3)
#'
to.ggig<-function(go,coors=NULL){
  list(n_index=get_n(go),v_index=get_v(go),e_index=get_e(go))->net_par_res

  if(is.null(coors))coors=c_net_lay(co_net)
  # coors$X=max(coors$X)-coors$X
  # coors$Y=max(coors$Y)-coors$Y

  #add coors
  coors=coors[,1:3]%>%na.omit()
  net_par_res$v_index%<>%dplyr::left_join(.,coors,by="name",suffix = c("", ".1"))
  net_par_res$e_index%<>%dplyr::left_join(.,coors,by=c("from"="name"))%>%rename(X1=X,Y1=Y)%>%
    dplyr::left_join(.,coors,by=c("to"="name"))%>%rename(X2=X,Y2=Y)

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
plot.ggig<-function(ggig,coors = NULL,...,labels_num=0,
                    legend_number=F,legend=T,
                    lty_legend=F,lty_legend_title="Edge class",
                    size_legend=F,size_legend_title="Node Size",
                    edge_legend=T,edge_legend_title="Edge type",edge_legend_order=NULL,
                    width_legend=F,width_legend_title="Edge width",
                    col_legend=T,col_legend_order=NULL,
                    class_legend_title="Node class",group_legend_order=NULL){
  lib_ps("ggplot2","ggnewscale")
  ggig$v_index->tmp_v
  ggig$e_index->tmp_e
  set.seed(1234)
  #get coordinates
  if (!is.null(coors)){
    tmp_v$X=tmp_v$Y=NULL
    tmp_e$X1=tmp_e$X2=tmp_e$Y1=tmp_e$Y2=NULL
    #add coors
    tmp_v%<>%dplyr::left_join(.,coors,by="name",suffix = c("", ".1"))
    tmp_e%<>%dplyr::left_join(.,coors,by=c("from"="name"))%>%rename(X1=X,Y1=Y)%>%
      dplyr::left_join(.,coors,by=c("to"="name"))%>%rename(X2=X,Y2=Y)
  }

  #get network type
  main="Network"
  if(!is.null(ggig$n_index$n_type)){
    switch (ggig$n_index$n_type,
            "single" = {main="Correlation network"},
            "bipartite" = {main="Bipartite network"},
            "multi_full" ={main="Multi-omics network"},
            "skeleton" ={main=paste0(ggig$n_index$skeleton," skeleton network")},
            default={main="Network"}
    )
  }
  #scale the size and width
  {
    node_size_text1=c()
    node_size_text2=c()
    for(i in unique(tmp_v$v_group)){
      node_size_text1=c(node_size_text1,min(tmp_v[tmp_v$v_group==i,"size"]))
      node_size_text2=c(node_size_text2,max(tmp_v[tmp_v$v_group==i,"size"]))
      tmp_v[tmp_v$v_group==i,"size"]%<>%mmscale(.,2,10)
    }
    node_size_text=c(paste(node_size_text1,collapse = "/"),paste(node_size_text2,collapse = "/"))
    edge_width_text=c(min(tmp_e$width),max(tmp_e$width))
    tmp_e$width=mmscale(tmp_e$width,0.5,1)
    #new shapes
    tmp_v$shape=tidai(tmp_v$v_group,21:26)

    #some custom parameters
    params=list(...)
    params_name=names(params)
    if("vertex.size"%in%params_name)tmp_v$size=params[["vertex.size"]]
    if("vertex.color"%in%params_name)tmp_v$color=condance(data.frame(tmp_v$color,tidai(tmp_v$v_class,params[["vertex.color"]])))
    if("vertex.shape"%in%params_name)tmp_v$shape=condance(data.frame(tmp_v$shape,tidai(tmp_v$v_group,params[["vertex.shape"]])))
    if("edge.color"%in%params_name)tmp_e$color=condance(data.frame(tmp_e$color,tidai(tmp_e$e_type,params[["edge.color"]])))
    if("edge.lty"%in%params_name)tmp_e$lty=condance(data.frame(tmp_e$lty,tidai(tmp_e$e_class,params[["edge.lty"]])))
    if("edge.width"%in%params_name)tmp_e$width=params[["edge.width"]]
    if("vertex.label"%in%params_name)tmp_e$label=params[["vertex.label"]]
  }
  #show labels
  if(labels_num>=1)tmp_v%>%dplyr::top_n(labels_num,size)%>%dplyr::pull(name)%>%head(labels_num)->toplabel
  else tmp_v%>%dplyr::top_frac(labels_num,size)%>%dplyr::pull(name)%>%head(ceiling(labels_num*nrow(tmp_v)))->toplabel
  tmp_v$label=ifelse(tmp_v$name%in%toplabel,tmp_v$label,NA)

  if(T){
    tmp_e$e_type=change_fac_lev(tmp_e$e_type,edge_legend_order)
    edges=levels(tmp_e$e_type)
    edge_cols=distinct(tmp_e,color,e_type)
    edge_cols=setNames(edge_cols$color,edge_cols$e_type)
    if(legend_number){
      eee=table(tmp_e$e_type)
      edge_text= paste(edges,eee[edges],sep = ": ")
    }
    else edge_text= edges}

  if(T){
    edges1=levels(factor(tmp_e$e_class))
    edge_ltys=distinct(tmp_e,lty,e_class)
    edge_ltys=setNames(edge_ltys$lty,edge_ltys$e_class)

    if(legend_number){
      eee=table(tmp_e$e_class)
      lty_text= paste(edges1,eee[edges1],sep = ": ")
    }
    else lty_text= edges1}

  if(T){
    pchls=c("circle"=21,"square"=22)

    vgroups=change_fac_lev(tmp_v$v_group,group_legend_order)
    new_f=c()
    for(g_i in levels(vgroups)){
      tmp_v1=tmp_v[tmp_v$v_group==g_i,c("v_class","color","shape")]
      tmp_f=change_fac_lev(tmp_v1$v_class,col_legend_order)
      new_f=c(new_f,levels(tmp_f))
    }

    tmp_v$v_class=change_fac_lev(tmp_v$v_class,new_f)
    vclass=levels(tmp_v$v_class)

    node_cols=distinct(tmp_v,color,v_class)
    node_cols=setNames(node_cols$color,node_cols$v_class)
    node_shapes=distinct(tmp_v,shape,v_class)
    node_shapes=setNames(node_shapes$shape,node_shapes$v_class)

    if(legend_number){
      eee=table(tmp_v$v_class)
      le_text= paste(vclass,eee[vclass],sep = ": ")
    }
    else le_text= vclass
  }

  p <- ggplot() +
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = e_type,
                     linewidth = width,linetype=e_class),data = tmp_e,alpha=0.7) +#draw edges
    scale_color_manual(name=edge_legend_title,values =edge_cols,
                       label=edge_text,guide=ifelse(edge_legend,"legend","none"))+#edge colors
    scale_linetype_manual(name=lty_legend_title,values =edge_ltys,
                          label=lty_text,guide=ifelse(lty_legend,"legend","none"))+#edge linetype
    scale_linewidth(name = width_legend_title,breaks =c(min(tmp_e$width),max(tmp_e$width)),range = c(0.5,1),
                    labels =edge_width_text,guide=ifelse(width_legend,"legend","none"))

  p1=p+
    geom_point(aes(X, Y,fill = v_class,size = size,shape=v_class), data = tmp_v) +#draw nodes
    #scale_shape_manual(values =setNames(pchls[node_shapes],vclass))+#node shape
    scale_shape_manual(values =node_shapes)+
    scale_fill_manual(name=class_legend_title,values =node_cols[vclass],
                      labels=le_text,guide=ifelse(col_legend,"legend","none"))+#node color
    scale_size(name=size_legend_title,breaks = c(min(tmp_v$size),max(tmp_v$size)),
               labels = node_size_text,guide=ifelse(size_legend,"legend","none"))+#node size

    ggnewscale::new_scale("size")+
    geom_text(aes(X, Y,size = size,label=label),col="black",data = tmp_v,show.legend = F)+
    scale_size(range = c(1,3),guide = "none")+
    guides(fill = guide_legend(override.aes = list(shape = node_shapes[vclass])),
           shape="none")

  p2=p1+labs(title = main)+
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    coord_fixed(ratio = 1)+
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA),
          legend.box.background = element_rect(colour = NA),
          legend.key = element_rect(fill = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

  if(!legend)return(p2+theme(legend.position = "none"))
  p2
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
  get_v(gephi)->tmp_v
  #extract coors
  if(!all(c("x","y","r","g","b","id")%in%colnames(tmp_v)))stop("The file is not exported by Gephi, please use igraph::read.graph()")
  coors<-tmp_v[,c("x","y")]
  coors<-data.frame(name = tmp_v$name,X=coors[,1],Y=coors[,2])
  coors%>%mutate(X=mmscale(X,-40,40),Y=mmscale(Y,-40,40))->coors
  #transform color
  pcutils::rgb2code(tmp_v[,c("r","g","b")])%>%dplyr::pull(code)->tmp_v$color
  E(gephi)$color=ifelse(E(gephi)$cor > 0, "#48A4F0", "#E85D5D")
  #scale size
  tmp_v$size=mmscale(tmp_v$size,1,5)
  E(gephi)$width=mmscale(E(gephi)$width,0.05,0.2)
  #delete
  tmp_v%>%dplyr::select(-c("label","x","y","r","g","b","id"))%>%as.list()->vertex.attributes(gephi)
  edge.attributes(gephi)["Edge Label"]=edge.attributes(gephi)["id"]=NULL
  gephi=c_net_update(gephi)
  return(list(go=gephi,coors=coors))
}

#' plot use networkD3
#'
#' @param go metanet
#' @param v_class which attributes use to be v_class
#' @param ... see \code{\link[networkD3]{forceNetwork}}
#'
#' @export
#'
#' @examples
#' plot(co_net2)
#' netD3plot(co_net2)
netD3plot<-function(go,v_class="v_class",...){
  flag="y"
  if(length(V(go))>200){
    print("Too big network, recommend using Gephi to layout,still use networkD3?")
    flag=readline("yes/no(y/n):")}
  if(tolower(flag)%in%c("yes","y")){
    lib_ps("networkD3")
    go=c_net_set(go,vertex_class = v_class)
    get_v(go)->tmp_v
    nodes=tmp_v[,c("name","v_class","size","color")]
    colnames(nodes)=c("name","group","size","color")
    nodes$size=pcutils::mmscale(nodes$size,2,40)

    colors=unique(nodes$color)

    get_e(go)->tmp_e
    links=tmp_e[,c("from","to","width","color")]
    links$width=pcutils::mmscale(links$width,0.5,1.5)
    #give ids
    links$IDsource <- match(links$from, nodes$name)-1
    links$IDtarget <- match(links$to, nodes$name)-1
    # Create force directed network plot
    networkD3::forceNetwork(Links = links, Nodes = nodes,
                            Source = 'IDsource', Target = 'IDtarget',linkColour = links$color,linkDistance=20,
                            linkWidth= JS("function(d) { return (d.width); }"),charge = -5,
                            NodeID = 'name', Group = 'group',Nodesize = 'size',
                            colourScale = JS(paste0("d3.scaleOrdinal([`",paste(colors,collapse = '`,`'),"`])")),legend = T,...)
  }
}

