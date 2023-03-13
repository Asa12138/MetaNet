#========6.modules=========

#' Generate a n-modules network
#'
#' @param n_modu number of modules
#' @param n_node_in_modu number of nodes in each modules
#' @param intra_modu_density intra_modu_density, recommend bigger than 20*inter_modu_density, default:0.3
#' @param inter_modu_density inter_modu_density, default:0.01
#'
#' @description this is just a random generation method, the module number of result is not exactly the n_modu, you can change the inter_modu_density and intra_modu_density to get the proper result.
#' @return n-modules metanet
#' @export
#'
#' @examples
#' g1=modu_net()
#' get_n(g1)
#' plot(g1)
#' plot(g1,coors=g_lay(g1,zoom2 = 20))
#' plot(g1,coors=g_lay_polyarc(g1,group = "module"))
#' plot(g1,coors=g_lay_polygon(g1,group = "module"))
modu_net<-function(n_modu=3,n_node_in_modu=30,
                   intra_modu_density=0.3,
                   inter_modu_density=0.01){
  n_node_in_modu=rep(n_node_in_modu,length=n_modu)
  mat=matrix(0,nrow = sum(n_node_in_modu),ncol = sum(n_node_in_modu))
  start=c(1,cumsum(n_node_in_modu[-length(n_node_in_modu)])+1)
  end=cumsum(n_node_in_modu)
  seqls=lapply(1:n_modu, \(i)(start[i]:end[i]))
  #generate intra_module_edges
  for (i in seq_len(n_modu)){
    seq=seqls[[i]]
    mat[seq,seq]=as.matrix(get.adjacency(erdos.renyi.game(n_node_in_modu[i], intra_modu_density)))
  }
  #generate inter_module_edges
  idls=combn(1:n_modu,2)%>%split(.,col(.))
  for (ids in idls) {
    seq1=seqls[[ids[1]]]
    seq2=seqls[[ids[2]]]
    mat[seq1,seq2]=sample(c(1,0),length(seq1)*length(seq2),replace = T,
           prob = c(inter_modu_density,1-inter_modu_density))
  }

  g <- graph.adjacency(mat, mode = "undirected")
  # plot(g)
  c_net_update(g)->g1
  g1=modu_dect(g1)
  g1=to_module_net(g1)
  g1
}

#' Detect the modules
#'
#' @param go a igraph object
#' @param method cluster_method:("cluster_walktrap","cluster_edge_betweenness","cluster_fast_greedy","cluster_spinglass")
#' @param n_modu transfer the modules less than n_modu to "others"
#' @return a igraph object
#' @export
#'
#' @examples
#' data("c_net")
#' modu_dect(co_net) -> co_net_modu
modu_dect <- function(go, method = "cluster_fast_greedy",n_modu=0,delete=F) {
  lib_ps("igraph")
  stopifnot(is_igraph(go))
  if("origin_module"%in%vertex_attr_names(go))message("'module' already exsited, start a new module detection!")
  ms <- c("cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass")
  method=match.arg(method,ms)
  if("weight"%in%graph_attr_names(go)) E(go)$weight=abs(igraph::E(go)$weight)

  switch (method,
          "cluster_walktrap" = {wc <- igraph::cluster_walktrap(go,weights =NULL)},
          "cluster_edge_betweenness" = {wc <- igraph::cluster_edge_betweenness(go,weights =NULL)},
          "cluster_fast_greedy" = {wc <- igraph::cluster_fast_greedy(go,weights =NULL)},
          "cluster_spinglass" = {wc <- igraph::cluster_spinglass(go,weights =NULL)})

  #add components
  V(go)$components=igraph::components(go)$membership%>%as.character()
  V(go)$module <- igraph::membership(wc)%>%as.character()
  V(go)$origin_module=V(go)$module
  go=filter_n_modu(go,n_modu =n_modu,delete = delete)
  go=c_net_update(go)

  graph_attr(go)$communities=wc
  graph_attr(go)$modularity<- modularity(wc)
  rand.g <- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
  rand_m=modularity(cluster_fast_greedy(rand.g))
  relative_modularity=(modularity(wc)-rand_m)/rand_m #relative modularity
  graph_attr(go)$relative_modularity=relative_modularity

  return(go)
}

as_module <- \(x){
  x <- list(membership = x)
  vids <- names(x$membership)
  modus <- tapply(vids, x$membership, simplify = FALSE, function(x) x)
  return(modus)
}

#' Filter some modules as others
#'
#' @param go_m module igraph
#' @param n_modu transfer the modules less than n_modu to "others"
#' @param keep_id keep modules ids, will not be "others"
#' @param delete logical, delete others modules? default:F, the others module will be "others".
#'
#' @export
#'
filter_n_modu<-function(go_m,n_modu=0,keep_id=NULL,delete=F){
  members=V(go_m)$origin_module
  table(members)%>%sort(decreasing = T)->s_members

  #filter modules whose nodes bigger than n_modu
  keep_id1=names(s_members[s_members>=n_modu])
  keep_id=union(as.character(keep_id),as.character(keep_id1))
  new_module=ifelse(members%in%keep_id,members,"others")
  V(go_m)$module=new_module
  if(delete)go_m=c_net_filter(go_m,module="others",exclude = T)
  return(go_m)
}

#' Combine or cut modules to modu_num
#'
#' @param go_m network
#' @param modu_num number of modules
#'
#' @export
#'
combine_n_modu<-function(go_m,modu_num=5){
  get_community(go_m)->comm
  igraph::cut_at(comm,modu_num)->new_modu
  V(go_m)$module=new_modu
  graph.attributes(go_m)$communities$membership=as.numeric(new_modu)
  go_m
}

#' Transformation a network to a module network
#'
#' @export
#'
to_module_net=function(go){
  if(!"module"%in%vertex_attr_names(go))stop("no modules, please modu_dect() first")
  go=anno_edge(go,get_v(go)[,c("name","module")])
  tmp_e=igraph::edge.attributes(go)
  E(go)$e_type=ifelse(tmp_e$module_from==tmp_e$module_to, "intra-module", "inter-module")
  go=c_net_set(go,vertex_class = "module")
  V(go)$color=ifelse(V(go)$module=="others","grey",V(go)$color)
  n_mod=unique(V(go)$module)
  igraph::graph.attributes(go)$n_type="module"
  igraph::graph.attributes(go)$n_modules=length(n_mod[n_mod!="others"])
  go
}

#' @export
get_community=function(go){
  if(is.null(igraph::graph_attr(go)$communities))stop("No community find, please do modu_net() first.")
  igraph::graph_attr(go)$communities
}

#' @export
get_module=function(go_m){
  if(!"module"%in%vertex_attr_names(go_m))stop("no modules, please modu_dect() first")
  setNames(V(go_m)$module,V(go_m)$name)
}

#' @export
get_module_eigen=function(go_m){
  graph_attr(co_net_pos_modu,"module_eigen")
}

#' @export
summ_module=function(go_m,var="v_class",module="module",...){
  tmp_v=get_v(go_m)
  if((length(module)>1)|(length(var)>1))stop("var or module should be one column!")
  a=tmp_v[,c(module,var)]
  colnames(a)[1]="module"

  i=var
  if(is.numeric(a[,i]))group_box(a[i],group = module,metadata = a,...)
  else {
    table(a[,i],a$module)%>%as.data.frame()%>%reshape2::acast(Var1~Var2,value.var = "Freq")%>%
      as.data.frame()->tab
    stackplot(tab,legend_title=var,...)+labs(x="Module")
  }
}

#' @export
plot_module_tree=function(go_m,class=c("module"),community=NULL,label.size=2){
  tmp_v=get_v(go_m)
  mdata=tmp_v[,c("name",class)]

  #modules tree
  if(is.null(community))get_community(go_m)%>%igraph::as_phylo()->mcl
  else community%>%igraph::as_phylo()->mcl

  lib_ps("ggtree")
  mcl=left_join(mcl,mdata,by=c("label"="name"))
  p=ggtree(mcl,size=0.3) + geom_tiplab(aes(color=module),show.legend=F,size=label.size)+
    scale_color_manual(values = pcutils::get_cols(length(unique(mdata$module)),"col3"))

  ggtree::gheatmap(p, mdata%>%tibble::column_to_rownames("name"))+
    scale_fill_manual(values = pcutils::get_cols(length(unique(mdata$module)),"col3"),name=NULL)
}

#' Calculate the eigenvalue of each module and correlation of nodes and eigenvalue (node_eigen_cor).
#'
#' @param go_m network
#' @param totu original abundance table
#'
#' @export
#'
module_eigen=function(go_m,totu){
  modules=get_module(go_m)
  totu=totu[,names(modules)]

  loop=\(i){
    if(i=="others")return(NULL)
    totu1=totu[,modules==i]
    #PCA
    # pc <- prcomp(totu1)
    # pc$x[, 1]
    # cor(pc$x[, 1],rowMeans(totu1))

    #SVD
    totu_scale=scale(totu1)
    svd=svd(totu_scale)
    if(cor(svd$u[,1],rowMeans(totu_scale))<0){
      return(-svd$u[,1])
    }
    else return(svd$u[,1])
  }

  res <-lapply(levels(factor(modules)), loop)

  #simplify method
  eigen_res=do.call(cbind,res)%>%data.frame(.,check.names = F)
  #colnames(eigen_res)=paste0("ME_",colnames(eigen_res))
  rownames(eigen_res)=rownames(totu)

  lapply(colnames(eigen_res), \(i)cor(totu[,modules==i],eigen_res[,i]))%>%
    do.call(rbind,.)%>%as.data.frame()->node_eigen_cor
  colnames(node_eigen_cor)="node_eigen_cor"
  go_m=anno_vertex(go_m,node_eigen_cor)
  graph.attributes(go_m)$module_eigen=eigen_res
  go_m
}

#' Plot the expression of each modules
#'
#' @param go_m network
#' @param totu original abundance table used for module_eigen().
#' @param group group variable for totu
#' @param cor the threshold for node_eigen_cor
#'
#' @export
#'
#' @examples
#' module_eigen(co_net_modu,totu)->co_net_modu
#' module_expression(co_net_modu,totu)
module_expression=function(go_m,totu,group=NULL,cor=0.6){
  if(is.null(graph_attr(go_m,"module_eigen")))stop("Please do module_eigen() first")

  # graph_attr(go_m,"module_eigen")%>%rownames_to_column()%>%melt()%>%
  #   ggplot(.,aes(x=rowname,y=value,col=variable,group=variable))+geom_point()+geom_line()

  get_v(go_m)->tmp_v

  if(!is.null(group))totu=hebing(totu,group,1)

  totu_scale=scale(totu)
  pdat=cbind(tmp_v[,c("name","module","node_eigen_cor")],t(totu_scale[,tmp_v$name]))
  pdat=filter(pdat,node_eigen_cor>cor)
  pdat_m=melt(pdat,id.vars = c("name","module","node_eigen_cor"))
  pdat_m$module=factor(pdat_m$module)

  ggplot(data=pdat_m, aes(x=variable, y=value, group=name, color=module,
                          size=node_eigen_cor,alpha=node_eigen_cor)) +
    geom_line(size=0.8)+ theme_pubr(base_size = 12,legend = 'right') +
    facet_wrap(~module)+
    scale_x_discrete(expand=c(0,0))+theme(plot.margin=unit(c(1,2,1,1),'lines'))+
    scale_color_manual(values = get_cols(nlevels(pdat_m$module)))
}

#' Zi-Pi calculate
#'
#' @param go igraph object after modu_dect()
#' @param mode use 7-group (mode=1) or 4-group (mode=2)
#' @param use_origin use origin_module, default:T, if F, use module
#'
#' @return igraph
#' @export
#'
#' @references 1. Guimerà, R. & Amaral, L. Functional cartography of complex metabolic networks. (2005) doi:10.1038/nature03288.
#' @examples
#' data("c_net")
#' modu_dect(co_net) -> co_net_modu
#' zp_analyse(co_net_modu)->co_net_modu
#' zp_plot(co_net_modu)
#' zp_plot(co_net_modu,mode=3)
zp_analyse<-function(go,mode=2,use_origin=T){
  go->go1
  v_index=get_v(go)
  if(!"module"%in%names(v_index))stop("no modules, please modu_dect() first")
  if("roles"%in%names(v_index))message("areadly has roles, overwrite!")

  #use origin_module to do zp_analyse
  if(use_origin){go1=anno_vertex(go1,data.frame(row.names = V(go1)$name,module=V(go1)$origin_module%>%as.numeric()))}%>%suppressMessages()
  #or consider others as a big module?
  else {
    if("others"%in%v_index$module)message("Consider others as one module!")
    {go1=anno_vertex(go1,data.frame(row.names = V(go1)$name,
                                       module=tidai(v_index$module,seq_along(unique(v_index$module)))))}%>%suppressMessages()}

  within <- within_module_deg_z_score(go1)
  v_index$Ki=within$Ki;v_index$Zi=within$Zi
  pc <- part_coeff(go1)
  v_index$Pi=pc$Pi

  if(mode==1){
    lab <- c("Ultra-peripherals","Peripherals","Non-hub connectors","Non-hub kinless nodes","Provincial hubs", "Connector hubs","Kinless hubs")
    backs=data.frame(
      x1= c(0,0.05, 0.62,0.8, 0,0.3,0.75),
      x2= c(0.05,0.62,0.8, 1, 0.3,0.75, 1),
      y1= c(-Inf,-Inf,-Inf,-Inf, 2.5, 2.5,2.5),
      y2= c(2.5, 2.5, 2.5,2.5, Inf, Inf, Inf),
      lab=factor(lab,levels = lab)
    )
  }
  else if(mode==2){
    lab <- c("Peripherals", "Network hubs", "Module hubs","Connectors")
    backs=data.frame(
      x1= c(0, 0.62, 0, 0.62),
      x2= c(0.62, 1, 0.62, 1),
      y1= c(-Inf, 2.5, 2.5, -Inf),
      y2= c(2.5, Inf, Inf, 2.5),
      lab=factor(lab,levels = lab)
    )
  }
  deter_role=\(x,y,backs=backs){
    for (i in 1:nrow(backs)){
      if(dplyr::between(x,backs$x1[i],backs$x2[i])&dplyr::between(y,backs$y1[i],backs$y2[i])){
      #if((backs$x1[i]<=x)&(backs$x2[i]>=x)&(backs$y1[i]<=y)&(backs$y2[i]>=y)){
        role=backs$lab[i]
        break
      }
    }
    return(role)
  }
  v_index$roles=apply(v_index, 1, \(x)deter_role(x["Pi"],x["Zi"],backs))
  vertex.attributes(go)<-as.list(v_index)
  return(go)

}

# calculate Zi
within_module_deg_z_score <- function(g, A = NULL, weighted = FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse = FALSE, names = TRUE, attr = "weight")
    } else {
      A <- as_adj(g, sparse = FALSE, names = TRUE)
    }
  }
  memb <- vertex_attr(g, "module")
  N <- max(memb)
  nS <- tabulate(memb)
  z <- Ki <- rep.int(0, dim(A)[1L])
  Ksi <- sigKsi <- rep.int(0, N)
  names(z) <- names(Ki) <- rownames(A)
  for (S in seq_len(N)) {
    x <- rowSums(as.matrix(A[memb == S, memb == S]))
    Ki[memb == S] <- x
    Ksi[S] <- sum(x) / nS[S]
    sigKsi[S] <- sqrt(sum((x - Ksi[S])^2) / (nS[S] - 1))
  }
  z <- (Ki - Ksi[memb]) / sigKsi[memb]
  z[is.infinite(z)] <- 0
  z[is.nan(z)] <- 0
  Zi <- z
  df <- data.frame(Ki, Zi, row.names = names(Ki))
  return(df)
}

# calculate Pi
part_coeff <- function(g, A = NULL, weighted = FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse = FALSE, attr = "weight")
    } else {
      A <- as_adj(g, sparse = FALSE)
    }
  }
  memb <- vertex_attr(g, "module")
  Ki <- colSums(A)
  N <- max(memb)
  Kis <- t(rowsum(A, memb))
  Pi <- 1 - ((1 / Ki^2) * rowSums(Kis^2))
  names(Pi) <- rownames(A)
  Pi <- data.frame(Pi)
  return(Pi)
}

#' Zi-Pi plot of vertexes
#'
#' @param go igraph object after zp_analyse()
#' @param label show label or not
#' @param mode 1,2,3
#'
#' @return a ggplot object
#' @export
#' @import ggpubr ggrepel
#' @seealso \link{zp_analyse}
zp_plot <- function(go,label=T,mode=1) {
  lib_ps("ggrepel","ggpubr","igraph")
  as.data.frame(vertex.attributes(go))->taxa.roles
  if(!"roles"%in%names(taxa.roles))stop("no roles, please zp_analyse() first")
  if(mode==3){
    reshape2::melt(taxa.roles,measure.vars = c("Zi","Pi"))->taxa.roles1
    p=ggplot(taxa.roles1,aes(x=v_class,y=value,col=v_class,size=size,shape=roles))+geom_point()+
      facet_grid(variable~.,scales = "free_y")+theme_bw()+labs(x=NULL,y=NULL)+guides(col="none")+
      theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
      scale_color_manual(values =get_cols(length(unique(taxa.roles1$v_class)),"col3"))
    return(p)
  }

  mode=ifelse(nlevels(V(go)$roles)==7,1,2)
  if(mode==1){
    lab <- c("Ultra-peripherals","Peripherals","Non-hub connectors","Non-hub kinless nodes","Provincial hubs", "Connector hubs","Kinless hubs")
    CPCOLS<-c("#FCF6EFFC","#EEBCF5", "#EDEDA4", "#FAA371","#FC5D6096" ,"#9BC799B9" ,"#94CCF2AC")
    names(CPCOLS)<-lab
    backs=data.frame(
      x1= c(0,0.05, 0.62,0.8, 0,0.3,0.75),
      x2= c(0.05,0.62,0.8, 1, 0.3,0.75, 1),
      y1= c(-Inf,-Inf,-Inf,-Inf, 2.5, 2.5,2.5),
      y2= c(2.5, 2.5, 2.5,2.5, Inf, Inf, Inf),
      lab=factor(lab,levels = lab)
    )
  }
  else if(mode==2){
    lab <- c("Peripherals", "Network hubs", "Module hubs","Connectors")
    CPCOLS<-c("#FCF6EFFC" ,"#FC5D6096" ,"#9BC799B9" ,"#94CCF2AC")
    names(CPCOLS)<-lab
    backs=data.frame(
      x1= c(0, 0.62, 0, 0.62),
      x2= c(0.62, 1, 0.62, 1),
      y1= c(-Inf, 2.5, 2.5, -Inf),
      y2= c(2.5, Inf, Inf, 2.5),
      lab=factor(lab,levels = lab)
    )
  }
  p <- ggplot() +
    geom_rect(data = backs, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = lab)) +
    guides(fill = guide_legend(title = "Topological roles")) +
    scale_fill_manual(values = CPCOLS) +
    geom_point(data = taxa.roles, aes(x = Pi, y = Zi, color = factor(v_class))) +
    scale_color_manual(values =get_cols(length(unique(taxa.roles$v_class)),"col3"))+
    ggpubr::theme_pubr(legend = "right") +
    guides(colour = "none") +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Participation Coefficient") +
    ylab("Within-module connectivity z-score")

  if(label){
    p=p+ggrepel::geom_text_repel(
      data =taxa.roles[!taxa.roles$roles%in%c("Peripherals","Ultra-peripherals"),],
      aes(x = Pi, y = Zi, label = name),size=3
    )
  }
  return(p)
}
