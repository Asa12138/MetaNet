#=========2.build======

#' Multi-omics network build
#'
#' @param ... some omics abundance tables
#' @param mode "full"
#' @param method spearman, pearson
#' @param filename the prefix of saved files
#' @param p.adjust.method see \code{\link[stats]{p.adjust}}
#'
#' @return metanet
#' @export
#'
#' @examples
#' data("multi_test")
#' multi_net_build(micro,metab,transc)->multi1
#' multi1=c_net_set(multi1,micro_g,metab_g,transc_g,vertex_class = c("Phylum","kingdom","type"))
#' multi1=c_net_set(multi1,data.frame("Abundance1"=colSums(micro)),data.frame("Abundance2"=colSums(metab)),data.frame("Abundance3"=colSums(transc)),vertex_size =paste0("Abundance",1:3))
#' c_net_plot(multi1)
#' multi1=c_net_set(multi1,vertex_class = c("kingdom","type"))
#' c_net_plot(multi1)
#'
multi_net_build<-function(...,mode="full",method = "spearman",
                          filename = "multi-net",p.adjust.method=NULL,
                          r_thres = 0.6, p_thres = 0.05,use_p_adj=T, del_single = T){
  tables=list(...)
  if(mode=="full"){
    all_totu=do.call(cbind,tables)
    all_corr=c_net_cal(all_totu,method = method,filename = filename,p.adjust.method = p.adjust.method)

    c_net_build(all_corr, r_thres = r_thres, p_thres = p_thres,use_p_adj=use_p_adj, del_single =del_single)->multi_net
    graph.attributes(multi_net)$n_type="multi_full"
    get_v(multi_net)->tmp_v

    position=rep(paste0("omic",seq_len(length(tables))),sapply(tables, ncol))
    names(position)=lapply(tables, colnames)%>%do.call(c,.)

    tmp_v$v_class=tmp_v$v_group=sapply(tmp_v$name, \(x)position[x])
    tmp_v%>%as.list()->vertex_attr(multi_net)
    # set edges from_to
    suppressMessages(anno_edge(multi_net,tmp_v[,c("name","v_group")])->multi_net)
    # set edge intra-inter
    E(multi_net)$e_class <- get_e(multi_net)[,c("v_group_from","v_group_to")] %>%
      apply(., 1, \(x)ifelse(all_same(x), "intra", "inter"))
    c_net_update(multi_net)->multi_net

  }
  multi_net
}

#' Construct a network from correlation table
#'
#' @param cor_res result of c_net_cal()
#' @param r_thres r_threshold (default:>0.6)
#' @param p_thres p_threshold (default:<0.05)
#' @param use_p_adj use the p.adjust instead of p-value (default: T)
#' @param del_single should delete single vertexes?
#'
#' @return an igraph object
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab) -> totu
#' metadata[,3:10] -> env
#' c_net_cal(totu) -> corr
#' c_net_build(corr,r_thres=0.65) -> co_net
#' c_net_cal(totu,env) -> corr2
#' c_net_build(corr2) -> co_net2
c_net_build <- function(cor_res, r_thres = 0.6, p_thres = 0.05,use_p_adj=T, del_single = T) {
  suppressMessages(lib_ps("igraph"))
  # set thresholds to construct
  if("r"%in%names(cor_res))  occor.r <- cor_res$r
  else stop("No r in the input object")
  if(("p.adjust"%in%names(cor_res))&use_p_adj)occor.p <- cor_res$p.adjust
  else {
    if("p.value"%in%names(cor_res)){
      occor.p <- cor_res$p.value;
      message("Have not do p-value adjust! use the p.value to build network")}
    else {
      occor.p =occor.r;occor.p[occor.p!=0]=0
      message("No p.value given, just use r threshold to build network!")}
  }

  occor.r[occor.p > p_thres | abs(occor.r) < r_thres] <- 0
  corr <- occor.r

  # make igraph
  # flag determined by the correlation table from one table or two
  flag <- \(corr){
    if (!nrow(corr) == ncol(corr)) {
      return(FALSE)
    }
    if (!all(rownames(corr) == colnames(corr))) {
      return(FALSE)
    }
    return(TRUE)
  }

  if (flag(corr)) {
    go <- igraph::graph_from_adjacency_matrix(as.matrix(corr), mode = "undirected", weighted = T, diag = F)
    igraph::graph.attributes(go)$n_type="single"
  }
  else {
    inset=intersect(rownames(corr),colnames(corr))
    if(length(inset)>0){
      warning("some nodes has same name, please check: ",paste0(inset,collapse = ","))
    }
    go <- igraph::graph_from_incidence_matrix(as.matrix(corr), directed = F, weighted = T)
    igraph::graph.attributes(go)$n_type="bipartite"
  }

  # delete single vertexes?
  if (del_single) go <- igraph::delete.vertices(go, V(go)[igraph::degree(go) == 0])

  # set vertex attributes
  # set vertices shape
  V(go)$v_group <- ifelse(V(go)$name %in% rownames(corr), "v_group1", "v_group2")
  V(go)$v_class <- V(go)$v_group
  V(go)$size <- ceiling(60/sqrt(length(V(go))))+1
  V(go)$label=V(go)$name

  # abs edges weight
  go.weight <- E(go)$weight
  E(go)$cor <- go.weight
  E(go)$weight <- abs(go.weight)
  if("p.adjust"%in%names(cor_res))E(go)$p.adjust <-get_e(go)%>%select(from,to)%>%apply(., 1, \(x)cor_res$p.adjust[x[1],x[2]])
  if("p.value"%in%names(cor_res))E(go)$p.value <-get_e(go)%>%select(from,to)%>%apply(., 1, \(x)cor_res$p.value[x[1],x[2]])

  # set edges type
  E(go)$e_type<-ifelse(go.weight > 0, "positive", "negative")
  # set edges width
  E(go)$width <- E(go)$weight

  # set edges from_to
  anno_edge(go,get_v(go)[,c("name","v_group")])->go
  # set edge intra-inter
  E(go)$e_class <- get_e(go)[,c("v_group_from","v_group_to")] %>%
    apply(., 1, \(x)ifelse(all_same(x), "intra", "inter"))
  c_net_update(go)->go

  return(go)
}

#' @export
c_net_update<-function(go){
  #name
  if(!"name"%in%vertex_attr_names(go))V(go)$name=paste0("n",seq_len(length(go)))
  get_v(go)->tmp_v
  #v_size
  if(!"size"%in%colnames(tmp_v)){tmp_v$size=ceiling(60/sqrt(length(V(go))))+1}

  #v_shape
  if("v_group"%in%colnames(tmp_v)){
    tmp_v$shape <- tidai(tmp_v$v_group,c("circle","square"))
  }
  else {tmp_v$v_group="v_group1";tmp_v$shape="circle"}

  #v_color

  if(!"v_class"%in%colnames(tmp_v)){tmp_v$v_class="v_class1"}
  tmp_col=paste0(tmp_v$v_group,"-",tmp_v$v_class)
  tmp_v$color =tidai(tmp_col,pcutils::get_cols(nlevels(factor(tmp_col)),"col3"))

  as.list(tmp_v)->vertex_attr(go)

  #e_color
  if(!"e_type"%in%edge_attr_names(go)){E(go)$e_type="e_type1"}
  edge.color <- droplevels(as.factor(E(go)$e_type))

  if(all(levels(edge.color)%in%c("negative","positive")))ncols=c(negative="#E85D5D",positive="#48A4F0")
  else if(all(levels(edge.color)%in%c("inter","intra")))ncols=c(inter="#FA789A",intra="#A6CEE3")
  else ncols=pcutils::get_cols(nlevels(edge.color),"col2")
  E(go)$color=tidai(E(go)$e_type,ncols)

  #e_lty
  if(!"e_class"%in%edge_attr_names(go)){E(go)$e_class="e_class1"}
  E(go)$lty=tidai(E(go)$e_class,1:4)

  #e_width
  if(!"width"%in%edge_attr_names(go)){E(go)$width=1}

  class(go)=c("metanet","igraph")
  return(go)
}


#' Construct a network from Edge_list dataframe
#'
#' @param edgelist first is source, second is target, others are annotation
#' @param direct logical
#' @param e_type set e_type
#' @param e_class set e_class
#'
#' @export
c_net_from_edgelist<-function(edgelist,direct=F,e_type=NULL,e_class=NULL){
  lib_ps("igraph")
  go=igraph::graph_from_data_frame(edgelist,directed = direct)
  if(!is.null(e_type))E(go)$e_type=edgelist[,e_type]
  if(!is.null(e_class))E(go)$e_class=edgelist[,e_class]
  go=c_net_update(go)
  go
}


#' Set basic attributes from totu table
#'
#' @param go a igraph object
#' @param ... some dataframes to annotate go
#' @param vertex_class choose which column to be vertex_class (map to vertex_color)
#' @param vertex_size choose which column to be vertex_size (map to vertex_size)
#'
#' @return a igraph object
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab) -> totu
#' metadata[,3:10] -> env
#'
#' data("c_net")
#' co_net<-c_net_set(co_net,taxonomy,data.frame("Abundance"=colSums(totu)),vertex_class="Phylum",vertex_size="Abundance")
#' co_net_rmt<-c_net_set(co_net_rmt,taxonomy,data.frame("Abundance"=colSums(totu)),vertex_class="Phylum",vertex_size="Abundance")
#' co_net2<-c_net_set(co_net2,taxonomy,data.frame(name=colnames(env),env=colnames(env)),vertex_class=c("Phylum","env"))
#' co_net2<-c_net_set(co_net2,data.frame("Abundance"=colSums(totu)),vertex_size="Abundance")
c_net_set <- function(go,...,vertex_group="v_group",vertex_class="v_class",vertex_size="size",
                      edge_type="e_type",edge_class="e_class",edge_width="width") {
  lib_ps("igraph")
#annotation vertex
  anno_dfs<-list(...)
  if(length(anno_dfs)>0){
    anno_dfs2=list()
    for (i in seq_len(length(anno_dfs))) {
      x=anno_dfs[[i]]
      if ("name" %in% colnames(x)) {rownames(x)<-x$name;x=select(x,-name)}
      anno_dfs2[[i]]=x
    }

    if(any(duplicated(sapply(anno_dfs2, names))))stop("duplicated name in your annotation data.frame!")

    Reduce(\(x,y)merge(x,y,by = "row.names", all = T)%>%tibble::column_to_rownames("Row.names"),
           anno_dfs2)->all_anno

    anno_vertex(go,all_anno)->go
  }
  get_v(go)->v_index
  get_e(go)->e_index
#set something
  vec_equal=\(x,y){
    if (length(x)!=length(y))return(F)
    all(x==y)
  }

  if(!vec_equal(vertex_group,"v_group"))select(v_index,v_group,!!vertex_group)%>%condance->V(go)$v_group
  if(!vec_equal(vertex_class,"v_class"))select(v_index,v_class,!!vertex_class)%>%condance->V(go)$v_class
  if(!vec_equal(vertex_size,"size"))select(v_index,size,!!vertex_size)%>%condance->V(go)$size
  if(!vec_equal(edge_type,"e_type"))select(e_index,e_type,!!edge_type)%>%condance->E(go)$e_type
  if(!vec_equal(edge_class,"e_class"))select(e_index,e_class,!!edge_class)%>%condance->E(go)$e_class
  if(!vec_equal(edge_width,"width"))select(e_index,width,!!edge_width)%>%condance->E(go)$width

  c_net_update(go)->go
  return(go)
}


tidai=\(x,y){
  tmp=y
  if(is.null(names(tmp))){
    tmp=rep(tmp,len=length(unique(x)))
    names(tmp)=unique(x)
  }
  return(unname(tmp[x]))
}
all_same <- \(x){
  return(all(x == x[1]))
}
condance<-\(aa){
  apply(aa, 1, \(x){
  tmp=x[!is.na(x)]
  if(length(tmp)==0)return(NA)
  return(tmp[length(tmp)])
})}
change_fac_lev<-\(x,level=NULL){
  ordervec=factor(x)
  if(!is.null(level)){
    level=intersect(level,levels(ordervec))
    shunxu=c(level,setdiff(levels(ordervec),level))
    ordervec=factor(ordervec,levels = shunxu)
  }
  ordervec
}

#==========2.1manipulate========
#' @export
get_v=function(go){
  as.data.frame(igraph::vertex.attributes(go))
}
#' @export
get_e=function(go){
  tmp_e=cbind_new(igraph::as_data_frame(go),as.data.frame(igraph::edge_attr(go)))
  tmp_e=cbind_new(tmp_e,data.frame(id=1:nrow(tmp_e)))
  dplyr::select(tmp_e,id,everything())
}
#' @export
get_n=function(go){
  gls=igraph::graph.attributes(go)
  gls=lapply(gls,\(x)ifelse(is.data.frame(x),"data.frame",
                        ifelse(length(x)>1,paste(x,collapse = ";"),x)))
  as.data.frame(gls)
}

#' Filter a network according to some attributes
#'
#' @param go igraph
#' @param ... some attributes of vertex and edge
#'
#' @return igraph
#' @export
#'
#' @examples
#'
#' c_net_filter(go,v_group=c("omic1","omic2"),e_class="intra")
c_net_filter<-function(go,...){
  f_ls=list(...)
  v_ls=list();e_ls=list()
  for (i in names(f_ls)) {
    if((i %in% igraph::vertex_attr_names(go))){
      if((i %in% igraph::edge_attr_names(go))){warning("both vertex and edge has attributes ",i," use `filter_e()` to filter edge")}
      v_ls[[i]]=f_ls[[i]]
    }
    else if(i %in% igraph::edge_attr_names(go)){
      e_ls[[i]]=f_ls[[i]]
    }
  }
  go1=filter_v(go,v_ls)
  go1=filter_e(go1,e_ls)
  go1
}

#' @export
filter_v<-function(go,...){
  f_ls=list(...)
  if(is.list(f_ls[[1]]))f_ls=f_ls[[1]]
  get_v(go)->tmp_v
  for (i in names(f_ls)) {
    if(!i%in%colnames(tmp_v)){warning(i," do not in vertex_attributes, skip");next}
    tmp_v=tmp_v[tmp_v[[i]]%in%f_ls[[i]],]
  }

  tmp_v$name->vid
  igraph::subgraph(go,vid)->go1
  go1
}

#' @export
filter_e<-function(go,...){
  f_ls=list(...)
  if(is.list(f_ls[[1]]))f_ls=f_ls[[1]]
  get_e(go)->tmp_e
  for (i in names(f_ls)) {
    if(!i%in%colnames(tmp_e)){warning(i," do not in edge_attributes, skip");next}
    tmp_e=tmp_e[tmp_e[[i]]%in%f_ls[[i]],]
  }

  tmp_e$id->eid
  igraph::subgraph.edges(go,eid)->go1
  go1
}


#' Use dataframe to annotate vertexes of a igraph
#'
#' @param go a igraph object
#' @param anno_tab a dataframe using to annotate (with rowname or a name column)
#' @return a annotated igraph object
#' @export
#'
#' @examples
#' data("c_net")
#' data("otutab")
#' anno_vertex(co_net, taxonomy)->a
anno_vertex <- function(go, anno_tab) {
  lib_ps("igraph")
  if(is.null(anno_tab))return(go)
  get_v(go) -> v_atr
  if (!"name" %in% colnames(anno_tab)) rownames(anno_tab) -> anno_tab$name
  v_atr <- dplyr::left_join(v_atr, anno_tab, by = "name", suffix = c(".x", ""))
  grep(".x",colnames(v_atr),value = T)%>%gsub(".x","",.)->du
  if(length(du)>0)message(length(du),(" attributes will be overwrited:\n"),paste0(du,collapse = ","),"\n")
  v_atr %>% dplyr::select(!ends_with(".x")) -> v_atr

  as.list(v_atr) -> igraph::vertex.attributes(go)
  return(go)
}

#' Use dataframe to annotate edges of a igraph
#'
#' @param go a igraph object
#' @param anno_tab a dataframe using to annotate (with rowname or a name column)
#' @return a annotated igraph object
#' @export
#'
#' @examples
#' data("c_net")
#' data("otutab")
#' anno_edge(co_net, taxonomy)->a
anno_edge <- function(go, anno_tab) {
  lib_ps("igraph")
  if(is.null(anno_tab))return(go)
  get_e(go) -> e_atr
  if (!"name" %in% colnames(anno_tab)) rownames(anno_tab) -> anno_tab$name
  anno_tab%>%select(name,everything())->anno_tab
  #from
  tmp=anno_tab
  colnames(tmp)=paste0(colnames(anno_tab),"_from")
  e_atr <- dplyr::left_join(e_atr, tmp, by = c("from"="name_from"), suffix = c(".x", ""))
  grep(".x",colnames(e_atr),value = T)%>%gsub(".x","",.)->du
  if(length(du)>0)message(length(du),(" attributes will be overwrited:\n"),paste0(du,collapse = ","),"\n")
  e_atr %>% dplyr::select(!ends_with(".x")) -> e_atr
  #to
  tmp=anno_tab
  colnames(tmp)=paste0(colnames(anno_tab),"_to")
  e_atr <- dplyr::left_join(e_atr, tmp, by = c("to"="name_to"), suffix = c(".x", ""))
  grep(".x",colnames(e_atr),value = T)%>%gsub(".x","",.)->du
  if(length(du)>0)message(length(du),(" attributes will be overwrited:\n"),paste0(du,collapse = ","),"\n")
  e_atr %>% dplyr::select(!ends_with(".x")) -> e_atr

  as.list(e_atr) -> igraph::edge.attributes(go)
  return(go)
}

#' Save network file
#'
#' @param go network
#' @param filename filename
#' @param format "data.frame"
#'
#' @export
#'
c_net_save<-function(go,filename="net",format="data.frame"){
  if(format=="data.frame"){
    get_v(go)%>%write.csv(.,paste0(filename,"_nodes.csv"),row.names = F)
    get_e(go)%>%select(-1)%>%write.csv(.,paste0(filename,"_edges.csv"),row.names = F)
  }
  if(format=="graphml"){
    igraph::write.graph(go,paste0(file,".graphml"),format = "graphml")
  }
  print(paste0(file," saved sucessfully!"))
}


#' Summaries two columns information
#' @param df dataframe
#' @param from first column name or index
#' @param to second column name or index
#' @param count (optional) weight column, if no, each equal to 1
#' @param direct consider direct? defualt: F
#'
#' @export
#' @examples
#' test=data.frame(a=sample(letters[1:4],10,replace = T),b=sample(letters[1:4],10,replace = T))
#' summ_2col(test,direct=T)
#' summ_2col(test,direct=F)
#' summ_2col(test,direct=T)%>%my_sankey()
summ_2col<-function(df,from=1,to=2,count=3,direct=F){
  lib_ps("dplyr")
  if(ncol(df)<2)stop("need at least two columns")
  if(ncol(df)==2)tmp=cbind(df,count=1)
  else  tmp=select(df,!!from,!!to,!!count)
  cols=colnames(tmp)
  colnames(tmp)=c("from","to","count")

  if(direct){
    tmp=(group_by(tmp,from,to)%>%summarise(count=sum(count)))
    colnames(tmp)=cols
    return(as.data.frame(tmp))
  }

  com=\(group1,group2,levels){
    factor(c(group1,group2),levels = levels)%>%sort
  }
  group=factor(c(tmp[,1],tmp[,2]))
  tmp1<-apply(tmp,1,function(x)com(x[1],x[2],levels(group)))%>%t%>%as.data.frame()
  tmp1=cbind(tmp1,tmp$count)
  colnames(tmp1)=c("from","to","count")
  tmp1=group_by(tmp1,from,to)%>%summarise(count=sum(count))
  colnames(tmp1)=cols
  return(as.data.frame(tmp1))
}


#=========2.2RMT optimize=====
#' Get RMT threshold for a correlation matrix
#'
#' @param occor.r a correlation matrix
#' @param min_threshold min_threshold
#' @param max_threshold max_threshold
#' @param step step
#' @param gif render .gif file?
#'
#' @return a r-threshold
#' @export
#' @references 1.  J. Zhou, Y. Deng, F. Luo, Z. He, Q. Tu, X. Zhi, Functional Molecular Ecological Networks (2010), doi:10.1128/mBio.00169-10.
#' @references 2. \link{https://matstat.org/content_en/RMT/RMThreshold_Intro.pdf}
#' @examples
#' t(otutab) -> totu
#' c_net_cal(totu) -> corr
#' rmt(corr$r)
#' #0.69
#' c_net_build(corr,r_thres=0.69) -> co_net_rmt
#' RMT_threshold(corr$r,gif=T)->rmt_res
RMT_threshold = function(occor.r, min_threshold = 0.5, max_threshold = 0.9, step = 0.02,
                         gif=F,quite=F){
  nwd=getwd()
  if(!dir.exists("./RMT_temp"))dir.create("./RMT_temp")

  if(max_threshold>=max(abs(occor.r)))max_threshold=(max(abs(occor.r))-step)
  if(min_threshold>=max_threshold)min_threshold=max_threshold-10*step

  thres_seq=seq(min_threshold, max_threshold, step)

  res <- data.frame()
  for(i in seq_len(length(thres_seq))){
    threshold=thres_seq[i]
    if(!quite)dabiao(paste0("Calculating",i,":  threshold =", signif(threshold,3)))
    corr_r1 <- occor.r
    corr_r1[abs(corr_r1) < threshold] <- 0
    #calculate eigenvalues
    rand.mat=corr_r1
    eigenvalues = eigen(rand.mat, only.values = T)$values
    eigenvalues = eigenvalues[order(eigenvalues)]/max(abs(eigenvalues))

    eigenvalues = remove.outliers(unique(eigenvalues))
    #get the NNDS
    {#uf <- rm.unfold.gauss(eigenvalues,pop.up = T)
      dens = density(eigenvalues,kernel = "gaussian")
      midpoints <- \(x)(x[-length(x)] + 0.5 * diff(x))
      scale.function = approx(dens$x, dens$y, xout = midpoints(eigenvalues))
      ev.spacing = diff(eigenvalues)
      ev.spacing = ev.spacing * scale.function$y
      ev.spacing = ev.spacing/mean(ev.spacing)
    }
    if(F){ #uf <- rm.unfold.spline(eigenvalues,pop.up = T)
      cumulative.distrib = ecdf(eigenvalues)
      support = seq(min(eigenvalues), max(eigenvalues), len = nr.fit.points)
      D = splinefun(support, cumulative.distrib(support), method = "hyman")
      uf.ev = D(eigenvalues)
      ev.spacing = diff(uf.ev)
      ev.spacing = ev.spacing/mean(ev.spacing)
    }
    if(F){
      ssp <- smooth.spline(eigenvalues, control.spar = list(low = 0, high = 3))
      ev.spacing=NNSD(ssp$y)/mean(NNSD(ssp$y))
    }
    ev.spacing = ev.spacing[ev.spacing <= 3]
    #test whether fit possion?
    p_ks_test = ks.test(unique(ev.spacing), "pexp", 1)$p.value
    #get sse
    #sse = rm.sse(ev.spacing)
    sse = get_sse(ev.spacing)
    log_sse=log(sse)

    if(F){nnsdw <- density(ev.spacing)
    poisson_d <- exp(seq(0, 3, len=1000))
    nnsdpois <- density(NNSD(poisson_d))
    chival <- sum((nnsdw$y - nnsdpois$y)^2/nnsdpois$y/1e3)}

    #maximum likehood
    evs = ev.spacing[ev.spacing != 0]
    N = length(evs)
    log_LE = -sum(evs)/N
    log_LW = log(pi/2) + sum(log(evs))/N - 0.25 * pi * sum(evs^2)/N

    #save pngs
    {
      histo <- hist(ev.spacing, breaks = seq(min(ev.spacing),max(ev.spacing), len = 51), plot = F)
      png(paste0("RMT_temp/rmt_nnsd",i,".png"),height = 600,width = 700,res = 130)
      nnsd_plot(histo = histo, title = "Eigenvalue spacing distribution (NNSD)", threshold = threshold,
                dis_GOE = log_LW, dis_possion = log_LE,p_ks_test = p_ks_test)
      dev.off()
    }
    res=rbind(res,data.frame(threshold,p_ks_test,log_sse,log_LW,log_LE))
  }
  #transfer to gif
  if(gif){
    lib_ps("gifski")
    gifski::gifski(paste0("RMT_temp/rmt_nnsd",seq_len(length(thres_seq)),".png"),gif_file = "RMT_temp/rmt_nnsd.gif")
  }
  message(paste("The Intermediate files saved in ./RMT_temp/ ."))
  class(res)=c("rmt_res",class(res))
  return(res)
}

#' Plot a rmt_res
#'
#' @return ggplot
#' @exportS3Method
plot.rmt_res=function(res){
  linedf=data.frame(variable=c("p_ks_test","log_sse","log_LW","log_LE"),
                    xi=c(res[which(res$p_ks_test==max(res$p_ks_test)),"threshold"],
                         res[which(res$log_sse==min(res$log_sse)),"threshold"],
                         res[which(res$log_LW==min(res$log_LW)),"threshold"],
                         res[which(res$log_LE==max(res$log_LE)),"threshold"]),
                    x=max(res$threshold)-min(res$threshold),
                    y=apply(res[,-1], 2, max))

  reshape2::melt(res,"threshold")->md

  #filter(threshold<0.77)%>%
  p=ggplot(md,aes(threshold,value))+geom_point(aes(col=variable))+
    geom_line(aes(col=variable))+scale_color_manual(values = get_cols(4,"col1"))+
    facet_wrap(.~variable,scales = "free_y")+theme_bw()+xlab(NULL)+
    geom_text(data = linedf,aes(x=xi-0.1*x,y=0.5*y,label=xi))+
    geom_vline(data = linedf,aes(xintercept=xi),linetype=2,col="red")+theme(legend.position = "none")
  print(p)
  print(paste("recommend r_threshold: ",mean(linedf$xi)))
  return(mean(linedf$xi))
}

nnsd_plot <- \(histo = histo, title = title, threshold = threshold,dis_GOE = dis_GOE, dis_possion = dis_possion,p_ks_test=p_ks_test) {
  plot(histo, freq = F, col = "#F4FCA1", main = title,font.main = 1, xlab = "eigenvalue spacing", ylab = "PDF of eigenvalue spacing")
  {
    actual.ymax = par("yaxp")[2]
    x0 = -log(actual.ymax * 0.98)
    possion_dis=\(x)exp(-x)
    curve(possion_dis, from = max(x0, min(histo$breaks)),
          to = max(histo$breaks), n = 1001, add = T, type = "l", lty = 1, col = "#EB34FF", lwd = 2)
  }
  {
    GOE=function (x)pi/2 * x * exp(-pi/4 * x^2)
    curve(GOE, from = min(histo$breaks),
          to = max(histo$breaks), n = 1001, add = T, type = "l",
          lty = 1, col = "blue", lwd = 2)
  }

  if ((!is.na(dis_GOE)) && (!is.na(dis_possion))) {
    mtext(side = 3, paste("Distance to GOE =",signif(dis_GOE, 3),
                          "\nDistance to Possion =", signif(dis_possion, 3),"; ks_test p.value for possion =",signif(p_ks_test, 3)), col = "#878787", cex = 0.6)}

  if (!is.na(threshold))mtext(side = 4, paste("threshold =", signif(threshold,4)))

  legend("topright", inset = 0.05, c("Possion", "GOE"), col = c("#EB34FF", "blue"), lty = 1, lwd = 2, cex = 0.8)
}
trapez=\(x,y){ind = 2:length(x);as.double((x[ind] - x[ind - 1]) %*% (y[ind] + y[ind - 1]))/2}
get_sse=\(ev.spacing){
  dens = density(ev.spacing)
  N=20
  x = seq(min(ev.spacing), max(ev.spacing), len = 1000)
  observed = approx(dens$x, dens$y, xout = x)$y
  expected = exp(-x)
  A = exp(-min(ev.spacing)) - exp(-max(ev.spacing))
  xs <- numeric(N + 1)
  xs[1] = min(ev.spacing)
  for (i in 1:N) xs[i + 1] = -log(exp(-xs[i]) - A/N)
  area = numeric(N)
  for (i in 1:N) {
    xsec = x[(x > xs[i]) & (x < xs[i + 1])]
    xsec = c(xs[i], xsec, xs[i + 1])
    ysec = approx(dens$x, dens$y, xout = xsec)$y
    area[i] = trapez(xsec, ysec)
  }
  sse = sum((area[i] - A/N)^2)
  sse
}

#' Get RMT threshold for a correlation matrix roughly
#'
#' @export
#'
#' @rdname RMT_threshold
rmt = function(occor.r, min_threshold = 0.5, max_threshold = 0.85, step = 0.01){
  lib_ps("ggplot2")
  NNSD=\(x)abs(diff(x))

  s <- seq(0, 3, 0.1)
  poisson_d <- exp(-s)
  nnsdpois <- density(NNSD(poisson_d))

  ps <- c()
  threshold <- c()

  for(i in seq(min_threshold, max_threshold, step)){
    corr_r1 <- occor.r
    corr_r1[abs(corr_r1) < i] <- 0

    {
      eigen_res <- sort(eigen(corr_r1)$value)
      #spline to eigen_res
      check <- tryCatch(ssp <- smooth.spline(eigen_res, control.spar = list(low = 0, high = 3)),
                        error = \(e) {TRUE})
      if(rlang::is_true(check))next
      nnsdw <- density(NNSD(ssp$y))
      chival <- sum((nnsdw$y - nnsdpois$y)^2/nnsdpois$y/1e3)
    }

    ps <- c(ps, chival)
    threshold <- c(threshold, i)
    if(((i*100) %% 5==0)){
      print(paste0("Calculating: ",i))
    }
  }

  res <- data.frame(threshold,ps)
  recommend_thres <- res[which.min(res[, 2]), 1]

  p=ggplot(res,aes(threshold,ps))+geom_point()+
    geom_vline(xintercept = recommend_thres,linetype=2,col="red")+
    geom_text(x=recommend_thres+0.01,y=0.5*max(res$ps),label=recommend_thres)+theme_bw(base_size = 15)
  print(p)
  res1=res[(res$threshold<(recommend_thres+0.05))&(res$threshold>(recommend_thres-0.05)),]
  p=ggplot(res1,aes(threshold,ps))+geom_point()+
    geom_vline(xintercept = recommend_thres,linetype=2,col="red")+
    geom_text(x=recommend_thres+0.01,y=0.5*max(res1$ps),label=recommend_thres)+theme_bw(base_size = 15)
  print(p)

  print(paste0("We recommend r-threshold: ",recommend_thres,", you can calculate again in a smaller region"))
  recommend_thres
}
