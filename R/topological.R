#========5.topological=======

#' Extract sub-network from the whole network
#'
#' @param a_net the whole network
#' @param otutab otutab, these columns will be extract
#' @param threads threads, default:4
#' @param save_net should save these sub_nets? F or a filename
#'
#' @return a dataframe contains all sub_net parameters
#' @export
#'
#' @examples
#' data(otutab)
#' extract_sub_net(co_net,otutab,save_net = "sub_net")->sub_net_pars
extract_sub_net<-function(a_net,otutab,threads=1,save_net=F){
  lib_ps("igraph")
  V(a_net)$name->v_name
  reps=ncol(otutab)

  print("extracting")
  sub_nets=lapply(1:reps, \(i){
    rownames(otutab)[otutab[,i]>0]->exist_sp
    subgraph(a_net,which(v_name%in%exist_sp))->spe_sub
    return(spe_sub)
  })
  names(sub_nets)<-colnames(otutab)
  print("calculating topological indexes")
  #parallel
  #main function
  loop=function(i){
    sub_nets[[i]]->spe_sub
    MetaNet::net_par(spe_sub,mode = "n")[["n_index"]]->indexs
    wc <- igraph::cluster_fast_greedy(spe_sub, weights = abs(igraph::E(spe_sub)$weight), )
    indexs$modularity <- igraph::modularity(wc)
    indexs
  }

  {if(threads>1){
    lib_ps("foreach","doSNOW")
    pb <- utils::txtProgressBar(max =reps, style = 3)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
    cl <- snow::makeCluster(threads)
    doSNOW::registerDoSNOW(cl)
    sub_net_pars <- foreach::foreach(i = 1:reps,.options.snow = opts) %dopar% {
                              loop(i)
                            }
    snow::stopCluster(cl)
    gc()
  }
    else {
      sub_net_pars <-lapply(1:reps, loop)
    }}
  #simplify method
  sub_net_pars=do.call(rbind,sub_net_pars)

  rownames(sub_net_pars)<-colnames(otutab)
  save(sub_nets,file = paste0(save_net,".rda"))

  sub_net_pars
}

#' Calculate natural_connectivity
#'
#' @param p a igraph object
#' @return natural_connectivity (numeric)
#' @export
#' @references \code{\link[ggClusterNet]{nc}}
#' @examples
#' make_ring(10) %>% nc()
nc <- function(p) {
  lib_ps("igraph")
  adj_matrix <- as.matrix(as_adj(p, sparse = F))
  adj_matrix[abs(adj_matrix) != 0] <- 1

  lambda <- eigen(adj_matrix, only.values = TRUE)$values
  lambda <- sort(lambda, decreasing = TRUE)

  lambda_sum <- 0
  N <- length(lambda)
  for (i in 1:N) lambda_sum <- lambda_sum + exp(lambda[i])
  lambda_average <- log(lambda_sum / N, base = exp(1))
  lambda_average
}


#' Calculate all indexes of a network
#'
#' @param go a igraph object
#' @param fast less indexes for faster calculate ?
#' @param mode calculate what? c("v", "e", "n", "all")
#'
#' @return a 3-elements list
#' \item{n_index}{indexs of the whole network}
#' \item{v_index}{indexs of each vertex}
#' \item{e_index}{indexs of each edge}
#' @export
#'
#' @examples
#' make_graph("Walther") %>% net_par()
#' get_net_par(co_net)->co_net_with_par
net_par <- function(go, mode = c("v", "e", "n", "all"),fast=T) {
  lib_ps("igraph","dplyr")
  stopifnot(is_igraph(go))
  if ("all" %in% mode) mode <- c("v", "e", "n")

  n_index <- NULL
  v_index <- NULL
  e_index <- NULL
  #non-weighted network
  up <- go
  if(!is.null(igraph::edge_attr(up)[["weight"]]))up<-igraph::delete_edge_attr(up,"weight")
  if ("n" %in% mode) {
    # Calculate Network Parameters
    n_index <- data.frame(
      num_nodes = length(igraph::V(go)), # number of nodes
      num_edges = length(igraph::E(go)), # number of edges
      edge_density = igraph::edge_density(go), # density of network, connectance
      neg_percent=ifelse(!is.null(E(go)$cor),sum(igraph::E(go)$cor<0)/length(igraph::E(go)),NA), # negative edges percentage
      ave_path_len = igraph::average.path.length(up), # Average path length
      #w_ave_path_len = ifelse(is.null(E(go)$weight), ave_path_len, average.path.length(go)) # weighted Average path length
      global_efficiency=igraph::global_efficiency(up),
      ave_degree = mean(igraph::degree(go)), # Average degree
      w_ave_degree = ifelse(is.null(igraph::E(go)$weight), mean(igraph::degree(go)), sum(igraph::E(go)$weight) / length(igraph::V(go))), # weighted degree
      diameter = igraph::diameter(up), # network diameter
      clusteringC = igraph::transitivity(go), # Clustering coefficient
      cen_betweenness = igraph::centralization.betweenness(go)$centralization, # Betweenness centralization
      nat_connectivity = MetaNet::nc(go) # natural
    )

    if(!fast){
      # mean_dist=mean_distance(go)#
      # w_mean_dist=ifelse(is.null(E(go)$weight),mean_dist,mean_distance(go))
      # v_conn= vertex.connectivity(go) #
      # e_conn= edge.connectivity(go) #
      # components= count_components(go) #
      modularity=igraph::modularity(igraph::cluster_fast_greedy(go))#
      rand.g <- igraph::erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
      rand_m=igraph::modularity(igraph::cluster_fast_greedy(rand.g))
      relative_modularity=(modularity-rand_m)/rand_m #

      n_index <- data.frame(
        n_index,
        modularity=modularity,
        relative_modularity=relative_modularity,
        cen_closeness = igraph::centralization.closeness(go)$centralization, # Closeness centralization
        cen_degree = igraph::centralization.degree(go)$centralization, # Degree centralization
        cen_evcent = igraph::centralization.evcent(go)$centralization # eigenvector centralization
      )
    }
    n_index <- apply(n_index, 1, FUN = \(x)replace(x, is.nan(x), 0)) %>%t() %>%as.data.frame()
    n_index <- cbind_new(get_n(go), n_index)
  }
  if ("v" %in% mode) {
    # Calculate Vertices Parameters
    v_index <- data.frame(
      degree = igraph::degree(go),
      clusteringC = igraph::transitivity(go, type = "local"), # local clustering coefficient
      betweenness = igraph::betweenness(go), # betweenness
      eccentricity = igraph::eccentricity(go),
      closeness=igraph::closeness(go),
      hub_score=igraph::hub_score(go)[["vector"]]
      #page_rank = page.rank(go)$vector
      #igraph::evcent(go)[["vector"]]
      #igraph::local_efficiency(go)
    )

    #weighted degree
    if(!is.null(E(go)$cor)){
      get_e(go)->edge_list
      edge_list%>%select(from,cor)%>%rbind(.,select(edge_list,to,cor)%>%rename(from=to))%>%
        group_by(from)%>%summarise(w_degree=sum(cor))->w_degree
      v_index$w_degree=w_degree[match(rownames(v_index),w_degree$from),"w_degree"]%>%unlist()
    }

    v_index <- apply(v_index, 1, FUN = \(x)replace(x, is.nan(x), 0)) %>%
      t() %>%as.data.frame()
    v_index <- cbind_new(get_v(go), v_index)
  }
  if ("e" %in% mode) {
    # Calculate Edges Parameters
    e_index <- get_e(go)
    # if(!(edge_attr(go)%>%unlist()%>%is.null()))e_index=data.frame(edge_attr(go),e_index)

  }
  return(list(n_index = n_index, v_index = v_index, e_index = e_index))
}

#' bind two df with same columns, the last df will replace first df
#' @export
#'
cbind_new<-\(df,df1){
  if(ncol(df)<1)return(df1)
  if(ncol(df1)<1)return(df)
  inter=intersect(colnames(df1),colnames(df))
  la=setdiff(colnames(df),inter)
  cbind(df[,la,drop=F],df1)
}

#'
#' @export
#'
#' @rdname net_par
get_net_par<-function(go,force=F){
  if(!force){if(!is.null(graph_attr(go)[["net_par"]]))stop("already calculate net_pars")}
  net_par(go,fast = F)->res
  graph_attr(go)<-as.list(res$n_index)
  graph_attr(go)[["net_par"]]=T
  vertex_attr(go)<-as.list(res$v_index)
  edge_attr(go)<-as.list(res$e_index)
  go
}

#' Get skeleton network according to a group
#'
#' @param go network
#' @param group vertex column name
#'
#' @return skeleton network
#' @export
#'
#' @examples
#' get_group_skeleton(co_net)->ske_net
#' skeleton_plot(ske_net)
#' x=tkplot(ske_net)
#' da <- tkplot.getcoords(x)
#' tkplot.close(x)
get_group_skeleton=function(go,Group ="v_class",count=NULL,direct=F){
  lib_ps("igraph")
  stopifnot(is_igraph(go))
  if(!Group%in%vertex_attr_names(go))stop("no Group named ",Group," !")
  get_v(go)->tmp_v
  tmp_v%>%dplyr::select(name,!!Group)->nodeGroup
  colnames(nodeGroup)=c("name","Group")
  nodeGroup$Group<-as.factor(nodeGroup$Group)
#summary edges counts in each e_type
  suppressMessages(anno_edge(go,nodeGroup)%>%get_e()->edge)
{  if(is.null(count))edge$count=1
  else edge$count=edge[,count]}
  bb=data.frame()
  for (i in unique(edge$e_type)) {
    bb=rbind(bb,data.frame(summ_2col(edge[edge$e_type==i,c("Group_from","Group_to","count")],direct = direct),e_type=i))
  }
  tmp_go=igraph::graph_from_data_frame(bb,directed = direct)
  nodeGroup=cbind_new(nodeGroup,data.frame(v_group=tmp_v$v_group))
  dplyr::distinct(nodeGroup,Group,v_group)%>%tibble::column_to_rownames("Group")->v_group_tab

  V(tmp_go)$v_group=v_group_tab[V(tmp_go)$name,"v_group"]
  V(tmp_go)$v_class=V(tmp_go)$name
  V(tmp_go)$size=aggregate(tmp_v$size,by=list(tmp_v[,Group]),sum)[["x"]]
  suppressWarnings({V(tmp_go)$count=tmp_v%>%dplyr::group_by_(Group)%>%dplyr::count()%>%pull(n)})

  tmp_go=c_net_update(tmp_go)
  get_e(tmp_go)->tmp_e

  E(tmp_go)$width=E(tmp_go)$label=tmp_e$count

  graph.attributes(tmp_go)$n_type="skeleton"
  graph.attributes(tmp_go)$skeleton=Group
  tmp_go
}

#' skeleton_plot
#'
#' @param tmp_go skeleton
#' @rdname get_group_skeleton
skeleton_plot<-function(tmp_go,...) {
  lib_ps("igraph")
  if(get_n(tmp_go)$n_type!="skeleton")stop("Not a skeleton network")
  get_e(tmp_go)->tmp_e
  get_v(tmp_go)->tmp_v
  c_net_lay(tmp_go)->coors

  #some
  params=list(...)
  params_name=names(params)
  legend_position_default=c(left_leg_x=-1.9,left_leg_y=1,right_leg_x=1.2,right_leg_y=1)
  if(!"legend_position"%in%params_name)legend_position=legend_position_default
  else {
    if(is.null(names(legend_position)))legend_position=setNames(legend_position,names(legend_position_default)[seq_along(legend_position)])
    legend_position=condance(data.frame(legend_position_default,tidai(names(legend_position_default),legend_position)))
  }
  if("vertex.color"%in%params_name)tmp_v$color=condance(data.frame(tmp_v$color,tidai(tmp_v$v_class,params[["vertex.color"]])))
  if("legend_cex"%in%params_name)legend_cex=parms[["legend_cex"]]
  else legend_cex=1

  for (i in unique(tmp_e$e_type)) {
    left_leg_x=legend_position["left_leg_x"]
    left_leg_y=legend_position["left_leg_y"]

    #main plot
    tmp_go1=c_net_filter(tmp_go,e_type=i)
    c_net_plot(tmp_go1,coors = coors,vertex.size=pcutils::mmscale(V(tmp_go1)$size,10,16),
               edge.width=pcutils::mmscale(E(tmp_go1)$width,0.5,8),col_legend = F,...,
               lty_legend = F,legend_number = F,size_legend = F)
    #color legend
    pchls=c("circle"=21,"square"=22)
    for(i in 1:length(unique(tmp_v$v_group))){
      g_i=unique(tmp_v$v_group)[i]
      tmp_v1=tmp_v[tmp_v$v_group==g_i,c("v_class","count","color","shape")]
      if(T){
        le_text= paste(tmp_v1$v_class,tmp_v1$count,sep = ": ")
      }
      else le_text= unique(tmp_v1$v_class)
      legend(left_leg_x, left_leg_y, cex = 0.7*legend_cex, adj = 0,
             legend = le_text,title.cex = 0.8*legend_cex,
             title =g_i,title.font = 2,title.adj = 0,
             col = "black", pt.bg = unique(tmp_v1$color), bty = "n", pch = pchls[unique(tmp_v1$shape)])
      left_leg_y=left_leg_y-(length(unique(tmp_v1$v_class))*0.12+0.2)*legend_cex
    }
  }
}

#' Link summary of the network
#'
#' @param go igraph
#' @param group summary which group of vertex attribution in names(vertex_attr(go))
#' @param inter "positive", "negative", "all"
#' @param topN topN of group, default:5
#'
#' @return plot
#' @export
#'
#' @examples
#' links_stat(co_net,topN=10)
#' modu_dect(co_net) -> co_net_modu
#' links_stat(co_net_modu,group="module")
#'
#' c_net_plot(co_net,g_lay_polyarc(co_net,"v_class"))
#'
links_stat<-function(go,group="v_class",inter="all",topN=6,direct=T,
                     legend_number=F,legend=T,legend_cex=1,
                     legend_position=c(left_leg_x=-1.6,left_leg_y=1,right_leg_x=1.2,right_leg_y=1),
                     col_legend_order=NULL,
                     group_legend_title=NULL,group_legend_order=NULL){
  pcutils::lib_ps("igraph","dplyr")
  go=c_net_set(go,vertex_class = group)
  get_v(go)->v_index
  v_index%>%select("name","v_class")->map

  suppressMessages(anno_edge(go,map)%>%get_e()->edge)
#statistics
  if(inter!="all")edge%>%filter(inter==!!inter)->edge
  summ_2col(edge[,paste0("v_class",c("_from","_to"))],direct = direct)->bb
  colnames(bb)=c("from","to","count")

  group_by(bb,from)%>%summarise(n=sum(count))%>%arrange(-n)%>%top_n(topN,n)%>%pull(from)->nnn
  filter(edge,v_class_from%in%nnn,v_class_to%in%nnn)%>%pull(id)->eid
  subgraph.edges(go,eid)->go1
  get_v(go1)->tmp_v
  colors=setNames(unique(tmp_v$color),unique(tmp_v$v_class))
#plot
  filter(bb,from%in%nnn,to%in%nnn)%>%pcutils::my_cicro(.,reorder = F,colors = colors)
#legend
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

  pchls=c("circle"=21,"square"=22)
  vgroups=change_fac_lev(tmp_v$v_group,group_legend_order)
  vgroups=levels(vgroups)
  if(is.null(group_legend_title)){
    if(length(vgroups)>1)group_legend_title=setNames(paste0(vgroups,"_",group),vgroups)
    else group_legend_title=setNames(group,vgroups)
    }
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
  }
}

#redundant，Deprecated
if(F){links_stat2<-function(go,group="v_class",inter="positive",topN=5){
    pcutils::lib_ps("igraph","dplyr")

    get_e(go)[,c("from","to","inter")]->edge
    if(inter!="all")edge%>%filter(inter==!!inter)->edge

    get_v(go)%>%select("name",!!group)->map
    edge=left_join(edge,map,by=c("from"="name"))%>%rename("f_g"=!!group)
    edge=left_join(edge,map,by=c("to"="name"))%>%rename("t_g"=!!group)

    c(edge$f_g,edge$t_g)%>%unique()->all_g
    expand.grid(all_g,all_g)->tab
    group_by(edge,f_g,t_g)%>%count()->link
    link=left_join(tab,link,by=c("Var1"="f_g","Var2"="t_g"))

    tab=reshape2::dcast(link, Var1~Var2, value.var = "n")%>%tibble::column_to_rownames("Var1")%>%as.matrix()
    tab[is.na(tab)]=0
    #tab=tab+t(tab)-diag(diag(tab))

    rowSums(tab)%>%sort(decreasing = T)%>%names()->s_name
    tab=tab[s_name,s_name]
    if(topN>ncol(tab))topN=ncol(tab)
    tab=tab[seq_len(topN),seq_len(topN)]

    pcutils::lib_ps("circlize")
    circlize::chordDiagram(tab,grid.col = pcutils::get_cols(topN))
    #chorddiag::chorddiag(tab)
    detach("package:circlize")
  }
}
#try synteny plot
if(F){
  #https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html

}

#' Fit power-law distribution for a igraph
#'
#' @param go igraph
#'
#' @return ggplot
#' @export
#'
#' @examples
#' fit_power(co_net)
fit_power<-function(go,pvalue=F){
  #igraph::degree distribution
  degree_dist <- table(igraph::degree(go))
  dat <- data.frame(degree = as.numeric(names(degree_dist)), count = as.numeric(degree_dist))
  #fit, set the original a & b
  mod <- stats::nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
  summary(mod)
  #extract the coefficient
  a <- round(coef(mod)[1], 3)
  b <- round(coef(mod)[2], 3)
  fit <- fitted(mod)
  SSre <- sum((dat$count-fit)^2)
  SStot <- sum((dat$count-mean(dat$count))^2)
  R2 <- round(1 - SSre/SStot, 3)

  #bootstrap t get pvalue
if(pvalue){  p_num <- 1
  dat_rand <- dat
  for (i in 1:999) {
    dat_rand$count <- sample(dat_rand$count)
    SSre_rand <- sum((dat_rand$count-fit)^2)
    SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
    R2_rand <- 1 - SSre_rand/SStot_rand
    if (R2_rand > R2) p_num <- p_num + 1
  }
  p_value <- p_num / (999+1)}

  p <- ggplot(dat, aes(x = degree, y = count)) +
    geom_point() +theme_bw() +
    stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
    labs(x = 'Degree', y = 'Count')

{if(pvalue){  label <- data.frame(
    x=0.8*max(dat$degree),
    y=c(0.9,0.85,0.8)*max(dat$count),
    formula = c(sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
                sprintf('italic(R^2) == %.3f', R2),
                sprintf('italic(P) < %.3f', p_value))
  )}
else{  label <- data.frame(
    x=0.8*max(dat$degree),
    y=c(0.9,0.85)*max(dat$count),
    formula = c(sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
                sprintf('italic(R^2) == %.3f', R2))
  )}}

  p + geom_text(aes(x=x,y=y,label = formula), data = label,parse = T)
}

###common_characteristic
#A network can be said "smallworld" if its smallworldness is higher than one (a stricter rule is smallworldness>=3; Humphries & Gurney, 2008)
# qgraph::smallworldness(g,B = 5)
# qgraph::smallworldIndex(g)
if(F){
  # degree_distribution,KS.p>0.05 indicated fit power-law distribution.
  fit_power_law(igraph::degree(go)+1,10)
  rand.g<- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
  fit_power_law(igraph::degree(rand.g)+1,10)
  fit_power(go)

  #smallworldness, the average path lengths which were close to logarithms of the total number of network nodes
  #生成具有 Nv=25 个节点，r=5 个邻居，重连概率 p=0.05 的网络
  g_rand_s <- watts.strogatz.game(dim = 1, size = 25, nei = 5, p = 0.05)
  transitivity(g_rand_s)
  average.path.length(g_rand_s)
  log10(length(V(g_rand_s)))
  #具有相同节点和边数量的经典随机图
  g_rand_ER <- erdos.renyi.game(n = vcount(g_rand_s), p = ecount(g_rand_s), type = 'gnm')
  transitivity(g_rand_ER)
  average.path.length(g_rand_ER)

  #modularity,values >0.4 suggest that the network has a modular structure; Newman, 2006
  modu_dect(go)%>%graph.attributes()

  #hierarchy, R2 values of the linear relationship between logarithms of clustering coefficients and the logarithms of connectivity
  data.frame(
    y=transitivity(go,"local"),
    x=log(igraph::degree(go)))->cc_k
  summary(lm(y~x,cc_k))

}


#' Degree distribution comparison with random network
#'
#' @param go igraph object
#'
#' @return ggplot
#' @export
#'
#' @examples
#' rand_net(co_net_rmt)
#' rand_net_par(co_net)->a
rand_net<-function(go = go){
  lib_ps("igraph")
  #generate a random network
  rand.g<- igraph::erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")

  data1 = data.frame(freq= igraph::degree_distribution(go),net = "Network",degree = 0:(length(degree_distribution(go))-1))

  data2 = data.frame(freq= igraph::degree_distribution(rand.g) ,net = "Random E-R",degree = 0:(length(degree_distribution(rand.g))-1))

  #if data1[1,1]=0, it'is delete single vertex
  if(data1[1,1]==0)data1=data1[-1,]
  data = rbind(data1,data2)
  p1 <- ggplot(data)+
    geom_point(aes(x = degree,y = freq,group =net,fill = net),pch = 21,size = 2) +
    geom_smooth(aes(x = degree,y = freq,group =net,color = net),se=F,method = 'loess',formula = 'y ~ x')+
    labs(x="Degree",y="Proportion")+scale_color_manual(values = c("#F58B8B","#7AADF0"))+
    scale_fill_manual(values = c("#F58B8B","#7AADF0"))+
    ggpubr::theme_pubr()+theme(legend.position = c(0.8,0.9),legend.title = element_blank())
  print(p1)
  return(rand.g)
}

#' Net_pars of many random network
#'
#' @param go igraph
#' @param reps simulation time
#' @param threads threads
#'
#' @export
#' @rdname compare_rand
rand_net_par<-function(go,reps=99,threads=1){
  #parallel
  #main function
  loop=function(i){
    #generate a random network
    rand.g<- igraph::erdos.renyi.game(length(igraph::V(go)), length(igraph::E(go)),type = "gnm")
    MetaNet::net_par(rand.g,mode = "n")[["n_index"]]->indexs
    wc <- igraph::cluster_fast_greedy(rand.g)
    indexs$modularity <- igraph::modularity(wc)
    indexs
  }
  {
    if(threads>1){
      lib_ps("foreach","doSNOW")
      pb <- utils::txtProgressBar(max =reps, style = 3)
      opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
      cl <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cl)
      res <- foreach::foreach(i = 1:reps,.options.snow = opts) %dopar% {
                                loop(i)
                              }
      snow::stopCluster(cl)
      gc()
    }
    else {
      res <-lapply(1:reps, loop)
    }}
  #simplify method
  rand_net_pars=do.call(rbind,res)

  if(F){
    ggplot(a,aes(x=clusteringC))+geom_histogram()+
      geom_vline(xintercept = transitivity(go), col = 'red')
    ggplot(a,aes(x=ave_path_len))+geom_histogram()+
      geom_vline(xintercept = average.path.length(go), col = 'red')+xlim(2,3)
  }
  rand_net_pars
}

#' Compare some indexes between your net with random networks
#'
#' @param pars your net pars resulted by net_pars()
#' @param randp random networks pars resulted by rand_net_par()
#' @param index compared indexes: "ave_path_len","clusteringC" or else
#'
#' @return ggplot
#' @export
#'
#' @examples
#' rand_net_par(co_net_rmt,reps = 30)->randp
#' net_par(co_net_rmt,fast = F)->pars
#' compare_rand(pars,randp)
compare_rand<-function(pars,randp,index=c("ave_path_len","clusteringC")){
  labss=t(pars$n_index[,index])%>%as.data.frame()
  rownames(labss)->labss$indexes
  pcutils::group_box(randp[,index])+
    geom_hline(data=labss,aes(yintercept = V1),linetype=2)+
    geom_text(data=labss,aes(x=1,y = V1*1.05,label=paste0("Net: ",round(V1,3))))+
    theme_classic()+
    theme(legend.position = "none",axis.text.x = element_blank())
}


#' Calculate small-world coefficient
#'
#' @param go igraph
#'
#' @export
#' @examples
#' smallworldness(co_net)
smallworldness<-function(go,reps=99,threads = 1){
  rand_net_par(go,reps = reps,threads = threads )->rands
  small_world_coefficient=(igraph::transitivity(go)/mean(rands$clusteringC))/
    (igraph::average.path.length(go)/mean(rands$ave_path_len))
  small_world_coefficient
}


