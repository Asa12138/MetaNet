#========5.topological=======
#' Extract sub-nwtwork from the whole network
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
#' c_net_cal(t(otutab))->corr
#' corr$p.value<-p.adjust.table(corr$p.value)
#' c_net_build(corr)->a_net
#' extract_sub_net(a_net,otutab,save_net = "testnet")
extract_sub_net<-function(a_net,otutab,threads=4,save_net=F){
  lib_ps("igraph")
  V(a_net)$name->v_name
  reps=ncol(otutab)

  lib_ps("parallel","doSNOW")
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  sub_nets <- foreach(i = 1:reps,.packages = c("igraph")) %dopar% {
    rownames(otutab)[otutab[,i]>0]->exist_sp
    induced.subgraph(a_net,which(v_name%in%exist_sp),impl = "copy_and_delete")->spe_sub
  }
  names(sub_nets)<-colnames(otutab)

  sub_net_pars <- foreach(i = 1:reps,.options.snow = opts,
                          .packages = c("igraph","pctax","pcutils"),.combine = "rbind") %dopar% {
                            sub_nets[[i]]->spe_sub
                            net_par(spe_sub,mode = "n")[["n_index"]]->indexs
                            wc <- igraph::cluster_fast_greedy(spe_sub, weights = abs(igraph::E(spe_sub)$weight), )
                            indexs$modularity <- modularity(wc)
                            indexs
                          }

  stopCluster(cl)
  rownames(sub_net_pars)<-colnames(otutab)
  if(is.character(save_net))save(sub_nets,file = paste0(save_net,".rda"))
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


#' Calculate all indexs of a network
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
net_par <- function(go, mode = c("v", "e", "n", "all"),fast=T) {
  pcutils::lib_ps("igraph","dplyr")
  stopifnot(is_igraph(go))
  if ("all" %in% mode) mode <- c("v", "e", "n")

  n_index <- NULL
  v_index <- NULL
  e_index <- NULL
  #non-weighted network
  up <- go
  if(!is.null(igraph::edge_attr(up)[["weight"]]))up<-delete_edge_attr(up,"weight")
  if ("n" %in% mode) {
    # Calculate Network Parameters
    n_index <- data.frame(
      num_nodes = length(V(go)), # number of nodes
      num_edges = length(E(go)), # number od edges
      edge_density = edge_density(go), # density of network, connectance
      neg_percent=sum(E(go)$cor<0)/length(E(go)), # negative edges percentage
      ave_path_len = average.path.length(up), # Average path length
      #w_ave_path_len = ifelse(is.null(E(go)$weight), ave_path_len, average.path.length(go)) # weighted Average path length
      global_efficiency=igraph::global_efficiency(up),
      ave_degree = mean(igraph::degree(go)), # Average degree
      w_ave_degree = ifelse(is.null(E(go)$weight), mean(igraph::degree(go)), sum(E(go)$weight) / length(V(go))), # weighted degree
      diameter = diameter(up), # network diameter
      clusteringC = transitivity(go), # Clustering coefficient
      cen_betweenness = centralization.betweenness(go)$centralization, # Betweenness centralization
      nat_connectivity = nc(go) # natural
    )

    if(!fast){
      # mean_dist=mean_distance(go)#
      # w_mean_dist=ifelse(is.null(E(go)$weight),mean_dist,mean_distance(go))
      # v_conn= vertex.connectivity(go) #
      # e_conn= edge.connectivity(go) #
      # components= count_components(go) #
      modularity=modularity(cluster_fast_greedy(go))#
      rand.g <- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
      rand_m=modularity(cluster_fast_greedy(rand.g))
      r_modularity=(modularity-rand_m)/rand_m #

      n_index <- data.frame(
        n_index,
        modularity=modularity,
        r_modularity=r_modularity,
        cen_closeness = centralization.closeness(go)$centralization, # Closeness centralization
        cen_degree = centralization.degree(go)$centralization, # Degree centralization
        cen_evcent = centralization.evcent(go)$centralization # eigenvector centralization
      )
    }

    n_index <- apply(n_index, 1, FUN = \(x)replace(x, is.nan(x), 0)) %>%
      t() %>%
      as.data.frame()
    if (!(graph_attr(go) %>% unlist() %>% is.null())) n_index <- data.frame(graph_attr(go), n_index)
  }
  if ("v" %in% mode) {
    # Calculate Vertices Parameters
    v_index <- data.frame(
      degree = degree(go),
      clusteringC = transitivity(go, type = "local"), # local clustering coefficient
      betweenness = betweenness(go), # betweenness
      eccentricity = eccentricity(go),
      closeness=closeness(go),
      hub_score=hub_score(go)[["vector"]]
      #page_rank = page.rank(go)$vector
      #igraph::evcent(g)[["vector"]]
      #igraph::local_efficiency(go)
    )

    if(!is.null(E(go)$cor)){
      igraph::as_data_frame(go)->edge_list
      edge_list%>%select(from,cor)%>%rbind(.,select(edge_list,to,cor)%>%rename(from=to))%>%
        group_by(from)%>%summarise(w_degree=sum(cor))->w_degree
      v_index$w_degree=w_degree[match(rownames(v_index),w_degree$from),"w_degree"]%>%unlist()
    }

    v_index <- apply(v_index, 1, FUN = \(x)replace(x, is.nan(x), 0)) %>%
      t() %>%
      as.data.frame()
    if (!(vertex_attr(go) %>% unlist() %>% is.null())) v_index <- data.frame(vertex_attr(go), v_index)
  }
  if ("e" %in% mode) {
    # Calculate Edges Parameters
    e_index <- data.frame(
      igraph::as_data_frame(go)
    )
    # if(!(edge_attr(go)%>%unlist()%>%is.null()))e_index=data.frame(edge_attr(go),e_index)

  }
  return(list(n_index = n_index, v_index = v_index, e_index = e_index))
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
fit_power<-function(go){
  #degree distribution
  degree_dist <- table(degree(go))
  dat <- data.frame(degree = as.numeric(names(degree_dist)), count = as.numeric(degree_dist))
  #fit, set the original a & b
  mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
  summary(mod)
  #extract the coefficient
  a <- round(coef(mod)[1], 3)
  b <- round(coef(mod)[2], 3)
  fit <- fitted(mod)
  SSre <- sum((dat$count-fit)^2)
  SStot <- sum((dat$count-mean(dat$count))^2)
  R2 <- round(1 - SSre/SStot, 3)

  #bootstrap t get pvalue
  p_num <- 1
  dat_rand <- dat
  for (i in 1:999) {
    dat_rand$count <- sample(dat_rand$count)
    SSre_rand <- sum((dat_rand$count-fit)^2)
    SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
    R2_rand <- 1 - SSre_rand/SStot_rand
    if (R2_rand > R2) p_num <- p_num + 1
  }
  p_value <- p_num / (999+1)

  p <- ggplot(dat, aes(x = degree, y = count)) +
    geom_point() +theme_bw() +
    stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
    labs(x = 'Degree', y = 'Count')

  label <- data.frame(
    x=0.8*max(dat$degree),
    y=c(0.9,0.85,0.8)*max(dat$count),
    formula = c(sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
                sprintf('italic(R^2) == %.3f', R2),
                sprintf('italic(P) < %.3f', p_value))
  )

  p + geom_text(aes(x=x,y=y,label = formula), data = label,parse = T)
}

###common_characteristic
#A network can be said "smallworld" if its smallworldness is higher than one (a stricter rule is smallworldness>=3; Humphries & Gurney, 2008)
# qgraph::smallworldness(g,B = 5)
# qgraph::smallworldIndex(g)
if(F){
  # degree_distribution,KS.p>0.05 indicated fit power-law distribution.
  fit_power_law(degree(go)+1,10)
  rand.g<- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
  fit_power_law(degree(rand.g)+1,10)
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
    x=log(degree(go)))->cc_k
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
#'
#' rand_net(co_net)
#' rand_net_par(co_net)->a
rand_net<-function(go = go){
  lib_ps("igraph")
  #generate a random network
  rand.g<- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")

  data1 = data.frame(freq= degree_distribution(go),net = "network",degree = 1:length(degree_distribution(go)))
  data2 = data.frame(freq = degree_distribution(rand.g) ,net = "random E-R",degree = 1:length(degree_distribution(rand.g)))
  data = rbind(data1,data2)
  p1 <- ggplot(data)+
    geom_point(aes(x = degree,y = freq,group =net,fill = net),pch = 21,size = 2) +
    geom_smooth(aes(x = degree,y = freq,group =net,color = net),se=F)+
    labs(x="Degree",y="Proportion")+scale_color_manual(values = c("#F58B8B","#7AADF0"))+
    scale_fill_manual(values = c("#F58B8B","#7AADF0"))+
    theme_pubr()+theme(legend.position = c(0.8,0.9),legend.title = element_blank())
  return(p1)
}

#' Net_pars of many random network
#'
#' @param go igraph
#' @param reps simulation time
#' @param threads threads
#'
#' @export
rand_net_par<-function(go,reps=99,threads=4){
  #parallel
  lib_ps("foreach","doSNOW")
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  rand_net_pars <- foreach(i = 1:reps,.options.snow = opts,
                           .packages = c("igraph"),.combine = "rbind") %dopar% {
                             #generate a random network
                             rand.g<- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
                             pctax::net_par(rand.g,mode = "n")[["n_index"]]->indexs
                             wc <- igraph::cluster_fast_greedy(rand.g)
                             indexs$modularity <- modularity(wc)
                             indexs
                           }
  stopCluster(cl)
  if(F){
    ggplot(a,aes(x=clusteringC))+geom_histogram()+
      geom_vline(xintercept = transitivity(go), col = 'red')
    ggplot(a,aes(x=ave_path_len))+geom_histogram()+
      geom_vline(xintercept = average.path.length(go), col = 'red')+xlim(2,3)
  }
  rand_net_pars
}

#' Calculate small-world coefficient
#'
#' @param go igraph
#'
#' @export
#' @examples
#' smallworldness(co_net)
smallworldness<-function(go){
  rand_net_par(go,reps = 99)->rands
  small_world_coefficient=(transitivity(go)/mean(rands$clusteringC))/
    (average.path.length(go)/mean(rands$ave_path_len))
  small_world_coefficient
}
