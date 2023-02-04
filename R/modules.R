#========6.modules=========
#' Detect the modules
#'
#' @param go a igraph object
#' @param method cluster_method:("cluster_walktrap","cluster_edge_betweenness","cluster_fast_greedy","cluster_spinglass")
#' @return a igraph object
#' @export
#'
#' @examples
#' data("c_net")
#' modu_dect(co_net) -> co_net_modu
modu_dect <- function(go, method = "cluster_fast_greedy") {
  lib_ps("igraph")
  stopifnot(is_igraph(go))
  ms <- c("cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass")
  if (!method %in% ms) stop("method should be one of ", paste(ms, collapse = ","))

  if (method == "cluster_walktrap") {
    wc <- igraph::cluster_walktrap(go, weights = abs(igraph::E(go)$weight))
  }
  if (method == "cluster_edge_betweenness") {
    wc <- igraph::cluster_edge_betweenness(go, weights = abs(igraph::E(go)$weight))
  }
  if (method == "cluster_fast_greedy") {
    wc <- igraph::cluster_fast_greedy(go, weights = abs(igraph::E(go)$weight), )
  }
  if (method == "cluster_spinglass") {
    wc <- igraph::cluster_spinglass(go, weights = abs(igraph::E(go)$weight))
  }
  V(go)$module <- membership(wc)

  graph_attr(go)$modularity<- modularity(wc)

  rand.g <- erdos.renyi.game(length(V(go)), length(E(go)),type = "gnm")
  rand_m=modularity(cluster_fast_greedy(rand.g))
  r_modularity=(modularity(wc)-rand_m)/rand_m #relative modularity
  graph_attr(go)$r_modularity=r_modularity
  return(go = go)
}

#' Zi-Pi calculate
#'
#' @param go igraph object after modu_dect()
#' @param mode use 7-group (mode=1) or 4-group (mode=2)
#'
#' @return igraph
#' @export
#'
#' @references 1. Guimerà, R. & Amaral, L. Functional cartography of complex metabolic networks. (2005) doi:10.1038/nature03288.
#' @examples
#' data("c_net")
#' modu_dect(co_net) -> co_net_modu
#' zp_analyse(co_net_modu,mode = 2)->co_net_modu
#' zp_plot(co_net_modu)
zp_analyse<-function(go,mode=1){
  v_index=igraph::as_data_frame(go,"vertices")
  if(!"module"%in%names(v_index))stop("no modules, please modu_dect() first")
  if("roles"%in%names(v_index))message("areadly has roles, overwrite!")
  within <- within_module_deg_z_score(go)
  v_index$Ki=within$Ki;v_index$Zi=within$Zi
  pc <- part_coeff(go)
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
      if(between(x,backs$x1[i],backs$x2[i])&&between(y,backs$y1[i],backs$y2[i])){
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

#' Plot a community (modularity)
#'
#' @param wc a community object or a list contain wc and go
#' @param go a igraph object
#' @param n_modu >n_modu modules will be ploted
#' @param ... additional
#' @return a plot
#' @export
#'
#' @examples
#' data("c_net")
#' modu_dect(co_net) -> co_net_modu
#' modu_plot(co_net_modu, n_modu = 20)
modu_plot <- function(go, coors=NULL,n_modu = 1, ...) {
  if(!"module"%in%vertex_attr_names(go))stop("no modules, please modu_dect() first")
  as_group <- \(x){
    x <- list(membership = x)
    vids <- names(x$membership)
    modus <- tapply(vids, x$membership, simplify = FALSE, function(x) x)
    return(modus)
  }

  set.seed(123)
  if (is.null(coors)) coors <- layout.fruchterman.reingold(go)
  else {
    coors<-coors[match(V(go)$name,coors$name),2:3]%>%as.matrix()
  }
  members=V(go)$module;names(members)=V(go)$name

  if (n_modu) {
    # 筛选成员大于n_modu的模块
    (which((members %>% table()) > n_modu))->cho_modu
    members[members %in%cho_modu ] %>% as_group() -> new_modu
    col=ifelse(members %in%cho_modu,members ,"grey")

    plot(go,...,
         main = "Modularity network", layout = coors,
         mark.groups = new_modu, vertex.color = col,
         vertex.label.font = 2, vertex.label.color = "black",
         vertex.label.cex = .05* V(go)$size, vertex.size = V(go)$size,
         vertex.label = ifelse(V(go)$size > 8, V(go)$name, NA),
         edge.lty = 1, edge.curved = TRUE, margin = c(0, 0, 0, 0)
    )
    leg1<-paste(names(new_modu),sapply(new_modu, length),"nodes",sep = " :")
    legend(1, 1, cex = .7, legend =leg1 ,
           fill = rainbow(length(new_modu),alpha = 0.3), bty = "n", title = "Modules")

  }
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
#'
#' @return a ggplot object
#' @export
#' @import ggpubr ggrepel
#' @seealso \link{zp_analyse}
zp_plot <- function(go,label=T) {
  lib_ps("ggrepel","ggpubr","igraph")
  as.data.frame(vertex.attributes(go))->taxa.roles
  if(!"roles"%in%names(taxa.roles))stop("no roles, please zp_analyse() first")
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
    geom_point(data = taxa.roles, aes(x = Pi, y = Zi, color = factor(class))) +
    ggpubr::theme_pubr(legend = "right") +
    guides(colour = "none") +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Participation Coefficient") +
    ylab("Within-module connectivity z-score")

  if(label){
    p=p+ggrepel::geom_text_repel(
      data = filter(taxa.roles, !roles%in%c("Peripherals","Ultra-peripherals")),
      aes(x = Pi, y = Zi, label = name)
    )
  }
  return(p)
}
