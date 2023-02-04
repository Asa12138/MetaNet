#=========2.build======
#' Construct a network for table corr
#'
#' @param r_thres r_threshold (default:>0.6)
#' @param p_thres p_threshold (default:<0.05)
#' @param del_single should delete single vertexes?
#' @param corr result from c_net_cal()
#'
#' @return an igraph object
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab) -> totu
#' metadata[,3:10] -> env
#' c_net_cal(totu) -> corr
#' c_net_build(corr,r_thres = 0.63) -> co_net
#' c_net_cal(totu,env) -> corr2
#' c_net_build(corr2) -> co_net2
c_net_build <- function(corr, r_thres = 0.6, p_thres = 0.05, del_single = T) {
  suppressMessages(lib_ps("igraph"))
  # set thresholds to construct
  occor.r <- corr$r
  occor.p <- corr$p.value
  occor.r[occor.p > p_thres | abs(occor.r) < r_thres] <- 0
  corr <- occor.r
  # make igraph
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
    go <- graph_from_adjacency_matrix(as.matrix(corr), mode = "undirected", weighted = T, diag = F)
  } else {
    go <- graph_from_incidence_matrix(as.matrix(corr), directed = F, weighted = T)
  }
  # delete single vertexes?
  if (del_single) go <- igraph::delete.vertices(go, V(go)[igraph::degree(go) == 0])
  # set vertex attributes
  # set vertices shape
  V(go)$group <- ifelse(V(go)$name %in% rownames(corr), "group1", "group2")
  V(go)$shape <- ifelse(V(go)$group == "group1", "circle", "square")
  V(go)$color <- ifelse(V(go)$group == "group1", "#A6CEE3","#B15928")
  V(go)$class <- ifelse(V(go)$group == "group1", "group1", "group2")
  V(go)$size <- 5
  # abs edges weight
  go.weight <- E(go)$weight
  E(go)$cor <- go.weight
  E(go)$weight <- abs(go.weight)
  # set edges color
  E(go)$color <- ifelse(go.weight > 0, "#48A4F0", "#E85D5D")
  E(go)$inter<-ifelse(go.weight > 0, "positive", "negative")
  # set edges width
  E(go)$width <- 1 * E(go)$weight
  return(go)
}


#' Get RMT threshold for a correlation matrix
#'
#' @param occor.r a correlation matrix
#'
#' @return a r-threshold
#' @export
#' @references 1.  J. Zhou, Y. Deng, F. Luo, Z. He, Q. Tu, X. Zhi, Functional Molecular Ecological Networks (2010), doi:10.1128/mBio.00169-10.
#' @examples
#' t(otutab) -> totu
#' c_net_cal(totu) -> corr
#' RMT_threshold(corr$r)
#' #0.65
RMT_threshold<-function(occor.r){
  #Random matrix theory
  lib_ps("RMThreshold")
  #Get threshold by RMT
  nwd=getwd()
  if(!dir.exists("./RMT_temp"))dir.create("./RMT_temp")
  setwd("./RMT_temp")

  RMThreshold::rm.matrix.validation(occor.r,unfold.method = "spline")

  res <- RMThreshold::rm.get.threshold(occor.r)
  print(paste("The Intermediate files saved in ",getwd()))
  setwd(nwd)
  #thre <- (res$sse.chosen + res$p.ks.chosen)/2
  thre <-res$chosen.thresholds
  thre
}


#' Use dataframe to annotate vertexes of a igraph
#'
#' @param go a igraph object
#' @param anno_tab a dataframe using to annotate(with rowname or a name column)
#' @return a annotated igraph object
#' @export
#'
#' @examples
#' data("c_net")
#' data("otutab")
#' anno_vertex(co_net, taxonomy)
anno_vertex <- function(go, anno_tab) {
  lib_ps("igraph")
  vertex.attributes(go) %>% as.data.frame() -> v_atr
  if (!"name" %in% colnames(anno_tab)) rownames(anno_tab) -> anno_tab$name
  v_atr <- dplyr::left_join(v_atr, anno_tab, by = "name", suffix = c(".x", ""))
  v_atr %>% dplyr::select(!ends_with(".x")) -> v_atr
  as.list(v_atr) -> vertex.attributes(go)
  return(go)
}


#' Set basic attributes from totu table
#'
#' @param go a igraph object
#' @param totu t(otu table)
#' @param anno_col a dataframe using to annotate(with rowname or a name column)
#' @param totu2 t(otu table) if it exist
#' @param anno_col2 a dataframe using to annotate(with rowname or a name column)
#'
#' @import  RColorBrewer
#' @return a igraph object
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab) -> totu
#' metadata[,3:10] -> env
#'
#' data("c_net")
#' co_net <- c_net_set(co_net, totu, taxonomy %>% dplyr::select(Phylum))
#'
#' c_net_set(co_net2, totu, taxonomy %>% dplyr::select(Phylum), env) -> co_net2
#'
c_net_set <- function(go, totu = NULL, anno_col = NULL, totu2 = NULL, anno_col2 = NULL) {
  lib_ps("igraph")
  # set vertex attributes
  # set vertices shape
  V(go)$group <- ifelse(V(go)$name %in% colnames(totu), "group1", "group2")
  V(go)$shape <- ifelse(V(go)$group == "group1", "circle", "square")

  # set vertices size
  if (!is.null(totu)) {
    otu_pro <- totu %>%
      colSums() %>%
      data.frame(name = names(.), abundant = ., size = mmscale(., 3, 12))
    if (!is.null(totu2)) {
      otu_pro2 <- totu2 %>%
        colSums() %>%
        data.frame(name = names(.), abundant = ., size = mmscale(., 3, 12))
      otu_pro <- rbind(otu_pro, otu_pro2)
    }
    anno_vertex(go, otu_pro) -> go
  }
  # set vertices color
  if (!is.null(anno_col)) {
    # stopifnot(ncol(anno_col)==1)
    anno_col <- data.frame(name = rownames(anno_col), class = anno_col[, 1])
    if (!is.null(anno_col2)) {
      # stopifnot(ncol(anno_col2)==1)
      anno_col <- data.frame(name = rownames(anno_col2), class = anno_col2[, 1])
      anno_col <- rbind(anno_col, anno_col2)
    }
    if (is.numeric(anno_col$class)) stop("anno is numeric!!!")
    go.col <- droplevels(as.factor(anno_col$class))
    levels(go.col) <- colorRampPalette(RColorBrewer::brewer.pal(7, "Dark2"))(nlevels(go.col))

    anno_col$color <- as.character(go.col)
    anno_vertex(go, anno_col) -> go
  }
  # fill NAs
  V(go)$color <- ifelse(is.na(V(go)$color), "#B15928", V(go)$color)
  V(go)$class <- ifelse(is.na(V(go)$class), "group2", V(go)$class)
  V(go)$size <- ifelse(is.na(V(go)$size), 5, V(go)$size)
  # set edge intra-inter
  same <- \(x){
    return(all(x == x[1]))
  }
  E(go)$class <- get.edgelist(go) %>% apply(., 1, \(x)ifelse(same(x %in% colnames(totu)), "intra", "inter"))
  if (length(E(go)$class %>% unique()) > 1) {
    E(go)$lty <- ifelse(E(go)$class == "intra", 1, 5)
  } else {
    E(go)$lty <- 1
  }

  return(go)
}
