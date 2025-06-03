#' Extract ego-centric subnetwork with preserved class attributes
#'
#' Wrapper around \code{igraph::make_ego_graph()} that ensures output retains
#' "metanet" and "igraph" class structure. Supports single or multiple center nodes.
#'
#' @param graph An igraph object with potential "metanet" class
#' @param nodes Center node(s) for subnetwork extraction (vertex IDs or names)
#' @param order Integer specifying the order of neighbors to include
#' @param mode Character scalar, either "in", "out" or "all" for directed networks
#'
#' @return metanet
#' @aliases c_net_neighbors
#'
#' @examples
#' library(igraph)
#' c_net_plot(co_net)
#' c_net_plot(c_net_ego(co_net, "s__Kribbella_catacumbae"))
#' nodes <- c("s__Kribbella_catacumbae", "s__Verrucosispora_andamanensis")
#' c_net_plot(c_net_ego(co_net, nodes))
#' @export
c_net_ego <- function(graph, nodes, order = 1, mode = "all") {
  # Input validation
  if (!inherits(graph, "igraph")) {
    stop("Input must be an igraph object")
  }

  # Extract ego networks
  ego_graphs <- make_ego_graph(
    graph = graph,
    order = order,
    nodes = nodes,
    mode = mode
  )

  # Merge if multiple centers
  if (length(nodes) > 1) {
    subg <- suppressMessages(Reduce(c_net_union, ego_graphs))
  } else {
    subg <- ego_graphs[[1]]
  }

  # Force class structure
  class(subg) <- c("metanet", "igraph")

  return(subg)
}


#' Highlight specific nodes in a network
#'
#' Adds highlight markers to specified nodes and grays out non-highlighted nodes.
#' Preserves all existing vertex/edge attributes and class structure.
#'
#' @param graph An igraph/metatnet object
#' @param nodes Vector of node names to highlight
#' @param gray_color Color for non-highlighted nodes (default: "gray80")
#' @param edges a data.frame of edges to highlight, colnames must be "from" and "to"
#'
#' @return metanet
#'
#' @examples
#' par(mfrow = c(1, 3))
#' nodes <- c("s__Kribbella_catacumbae", "s__Verrucosispora_andamanensis")
#' nodes <- V(c_net_ego(co_net, nodes))$name
#' g_hl <- c_net_highlight(co_net, nodes = nodes)
#' plot(g_hl) # Highlighted nodes keep colors, others turn gray
#' get_e(co_net) %>% head(20) -> hl_edges
#' g_hl2 <- c_net_highlight(co_net, edges = hl_edges[, 2:3])
#' c_net_plot(g_hl2)
#' g_hl3 <- c_net_highlight(co_net, nodes = nodes, edges = hl_edges[, 2:3])
#' c_net_plot(g_hl3)
#' @export
c_net_highlight <- function(graph, nodes = NULL, edges = NULL, gray_color = "gray80") {
  # 验证输入
  if (!inherits(graph, "igraph")) stop("Input must be an igraph object")
  if (!is_metanet(graph)) graph <- as.metanet(graph)
  # 备份原始颜色（如果不存在则创建）
  if ("color_original" %in% vertex_attr_names(graph)) {
    V(graph)$color <- V(graph)$color_original
  } else {
    V(graph)$color_original <- V(graph)$color
  }
  if ("color_original" %in% edge_attr_names(graph)) {
    E(graph)$color <- E(graph)$color_original
  } else {
    E(graph)$color_original <- E(graph)$color
  }

  E(graph)$v_highlight_from <- E(graph)$v_highlight_to <- V(graph)$v_highlight <- FALSE

  graph2 <- graph
  if (!is.null(nodes)) {
    # 转换节点标识为统一格式
    nodes <- which(V(graph2)$name %in% nodes)
    if (length(nodes) == 0) {
      message("No nodes found in the network.")
    }
    # 添加高亮标记
    V(graph2)$v_highlight[nodes] <- TRUE
    # 处理相应边
    graph2 <- suppressMessages(anno_edge(graph2, get_v(graph2)[, c("name", "v_highlight")]))
  }

  graph3 <- graph
  if (!is.null(edges)) {
    stopifnot(is.data.frame(edges))
    stopifnot(identical(c("from", "to"), colnames(edges)))
    if (!igraph::is_directed(graph3)) {
      edges2 <- edges
      edges2$from <- edges$to
      edges2$to <- edges$from
      edges <- rbind(edges, edges2)
    }
    edges$v_highlight_from <- TRUE
    edges$v_highlight_to <- TRUE
    graph3 <- suppressMessages(anno_edge(graph3, edges))
    E(graph3)$v_highlight_from[is.na(E(graph3)$v_highlight_from)] <- FALSE
    E(graph3)$v_highlight_to[is.na(E(graph3)$v_highlight_to)] <- FALSE

    nodes <- c(E(graph3)$from[E(graph3)$v_highlight_from], E(graph3)$to[E(graph3)$v_highlight_to])
    nodes <- which(V(graph3)$name %in% nodes)
    V(graph3)$v_highlight[nodes] <- TRUE
  }

  V(graph)$v_highlight <- V(graph2)$v_highlight | V(graph3)$v_highlight
  E(graph)$v_highlight_from <- E(graph2)$v_highlight_from | E(graph3)$v_highlight_from
  E(graph)$v_highlight_to <- E(graph2)$v_highlight_to | E(graph3)$v_highlight_to

  # 修改非高亮节点颜色
  V(graph)$color[!V(graph)$v_highlight] <- gray_color
  E(graph)$color[!(E(graph)$v_highlight_from & E(graph)$v_highlight_to)] <- gray_color

  return(graph)
}


#' Batch drawing multiple network diagrams
#'
#' @param graph_ls a list containing igraph objects
#' @param nrow nrow
#' @param ncol ncol
#' @param multi_params_list a list of parameters for each network
#'
#' @return No value
#'
#' @examples
#' plot_multi_nets(list(co_net, co_net2),
#'   multi_params_list = list(
#'     list(vertex.color = "skyblue"),
#'     list(vertex.color = "green3")
#'   )
#' )
#' @export
plot_multi_nets <- function(graph_ls,
                            nrow = NULL,
                            ncol = NULL,
                            multi_params_list = NULL) {
  # 参数验证
  if (!all(sapply(graph_ls, inherits, "igraph"))) {
    stop("All elements in graph_ls must be igraph objects")
  }

  # 计算布局行列数
  n <- length(graph_ls)
  if (is.null(nrow) && is.null(ncol)) {
    nrow <- floor(sqrt(n))
    ncol <- ceiling(n / nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(n / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(n / nrow)
  }

  # 保存原始图形参数
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # 设置多图布局
  par(mfrow = c(nrow, ncol), mar = c(1, 1, 2, 1))

  if (length(multi_params_list) == 0) multi_params_list <- vector("list", n)
  if (length(multi_params_list) != n) multi_params_list <- rep(multi_params_list, n)

  # 批量绘图
  for (i in seq_along(graph_ls)) {
    go <- graph_ls[[i]] # 注意这里只能用go，属于是参数环境问题
    c_net_plot(go, params_list = multi_params_list[[i]])
  }
}

#' Calculate Similarity Between Two Graphs via Adjacency Matrices
#'
#' Computes the similarity between two igraph objects using their adjacency matrices.
#' Supports Frobenius norm-based similarity and cosine similarity.
#'
#' @param g1 An igraph object representing the first graph.
#' @param g2 An igraph object representing the second graph.
#' @param method A character string specifying the similarity method:
#'               "frobenius" (default) or "cosine".
#' @return A numeric value between 0 (no similarity) and 1 (identical graphs).
#' @export
#' @examples
#' library(igraph)
#' g1 <- graph_from_edgelist(matrix(c(1, 2, 2, 3), ncol = 2, byrow = TRUE), directed = FALSE)
#' g2 <- graph_from_edgelist(matrix(c(1, 2, 2, 4), ncol = 2, byrow = TRUE), directed = FALSE)
#' adjacency_similarity(g1, g2, method = "frobenius") # Output: 0.5
#' adjacency_similarity(g1, g2, method = "cosine") # Output: 0.5
adjacency_similarity <- function(g1, g2, method = "frobenius") {
  if (!is_metanet(g1)) g1 <- as.metanet(g1)
  if (!is_metanet(g2)) g2 <- as.metanet(g2)
  # 获取邻接矩阵
  adj1 <- as.matrix(igraph::as_adjacency_matrix(g1))
  adj2 <- as.matrix(igraph::as_adjacency_matrix(g2))

  # 统一节点集合
  all_nodes <- union(rownames(adj1), rownames(adj2))

  # 初始化全零矩阵
  adj1_fixed <- matrix(0,
    nrow = length(all_nodes), ncol = length(all_nodes),
    dimnames = list(all_nodes, all_nodes)
  )
  adj2_fixed <- matrix(0,
    nrow = length(all_nodes), ncol = length(all_nodes),
    dimnames = list(all_nodes, all_nodes)
  )

  # 填充已知边
  adj1_fixed[rownames(adj1), colnames(adj1)] <- adj1
  adj2_fixed[rownames(adj2), colnames(adj2)] <- adj2

  # 计算相似性
  if (method == "frobenius") {
    diff_norm <- norm(adj1_fixed - adj2_fixed, "F")
    max_norm <- sqrt(nrow(adj1_fixed) * ncol(adj2_fixed))
    similarity <- 1 - diff_norm / max_norm
  } else if (method == "cosine") {
    similarity <- sum(adj1_fixed * adj2_fixed) /
      (norm(adj1_fixed, "F") * norm(adj2_fixed, "F"))
  } else {
    stop("Method must be 'frobenius' or 'cosine'.")
  }

  return(similarity)
}

#' Compare Two Networks
#'
#' @param g1 network1
#' @param g2 network2
#'
#' @returns A list containing the following elements:
#' - `g1`: The first network.
#' - `g2`: The second network.
#' - `g_union`: The union of the two networks.
#' - `g_inter`: The intersection of the two networks.
#' - `net_par_df`: A data frame containing the network parameters.
#' - `net_similarity`: A list containing the similarity metrics.
#'
#' @export
#'
#' @examples
#' data("c_net")
#' set.seed(12)
#' co_net_p1 <- c_net_filter(co_net, name %in% sample(V(co_net)$name, 300))
#' co_net_p2 <- c_net_filter(co_net, name %in% sample(V(co_net)$name, 300))
#' c_net_compare(co_net_p1, co_net_p2) -> c_net_comp
#' plot(c_net_comp)
c_net_compare <- function(g1, g2) {
  g_union <- suppressMessages(c_net_union(g1, g2))
  g_inter <- suppressMessages(c_net_intersect(g1, g2))

  lapply(list(g1, g2, g_union, g_inter), function(g) {
    net_par(g, mode = "n")$n_index
  }) %>%
    do.call(rbind, .) %>%
    pcutils::t2() -> net_par_df
  colnames(net_par_df) <- c("g1", "g2", "g_union", "g_inter")

  # 计算相似性
  net_similarity <- c(net_par_df[1:2, 4] / net_par_df[1:2, 3], adjacency_similarity(g1, g2))
  names(net_similarity) <- c("node_jaccard", "edge_jaccard", "adjacency_similarity")

  res <- list(
    g1 = g1, g2 = g2, g_union = g_union, g_inter = g_inter,
    net_par_df = net_par_df, net_similarity = net_similarity
  )
  class(res) <- c("metanet_compare")
  res
}

#' Plot a metanet_compare
#'
#' @param x metanet_compare object
#' @param ... add
#' @param coors_com coors object
#' @param mains a vector of two strings for the main titles of the two networks
#'
#' @return plot
#' @exportS3Method
#' @method plot metanet_compare
plot.metanet_compare <- function(x, coors_com = NULL, mains = NULL, ...) {
  c_net_comp <- x
  if (is.null(coors_com)) {
    coors_com <- c_net_layout(c_net_comp$g_union)
  }
  if (is.null(mains)) {
    mains <- c("g1", "g2")
  }

  c_net_highlight(c_net_comp$g1, V(c_net_comp$g_inter)$name) -> c_net_comp_g1
  c_net_highlight(c_net_comp$g2, V(c_net_comp$g_inter)$name) -> c_net_comp_g2
  plot_multi_nets(list(c_net_comp_g1, c_net_comp_g2),
    nrow = 1, ncol = 2,
    multi_params_list = list(
      list(coors = coors_com, main = mains[1], legend = FALSE),
      list(coors = coors_com, main = mains[2], legend = FALSE)
    )
  )
}
