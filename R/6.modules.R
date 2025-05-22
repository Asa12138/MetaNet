# ========6.modules=========

#' Generate a n-modules network
#'
#' @param module_number number of modules
#' @param n_node_in_module number of nodes in each modules
#' @param intra_module_density intra_module_density, recommend bigger than 20*inter_module_density, default:0.3
#' @param inter_module_density inter_module_density, default:0.01
#'
#' @description this is just a random generation method, the module number of result is not exactly the module_number, you can change the inter_module_density and intra_module_density to get the proper result.
#' @return n-modules metanet
#' @export
#'
#' @family module
#' @examples
#' g1 <- module_net()
#' get_n(g1)
#' plot(g1, mark_module = TRUE)
#' plot(g1, coors = g_layout(g1, zoom2 = 20))
#' plot(g1, coors = g_layout_polyarc(g1, group = "module"))
#' plot(g1, coors = g_layout_polygon(g1, group = "module"))
module_net <- function(module_number = 3, n_node_in_module = 30,
                       intra_module_density = 0.3,
                       inter_module_density = 0.01) {
  n_node_in_module <- rep(n_node_in_module, length = module_number)
  mat <- matrix(0, nrow = sum(n_node_in_module), ncol = sum(n_node_in_module))
  start <- c(1, cumsum(n_node_in_module[-length(n_node_in_module)]) + 1)
  end <- cumsum(n_node_in_module)
  seqls <- lapply(1:module_number, \(i)(start[i]:end[i]))
  # generate intra_module_edges
  for (i in seq_len(module_number)) {
    seq <- seqls[[i]]
    mat[seq, seq] <- as.matrix(igraph::get.adjacency(igraph::erdos.renyi.game(n_node_in_module[i], intra_module_density)))
  }
  # generate inter_module_edges
  idls <- combn(1:module_number, 2) %>% split(., col(.))
  for (ids in idls) {
    seq1 <- seqls[[ids[1]]]
    seq2 <- seqls[[ids[2]]]
    mat[seq1, seq2] <- sample(c(1, 0), length(seq1) * length(seq2),
      replace = TRUE,
      prob = c(inter_module_density, 1 - inter_module_density)
    )
  }

  g <- igraph::graph.adjacency(mat, mode = "undirected")
  # plot(g)
  c_net_update(g, initialize = TRUE, verbose = FALSE) -> g1
  g1 <- module_detect(g1)
  g1 <- to_module_net(g1)
  g1
}

#' Detect the modules
#'
#' @param go an igraph object
#' @param method cluster_method: "cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass"
#' @param n_node_in_module transfer the modules less than n_node_in_module to "others"
#' @param delete logical, delete others modules? default:FALSE, the others module will be "others".
#'
#' @return an igraph object
#' @export
#' @aliases c_net_module
#' @family module
#' @examples
#' data("c_net")
#' module_detect(co_net) -> co_net_modu
module_detect <- function(go, method = "cluster_fast_greedy", n_node_in_module = 0, delete = FALSE) {
  stopifnot(is_igraph(go))
  if (!is_metanet(go)) go <- c_net_update(go, initialize = TRUE, verbose = FALSE)

  if ("original_module" %in% vertex_attr_names(go)) message("'module' already exsited, start a new module detection!")
  ms <- c("cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass")
  method <- match.arg(method, ms)
  if ("weight" %in% edge_attr_names(go)) E(go)$weight <- abs(igraph::E(go)$weight)

  switch(method,
    "cluster_walktrap" = {
      wc <- igraph::cluster_walktrap(go, weights = NULL)
    },
    "cluster_edge_betweenness" = {
      wc <- igraph::cluster_edge_betweenness(go, weights = NULL)
    },
    "cluster_fast_greedy" = {
      wc <- igraph::cluster_fast_greedy(go, weights = NULL)
    },
    "cluster_spinglass" = {
      wc <- igraph::cluster_spinglass(go, weights = NULL)
    }
  )

  # add components
  V(go)$components <- igraph::components(go)$membership %>% as.character()
  V(go)$module <- igraph::membership(wc) %>% as.character()
  V(go)$original_module <- V(go)$module

  if (n_node_in_module > 0) {
    go <- filter_n_module(go, n_node_in_module = n_node_in_module, delete = delete)
    go <- c_net_update(go, initialize = TRUE, verbose = FALSE)
  }

  graph_attr(go)$communities <- wc
  graph_attr(go)$modularity <- modularity(wc)
  rand.g <- erdos.renyi.game(length(V(go)), length(E(go)), type = "gnm")
  rand_m <- modularity(cluster_fast_greedy(rand.g))
  relative_modularity <- (modularity(wc) - rand_m) / rand_m # relative modularity
  graph_attr(go)$relative_modularity <- relative_modularity

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
#' @param go_m metanet with modules
#' @param n_node_in_module transfer the modules less than n_node_in_module to "others"
#' @param keep_id keep modules ids, will not be "others"
#' @param delete logical, delete others modules? default:FALSE, the others module will be "others".
#' @family module
#' @return metanet with modules
#' @export
#' @examples
#' data("c_net")
#' module_detect(co_net) -> co_net_modu
#' filter_n_module(co_net_modu, n_node_in_module = 30) -> co_net_modu
#' if (requireNamespace("ggtree") && requireNamespace("treeio")) plot_module_tree(co_net_modu)
#' combine_n_module(co_net_modu, 20) -> co_net_modu1
#' if (requireNamespace("ggtree") && requireNamespace("treeio")) plot_module_tree(co_net_modu1)
filter_n_module <- function(go_m, n_node_in_module = 0, keep_id = NULL, delete = FALSE) {
  module <- NULL
  if (!"original_module" %in% vertex_attr_names(go_m)) stop("'module' do not exsited, please do a `module_detect` first!")
  members <- V(go_m)$original_module
  table(members) %>% sort(decreasing = TRUE) -> s_members

  # filter modules whose nodes bigger than n_node_in_module
  keep_id1 <- names(s_members[s_members >= n_node_in_module])
  keep_id <- base::union(as.character(keep_id), as.character(keep_id1))
  new_module <- ifelse(members %in% keep_id, members, "others")

  V(go_m)$module <- new_module
  if (delete) go_m <- c_net_filter(go_m, module != "others")
  return(go_m)
}

#' Combine or cut modules to module_number
#'
#' @param module_number number of modules
#' @export
#' @rdname filter_n_module
combine_n_module <- function(go_m, module_number = 5) {
  get_community(go_m) -> comm
  igraph::cut_at(comm, module_number) -> new_modu
  V(go_m)$module <- as.character(new_modu)
  graph.attributes(go_m)$communities$membership <- as.numeric(new_modu)
  go_m
}

#' Transformation a network to a module network
#'
#' @param edge_type "module", "module_from", "module_to"
#' @param go metanet
#'
#' @export
#' @family module
#' @return metanet with modules
to_module_net <- function(go, edge_type = c("module", "module_from", "module_to")[1]) {
  if (!"module" %in% vertex_attr_names(go)) stop("no 'module', please `module_detect()` first or set the V(net)$module.")
  edge_type <- match.arg(edge_type, c("module", "module_from", "module_to"))
  go <- anno_edge(go, get_v(go)[, c("name", "module")], verbose = FALSE)

  tmp_e <- igraph::edge.attributes(go)
  if (edge_type == "module") {
    E(go)$e_type <- ifelse(tmp_e$module_from == tmp_e$module_to, "intra-module", "inter-module")
  } else if (edge_type == "module_from") {
    E(go)$e_type <- tmp_e$module_from
  } else if (edge_type == "module_to") {
    E(go)$e_type <- tmp_e$module_to
  }

  # 刷新颜色
  # go=delete_edge_attr(go,"color")
  V(go)$v_class <- V(go)$module
  go <- c_net_update(go, initialize = TRUE, verbose = FALSE)
  V(go)$color <- ifelse(V(go)$module == "others", "grey", V(go)$color)

  n_mod <- unique(V(go)$module)
  igraph::graph.attributes(go)$n_type <- "module"
  igraph::graph.attributes(go)$n_modules <- length(n_mod[n_mod != "others"])
  go
}

#' Get community
#' @param go_m module metanet
#'
#' @export
#' @family module
#' @return community
get_community <- function(go_m) {
  if (is.null(igraph::graph_attr(go_m)$communities)) stop("No community find, please do module_net() first.")
  igraph::graph_attr(go_m)$communities
}

#' Get module
#' @param go_m module metanet
#'
#' @export
#' @family module
#' @return module
get_module <- function(go_m) {
  if (!"module" %in% vertex_attr_names(go_m)) stop("no modules, please `module_detect()` first")
  setNames(V(go_m)$module, V(go_m)$name)
}

#' Get module_eigen
#' @param go_m module metanet
#'
#' @export
#' @family module
#' @return module_eigen
get_module_eigen <- function(go_m) {
  graph_attr(go_m, "module_eigen")
}

#' Summary module index
#' @param go_m module metanet
#' @param var variable name
#' @param module which column name is module. default: "module"
#' @param ... add
#' @family module
#' @export
#' @return ggplot
#' @examples
#' data("c_net")
#' module_detect(co_net, n_node_in_module = 30) -> co_net_modu
#' summary_module(co_net_modu, var = "v_class", module = "module")
#' summary_module(co_net_modu, var = "Abundance", module = "module")
summary_module <- function(go_m, var = "v_class", module = "module", ...) {
  tmp_v <- get_v(go_m)
  if ((length(module) > 1) || (length(var) > 1)) stop("var or module should be one column!")
  a <- tmp_v %>% dplyr::select(!!module, !!var)
  colnames(a)[1] <- "module"

  i <- var
  if (is.numeric(a[, i])) {
    pcutils::group_box(a[i], group = "module", metadata = a, ...)
  } else {
    table(a[, i], a$module) %>%
      as.data.frame() %>%
      reshape2::acast(Var1 ~ Var2, value.var = "Freq") %>%
      as.data.frame() -> tab
    pcutils::stackplot(tab, legend_title = var, ...) + labs(x = "Module")
  }
}

#' Plot module tree
#' @param go_m module metanet
#' @param module which column name is module. default: "module"
#' @param community community object, default: NULL, use the community of go_m
#' @param label.size label.size
#' @return ggplot
#'
#' @export
#' @rdname filter_n_module
plot_module_tree <- function(go_m, module = "module", community = NULL, label.size = 2) {
  tmp_v <- get_v(go_m)
  mdata <- tmp_v[, c("name", module)]

  lib_ps("ggtree", "treeio", library = FALSE)
  # modules tree
  if (is.null(community)) {
    get_community(go_m) %>% treeio::as.phylo() -> mcl
  } else {
    community %>% treeio::as.phylo() -> mcl
  }
  mcl <- dplyr::left_join(mcl, mdata, by = c("label" = "name"))
  p <- ggtree::ggtree(mcl, size = 0.3) +
    ggtree::geom_tiplab(aes(color = module), show.legend = FALSE, size = label.size) +
    scale_color_manual(values = pcutils::get_cols(length(unique(mdata$module)), "col3"))

  ggtree::gheatmap(p, mdata %>% tibble::column_to_rownames("name")) +
    scale_fill_manual(values = pcutils::get_cols(length(unique(mdata$module)), "col3"), name = NULL)
}

#' Calculate the eigenvalue of each module and correlation of nodes and eigenvalue (node_eigen_cor).
#'
#' @param go_m module metanet
#' @param totu original abundance table
#' @param cor_method "pearson", "kendall", "spearman"
#'
#' @export
#' @return module metanet with module_eigen
#' @rdname module_expression
module_eigen <- function(go_m, totu, cor_method = "spearman") {
  modules <- get_module(go_m)
  totu <- totu[, names(modules)]

  res <- lapply(levels(factor(modules)), \(i){
    if (i == "others") {
      return(NULL)
    }
    totu1 <- totu[, modules == i]
    # PCA
    # pc <- prcomp(totu1)
    # pc$x[, 1]
    # cor(pc$x[, 1],rowMeans(totu1))

    # SVD
    totu_scale <- scale(totu1)
    svd <- svd(totu_scale)
    if (cor(svd$u[, 1], rowMeans(totu_scale)) < 0) {
      return(-svd$u[, 1])
    } else {
      return(svd$u[, 1])
    }
  })
  names(res) <- levels(factor(modules))
  # simplify method
  eigen_res <- do.call(cbind, res) %>% data.frame(., check.names = FALSE)
  # colnames(eigen_res)=paste0("ME_",colnames(eigen_res))
  rownames(eigen_res) <- rownames(totu)

  lapply(colnames(eigen_res), \(i)cor(totu[, modules == i], eigen_res[, i], method = cor_method)) %>%
    do.call(rbind, .) %>%
    as.data.frame() -> node_eigen_cor
  colnames(node_eigen_cor) <- "node_eigen_cor"
  go_m <- anno_vertex(go_m, node_eigen_cor, verbose = FALSE)
  graph.attributes(go_m)$module_eigen <- eigen_res
  go_m
}

#' Plot the expression of each modules
#'
#' @param go_m module metanet
#' @param totu original abundance table used for module_eigen().
#' @param group group variable for totu
#' @param r_threshold the threshold for node_eigen_cor, default: 0.6.
#' @param x_order order the x axis.
#' @param facet_param parameters parse to \code{\link[ggplot2]{facet_wrap}}, e.g. nrow=2.
#' @param plot_eigen plot the eigen value line.
#' @family module
#' @export
#'
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab) -> totu
#' data("c_net")
#' module_detect(co_net, n_node_in_module = 30) -> co_net_modu
#' module_eigen(co_net_modu, totu) -> co_net_modu
#' module_expression(co_net_modu, totu)
module_expression <- function(go_m, totu, group = NULL, r_threshold = 0.6,
                              x_order = NULL, facet_param = NULL, plot_eigen = FALSE) {
  node_eigen_cor <- variable <- value <- name <- module <- rowname <- NULL
  if (is.null(graph_attr(go_m, "module_eigen"))) stop("Please do module_eigen() first")

  graph_attr(go_m, "module_eigen") %>%
    rownames_to_column() %>%
    reshape2::melt(id.vars = "rowname", variable.name = "module") -> module_eigen

  get_v(go_m) -> tmp_v

  if (!is.null(group)) totu <- pcutils::hebing(totu, group, 1)

  totu_scale <- scale(totu)
  pdat <- cbind(tmp_v[, c("name", "module", "node_eigen_cor")], t(totu_scale[, tmp_v$name]))
  pdat <- dplyr::filter(pdat, node_eigen_cor > r_threshold)
  pdat_m <- reshape2::melt(pdat, id.vars = c("name", "module", "node_eigen_cor"))
  if (!is.null(x_order)) pdat_m$variable <- pcutils::change_fac_lev(pdat_m$variable, level = x_order)
  pdat_m$module <- factor(pdat_m$module)

  p1 <- ggplot() +
    geom_line(data = pdat_m, aes(
      x = variable, y = value, group = name, color = module,
      size = node_eigen_cor, alpha = node_eigen_cor
    ), size = 0.8) +
    MetaNet_theme +
    theme(plot.margin = unit(c(1, 2, 1, 1), "lines")) +
    do.call(facet_wrap, append(list(facets = ~module), pcutils::update_param(list(scales = "free_y", ncol = 2), facet_param))) +
    labs(x = NULL, y = NULL) +
    # facet_wrap(facets = ~module,nrow = nrow,scales = "free_y")+
    scale_alpha_continuous(range = c(0, 0.6)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_color_manual(values = get_cols(nlevels(pdat_m$module), default_v_color))

  if (plot_eigen) p1 <- p1 + geom_line(data = module_eigen, aes(x = rowname, y = value, col = module, group = module), size = 2, alpha = 1)

  return(p1)
}

#' Zi-Pi calculate
#'
#' @param go_m igraph object after `module_detect()`
#' @param mode use 7-group (mode=1) or 4-group (mode=2), default: mode=2
#' @param use_origin use original_module, default:TRUE, if FALSE, use module
#'
#' @return igraph
#' @export
#' @family module
#' @references 1. Guimerà, R. & Amaral, L. Functional cartography of complex metabolic networks. (2005) doi:10.1038/nature03288.
#' @examples
#' data("c_net")
#' module_detect(co_net) -> co_net_modu
#' zp_analyse(co_net_modu) -> co_net_modu
#' zp_plot(co_net_modu)
#' zp_plot(co_net_modu, mode = 3)
#'
zp_analyse <- function(go_m, mode = 2, use_origin = TRUE) {
  go_m -> go1
  v_index <- get_v(go_m)
  if (!"module" %in% names(v_index)) stop("no modules, please `module_detect()` first")
  if ("roles" %in% names(v_index)) message("areadly has roles, overwrite!")

  if (!"original_module" %in% names(v_index)) {
    if (use_origin) {
      use_origin <- FALSE
    }
  }
  # use original_module to do zp_analyse
  if (use_origin) {
    {
      go1 <- anno_vertex(go1, data.frame(row.names = V(go1)$name, module = V(go1)$original_module %>% as.numeric()), verbose = FALSE)
    } %>% suppressMessages()
  } else {
    if ("others" %in% v_index$module) message("Consider others as one module!")
    {
      go1 <- anno_vertex(go1, data.frame(
        row.names = V(go1)$name,
        module = tidai(v_index$module, seq_along(unique(v_index$module)))
      ), verbose = FALSE)
    } %>% suppressMessages()
  }

  within <- within_module_deg_z_score(go1)
  v_index$Ki <- within$Ki
  v_index$Zi <- within$Zi
  pc <- part_coeff(go1)
  v_index$Pi <- pc$Pi

  if (mode == 1) {
    lab <- c("Ultra-peripherals", "Peripherals", "Non-hub connectors", "Non-hub kinless nodes", "Provincial hubs", "Connector hubs", "Kinless hubs")
    backs <- data.frame(
      x1 = c(0, 0.05, 0.62, 0.8, 0, 0.3, 0.75),
      x2 = c(0.05, 0.62, 0.8, 1, 0.3, 0.75, 1),
      y1 = c(-Inf, -Inf, -Inf, -Inf, 2.5, 2.5, 2.5),
      y2 = c(2.5, 2.5, 2.5, 2.5, Inf, Inf, Inf),
      lab = factor(lab, levels = lab)
    )
  } else if (mode == 2) {
    lab <- c("Peripherals", "Network hubs", "Module hubs", "Connectors")
    backs <- data.frame(
      x1 = c(0, 0.62, 0, 0.62),
      x2 = c(0.62, 1, 0.62, 1),
      y1 = c(-Inf, 2.5, 2.5, -Inf),
      y2 = c(2.5, Inf, Inf, 2.5),
      lab = factor(lab, levels = lab)
    )
  }
  deter_role <- \(x, y, backs = backs){
    for (i in seq_len(nrow(backs))) {
      flag <- dplyr::between(as.numeric(x), backs$x1[i], backs$x2[i]) && dplyr::between(as.numeric(y), backs$y1[i], backs$y2[i])
      if (is.na(flag)) {
        return(NA)
      } else if (flag) {
        # if((backs$x1[i]<=x)&(backs$x2[i]>=x)&(backs$y1[i]<=y)&(backs$y2[i]>=y)){
        role <- backs$lab[i]
        break
      }
    }
    return(role)
  }
  v_index$roles <- apply(v_index, 1, \(x)deter_role(x["Pi"], x["Zi"], backs))

  vertex.attributes(go_m) <- as.list(v_index)
  return(go_m)
}


#' calculate Zi
#'
#' @param g igraph object
#' @param A adjacency matrix
#' @param weighted logical, default: FALSE
#'
#' @return within_module_deg_z_score
#' @noRd
#' @references https://github.com/cwatson/brainGraph/blob/master/R/vertex_roles.R
within_module_deg_z_score <- function(g, A = NULL, weighted = FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse = FALSE, names = TRUE, attr = "weight")
    } else {
      A <- as_adj(g, sparse = FALSE, names = TRUE)
    }
  }
  memb <- vertex_attr(g, "module") %>% as.numeric()
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
  memb <- vertex_attr(g, "module") %>% as.numeric()
  Ki <- colSums(A)
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
#' @param mode plot style, 1~3
#'
#' @return a ggplot object
#' @export
#' @rdname zp_analyse
zp_plot <- function(go, label = TRUE, mode = 1) {
  v_class <- value <- size <- roles <- x1 <- x2 <- y1 <- y2 <- Pi <- Zi <- name <- NULL
  lib_ps("ggrepel", library = FALSE)
  get_v(go) -> taxa.roles
  if (!"roles" %in% names(taxa.roles)) stop("no roles, please zp_analyse() first")

  if (mode == 3) {
    reshape2::melt(taxa.roles, measure.vars = c("Zi", "Pi")) -> taxa.roles1
    p <- ggplot(taxa.roles1, aes(x = v_class, y = value, col = v_class, size = size, shape = roles)) +
      geom_point() +
      facet_grid(variable ~ ., scales = "free_y") +
      theme_bw() +
      labs(x = NULL, y = NULL) +
      guides(col = "none") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      scale_color_manual(values = setNames(unique(taxa.roles$color), unique(taxa.roles$v_class)))
    return(p)
  }

  mode <- ifelse(nlevels(V(go)$roles) == 7, 1, 2)
  if (mode == 1) {
    lab <- c("Ultra-peripherals", "Peripherals", "Non-hub connectors", "Non-hub kinless nodes", "Provincial hubs", "Connector hubs", "Kinless hubs")
    CPCOLS <- c("#FCF6EFFC", "#EEBCF5", "#EDEDA4", "#FAA371", "#FC5D6096", "#9BC799B9", "#94CCF2AC")
    names(CPCOLS) <- lab
    backs <- data.frame(
      x1 = c(0, 0.05, 0.62, 0.8, 0, 0.3, 0.75),
      x2 = c(0.05, 0.62, 0.8, 1, 0.3, 0.75, 1),
      y1 = c(-Inf, -Inf, -Inf, -Inf, 2.5, 2.5, 2.5),
      y2 = c(2.5, 2.5, 2.5, 2.5, Inf, Inf, Inf),
      lab = factor(lab, levels = lab)
    )
  } else if (mode == 2) {
    lab <- c("Peripherals", "Network hubs", "Module hubs", "Connectors")
    CPCOLS <- c("#FCF6EFFC", "#FC5D6096", "#9BC799B9", "#94CCF2AC")
    names(CPCOLS) <- lab
    backs <- data.frame(
      x1 = c(0, 0.62, 0, 0.62),
      x2 = c(0.62, 1, 0.62, 1),
      y1 = c(-Inf, 2.5, 2.5, -Inf),
      y2 = c(2.5, Inf, Inf, 2.5),
      lab = factor(lab, levels = lab)
    )
  }
  p <- ggplot() +
    geom_rect(data = backs, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = lab), alpha = 0.7) +
    guides(fill = guide_legend(title = "Topological roles")) +
    scale_fill_manual(values = CPCOLS) +
    geom_point(data = taxa.roles, aes(x = Pi, y = Zi, color = factor(v_class))) +
    scale_color_manual(values = setNames(unique(taxa.roles$color), unique(taxa.roles$v_class))) +
    MetaNet_theme +
    guides(colour = "none") +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Participation Coefficient (Pi)") +
    ylab("Within-module connectivity (Zi)")

  if (label) {
    label_dat <- taxa.roles[!taxa.roles$roles %in% c("Peripherals", "Ultra-peripherals"), ]
    label_dat <- label_dat[!is.na(label_dat$roles), ]
    p <- p + ggrepel::geom_text_repel(
      data = label_dat,
      aes(x = Pi, y = Zi, label = name), size = 3
    )
  }
  return(p)
}
