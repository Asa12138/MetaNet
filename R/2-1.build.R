# =========2.build======

#' Construct a metanet from a corr object
#'
#' @param corr corr object from `c_net_calculate()` or `read_corr()`.
#' @param r_threshold r_threshold (default: >0.6).
#' @param p_threshold p_threshold (default: <0.05).
#' @param use_p_adj use the p.adjust instead of p.value (default: TRUE), if p.adjust not in the corr object, use p.value.
#' @param delete_single should delete single vertexes?
#'
#' @return an metanet object
#' @export
#'
#' @family build
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab) -> totu
#' metadata[, 3:10] -> env
#' c_net_calculate(totu) -> corr
#' c_net_build(corr, r_threshold = 0.65) -> co_net
#'
#' c_net_calculate(totu, env) -> corr2
#' c_net_build(corr2) -> co_net2
c_net_build <- function(corr, r_threshold = 0.6, p_threshold = 0.05, use_p_adj = TRUE, delete_single = TRUE) {
  stopifnot(inherits(corr, "corr"))
  # set thresholds to construct
  if ("r" %in% names(corr)) {
    occor.r <- corr$r
  } else {
    stop("No r in the input corr object.")
  }

  if (("p.adjust" %in% names(corr)) && use_p_adj) {
    occor.p <- corr$p.adjust
  } else {
    if ("p.value" %in% names(corr)) {
      occor.p <- corr$p.value
      message("Have not do p-value adjustment! use the p.value to build network.")
    } else {
      occor.p <- occor.r
      occor.p[occor.p != 0] <- 0
      message("No p.value or p.adjust in the input corr object, use the r to build network.")
    }
  }

  if (any(is.na(occor.r))) {
    occor.r[is.na(occor.r)] <- 0
    warning("NA values in the input corr object, set to 0.")
  }
  occor.r[occor.p > p_threshold | abs(occor.r) < r_threshold] <- 0

  # make igraph
  if (t_flag(occor.r)) {
    go <- igraph::graph_from_adjacency_matrix(as.matrix(occor.r),
      mode = "undirected", weighted = TRUE, diag = FALSE
    )
    igraph::graph.attributes(go)$n_type <- "single"
  } else {
    inset <- intersect(rownames(occor.r), colnames(occor.r))
    if (length(inset) > 0) {
      insets <- paste0(inset, collapse = ", ")
      warning("some nodes (row and column) has same name, please check: ", insets)
    }
    go <- igraph::graph_from_incidence_matrix(as.matrix(occor.r),
      directed = FALSE, weighted = TRUE
    )
    igraph::graph.attributes(go)$n_type <- "bipartite"
  }

  # delete single vertexes?
  if (delete_single) go <- igraph::delete.vertices(go, V(go)[igraph::degree(go) == 0])

  # set vertex attributes
  V(go)$v_group <- ifelse(V(go)$name %in% rownames(occor.r), "v_group1", "v_group2")
  V(go)$v_class <- ifelse(V(go)$name %in% rownames(occor.r), "v_class1", "v_class2")
  V(go)$size <- ceiling(60 / sqrt(length(V(go)))) + 1
  V(go)$label <- V(go)$name

  # abs edges weight
  go.weight <- E(go)$weight
  E(go)$cor <- go.weight
  E(go)$weight <- abs(go.weight)

  # time-consuming
  if (TRUE) {
    tmp_e <- get_e(go)
    if ("p.value" %in% names(corr)) {
      # E(go)$p.value <-get_e(go)%>%dplyr::select(from,to)%>%apply(., 1, \(x)occor.r$p.value[x[1],x[2]])
      tmp <- reshape2::melt(corr$p.value, varnames = c("from", "to"), value.name = "p.value", as.is = TRUE)
      tmp_e <- dplyr::left_join(tmp_e, tmp, by = c("from", "to"))
    }
    if ("p.adjust" %in% names(corr)) {
      # E(go)$p.adjust <-get_e(go)%>%dplyr::select(from,to)%>%apply(., 1, \(x)occor.r$p.adjust[x[1],x[2]])
      tmp <- reshape2::melt(corr$p.adjust, varnames = c("from", "to"), value.name = "p.adjust", as.is = TRUE)
      tmp_e <- dplyr::left_join(tmp_e, tmp, by = c("from", "to"))
    }
    # 应该直接expand，再left_join快很多
    igraph::edge.attributes(go) <- as.list(tmp_e)
  }

  # set edges type
  E(go)$e_type <- ifelse(go.weight > 0, "positive", "negative")
  # set edges width
  E(go)$width <- E(go)$weight

  if (!t_flag(occor.r)) {
    # set edges from_to
    anno_edge(go, get_v(go)[, c("name", "v_group")], verbose = FALSE) -> go
    # set edge intra-inter
    tmp_e <- igraph::edge.attributes(go)
    E(go)$e_class <- ifelse(tmp_e$v_group_from == tmp_e$v_group_to, "intra", "inter")
  }

  c_net_update(go, initialize = TRUE) -> go

  return(go)
}


#' Multi-omics network build
#'
#' @param ... some omics abundance tables
#' @param mode "full"
#' @param method "spearman" or "pearson"
#' @param filename the prefix of saved .corr file or FALSE
#' @param p.adjust.method see \code{\link[stats]{p.adjust}}
#' @param r_threshold r_threshold (default: >0.6)
#' @param p_threshold p_threshold (default: <0.05)
#' @param use_p_adj use the p.adjust instead of p-value (default: TRUE)
#' @param delete_single should delete single vertexes?
#'
#' @return metanet
#' @export
#' @family build
#' @examples
#' data("multi_test")
#' multi1 <- multi_net_build(list(Microbiome = micro, Metabolome = metab, Transcriptome = transc))
#' multi1 <- c_net_set(multi1, micro_g, metab_g, transc_g,
#'   vertex_class = c("Phylum", "kingdom", "type")
#' )
#' multi1 <- c_net_set(multi1, data.frame("Abundance1" = colSums(micro)),
#'   data.frame("Abundance2" = colSums(metab)), data.frame("Abundance3" = colSums(transc)),
#'   vertex_size = paste0("Abundance", 1:3)
#' )
#' c_net_plot(multi1)
multi_net_build <- function(..., mode = "full", method = "spearman",
                            filename = FALSE, p.adjust.method = NULL,
                            r_threshold = 0.6, p_threshold = 0.05, use_p_adj = TRUE, delete_single = TRUE) {
  tables <- list(...)
  if (all(class(tables[[1]]) == "list")) tables <- tables[[1]]
  tabs_name <- names(tables)

  tables <- check_tabs(tables)
  if (mode == "full") {
    all_totu <- do.call(cbind, tables)
    message("Calculating ", nrow(all_totu), " samples and ", ncol(all_totu), " features of ", length(tables), " groups.\n")
    all_corr <- c_net_calculate(all_totu, method = method, filename = filename, p.adjust.method = p.adjust.method)

    c_net_build(all_corr,
      r_threshold = r_threshold, p_threshold = p_threshold,
      use_p_adj = use_p_adj, delete_single = delete_single
    ) -> multi_net

    igraph::graph.attributes(multi_net)$n_type <- "multi_full"

    get_v(multi_net) -> tmp_v
    if (is.null(tabs_name)) tabs_name <- paste0("omic", seq_len(length(tables)))
    position <- rep(tabs_name, vapply(tables, ncol, numeric(1)))
    names(position) <- lapply(tables, colnames) %>% do.call(c, .)

    tmp_v$v_class <- tmp_v$v_group <- vapply(tmp_v$name, \(x)position[x], character(1))
    tmp_v %>% as.list() -> igraph::vertex_attr(multi_net)

    # set edges from_to
    suppressMessages(anno_edge(multi_net, tmp_v[, c("name", "v_group")], verbose = FALSE) -> multi_net)
    # set edge intra-inter
    tmp_e <- igraph::edge.attributes(multi_net)
    E(multi_net)$e_class <- ifelse(tmp_e$v_group_from == tmp_e$v_group_to, "intra", "inter")
    c_net_update(multi_net, initialize = TRUE, verbose = FALSE) -> multi_net
  }
  multi_net
}


default_v_color <- c(
  "#a6bce3", "#fb9a99", "#fdbf6f", "#1f78b4", "#b2df8a",
  "#cab2d6", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
  "#F8CC00", "#b15928"
)
default_v_color_num <- c(
  "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D",
  "#A50F15", "#67000D"
)
default_e_color <- c(
  "#a6cee3", "#78c679", "#c2a5cf", "#ff7f00", "#1f78b4",
  "#810f7c", "#F8CC00", "#006d2c", "#4d4d4d", "#8c510a",
  "#d73027", "#7f0000", "#41b6c4", "#e7298a", "#54278f"
)
default_e_color_num <- c(
  "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
  "#006D2C", "#00441B"
)

default_v_shape <- c("circle" = 21, "square" = 22, "diamond" = 23, "triangle1" = 24, "triangle2" = 25) # "triangle", "diamond", "hexagon", "star"

color_generate <- function(v_class, n_break = 5, mode = "v",
                           v_color = default_v_color, v_color_num = default_v_color_num,
                           e_color = default_e_color, e_color_num = default_e_color_num) {
  # ！！！现在的legend数字不好看，考虑pretty改进
  if (is.numeric(v_class)) {
    # n_color=n_break
    # n_label=seq(min(v_class, na.rm = TRUE), max(v_class, na.rm = TRUE), length.out = n_break)

    n_break <- pretty(v_class, n_break)
    n_color <- length(n_break) - 1
    n_label <- n_break[seq_len(n_color)] %>% format(., trim = TRUE)

    if (mode == "v") {
      return(as.character(cut(v_class,
        breaks = n_break,
        labels = pcutils::get_cols(n_color, v_color_num)
      )))
    } else if (mode == "e") {
      return(as.character(cut(v_class,
        breaks = n_break,
        labels = pcutils::get_cols(n_color, e_color_num)
      )))
    } else if (mode == "label") {
      return(as.character(cut(v_class,
        breaks = n_break,
        labels = n_label
      )))
    }
  } else {
    if (mode == "v") {
      return(pcutils::tidai(v_class, pcutils::get_cols(nlevels(factor(v_class)), v_color), fac = TRUE))
    } else if (mode == "e") {
      v_class <- droplevels(as.factor(v_class))
      if (all(levels(v_class) %in% c("negative", "positive"))) {
        ncols <- c(negative = "#E85D5D", positive = "#48A4F0")[levels(v_class)]
      } else if (all(levels(v_class) %in% c("inter-module", "intra-module"))) {
        ncols <- c("inter-module" = "#FA789A", "intra-module" = "#A6CEE3")[levels(v_class)]
      } else {
        ncols <- e_color
      }
      return(pcutils::tidai(v_class, pcutils::get_cols(nlevels(factor(v_class)), ncols), fac = TRUE))
    } else if (mode == "label") {
      return(v_class)
    }
  }
}

#' Update a metanet object or transform igraph object to metanet object
#'
#' @param go a metanet object or igraph object
#' @param node_break node_break if v_class is numeric, default: 5
#' @param initialize initialize?
#' @param edge_break edge_break if e_type is numeric, default: 5
#' @param verbose verbose?
#' @param uniq_v_class if TRUE, add prefix to v_class if multiple v_class belong to same v_group.
#'
#' @aliases as.metanet
#' @export
#' @family build
#' @return metanet
c_net_update <- function(go, node_break = 5, edge_break = 5, initialize = FALSE, verbose = TRUE, uniq_v_class = FALSE) {
  stopifnot(inherits(go, "igraph"))
  if (length(go) == 0) {
    message("Empty graph, return empty metanet object.")
    return(go)
  }
  v_class <- NULL
  # name
  if (!"name" %in% igraph::vertex_attr_names(go)) V(go)$name <- as.character(seq_len(length(go)))
  if (!"label" %in% igraph::vertex_attr_names(go)) V(go)$label <- V(go)$name
  get_v(go) -> tmp_v
  # v_size
  if (!"size" %in% colnames(tmp_v)) tmp_v$size <- 1

  # ！！！从这里开始仔细修改，考虑对应关系
  # v_shape
  if (!"v_group" %in% colnames(tmp_v)) {
    tmp_v$v_group <- "v_group1"
  } else {
    tmp_v$v_group <- as.character(tmp_v$v_group)
    if (length(unique(tmp_v$v_group)) > 10) warning("Too many groups, please check the 'v_group'.")
  }
  flag <- TRUE
  if ("shape" %in% colnames(tmp_v)) {
    if (initialize) {
      if (e_match(tmp_v[, c("v_group", "shape")])) {
        flag <- FALSE
      } else if (verbose) message("'v_group' and 'shape' not match one by one, update 'shape'.")
    } else {
      if (e_match(tmp_v[, c("v_group", "shape")], test = 1)) {
        flag <- FALSE
      } else if (verbose) message("a 'v_group' should not match multiple shapes, update 'shape'.")
    }
  }
  if (flag) tmp_v$shape <- pcutils::tidai(tmp_v$v_group, names(default_v_shape), fac = TRUE)

  # v_color
  if (!"v_class" %in% colnames(tmp_v)) {
    tmp_v$v_class <- "v_class1"
  }

  if (uniq_v_class) {
    # 如果有同样的v_class对应不同v_group，则加上前缀
    dplyr::distinct_all(tmp_v[, c("v_class", "v_group")]) %>% dplyr::count(v_class) -> group_class
    group_class <- setNames(group_class$n, group_class$v_class)
    if (any(group_class > 1)) {
      la <- which(group_class > 1)
      if (verbose) message("some 'v_class' belong to multiple 'v_group': \n", paste0(names(la), collapse = ", "))
      tmp_v$v_class %in% names(la) -> r_index
      tmp_v$v_class[r_index] <- paste0(tmp_v$v_group[r_index], "-", tmp_v$v_class[r_index])
    }
  }

  flag <- TRUE
  if ("color" %in% colnames(tmp_v)) {
    if (initialize) {
      if (e_match(tmp_v[, c("v_class", "color")])) {
        flag <- FALSE
      } else if (verbose) message("'v_class' and 'color' not match one by one, update 'color'.")
    } else {
      if (e_match(tmp_v[, c("v_class", "color")], test = 1)) {
        flag <- FALSE
      } else if (verbose) message("a 'v_class' should not match multiple colors, update 'color'.")
    }
  }

  if (flag) {
    tmp_v$color <- color_generate(tmp_v$v_class, node_break, mode = "v")
    tmp_v$v_class <- color_generate(tmp_v$v_class, node_break, mode = "label")
  }

  as.list(tmp_v) -> igraph::vertex_attr(go)

  # e_color
  if (!"e_type" %in% igraph::edge_attr_names(go)) {
    E(go)$e_type <- "e_type1"
  }
  flag <- TRUE
  if ("color" %in% igraph::edge_attr_names(go)) {
    if (initialize) {
      if (e_match(data.frame(a = E(go)$e_type, b = E(go)$color))) {
        flag <- FALSE
      } else if (verbose) message("'e_type' and 'color' not match one by one, update 'color'.")
    } else {
      if (e_match(data.frame(a = E(go)$e_type, b = E(go)$color), test = 1)) {
        flag <- FALSE
      } else if (verbose) message("a 'e_type' should not match multiple colors, update 'color'.")
    }
  }
  if (flag) {
    E(go)$color <- color_generate(E(go)$e_type, edge_break, mode = "e")
    E(go)$e_type <- color_generate(E(go)$e_type, edge_break, mode = "label")
  }

  # e_lty
  if (!"e_class" %in% igraph::edge_attr_names(go)) {
    E(go)$e_class <- "e_class1"
  }
  flag <- TRUE
  if ("lty" %in% igraph::edge_attr_names(go)) {
    if (initialize) {
      if (e_match(data.frame(a = E(go)$e_class, b = E(go)$lty))) {
        flag <- FALSE
      } else if (verbose) message("'e_class' and 'lty' not match one by one, update 'lty'.")
    } else {
      if (e_match(data.frame(a = E(go)$e_class, b = E(go)$lty), test = 1)) {
        flag <- FALSE
      } else if (verbose) message("a 'e_class' should not match multiple linetypes, update 'lty'.")
    }
  }
  if (flag) E(go)$lty <- pcutils::tidai(E(go)$e_class, 1:4, fac = TRUE)

  # e_width
  if (!"width" %in% igraph::edge_attr_names(go)) {
    E(go)$width <- 1
  }

  class(go) <- c("metanet", "igraph")
  return(go)
}


#' Clean a igraph object
#'
#' @param go igraph, metanet objects
#' @param direct direct?
#'
#' @return a igraph object
#' @export
clean_igraph <- function(go, direct = NULL) {
  stopifnot(inherits(go, "igraph"))
  if (is.null(direct)) direct <- is.directed(go)
  go <- igraph::graph_from_data_frame(
    d = get_e(go)[, c("from", "to")],
    directed = direct,
    vertices = get_v(go)["name"]
  )
  go
}


#' Construct a network from edge_list dataframe
#'
#' @param edgelist first is source, second is target, others are annotation
#' @param vertex_df vertex metadata data.frame
#' @param direct logical
#' @param e_type set e_type
#' @param e_class set e_class
#' @return metanet
#' @export
#' @family build
#' @examples
#' data(edgelist)
#' edge_net <- c_net_from_edgelist(arc_count, vertex_df = arc_taxonomy)
#' edge_net <- c_net_set(edge_net, vertex_class = "Phylum", edge_width = "n")
#' c_net_plot(edge_net)
c_net_from_edgelist <- function(edgelist, vertex_df = NULL, direct = FALSE, e_type = NULL, e_class = NULL) {
  node_type <- c("from", "to")
  if (!all(c("from", "to") %in% colnames(edgelist))) {
    message("No 'from' and 'to' in the colnames(edgelist), use the first two columns as the 'from' and 'to'.")
    node_type <- colnames(edgelist)[1:2]
    colnames(edgelist)[1:2] <- c("from", "to")
  }
  edgelist <- data.frame(edgelist[, c("from", "to"), drop = FALSE], edgelist[, setdiff(colnames(edgelist), c("from", "to")), drop = FALSE])
  name <- NULL
  vertices <- data.frame(
    name = c(unique(edgelist$from), unique(edgelist$to)),
    `_type` = c(rep(node_type[1], length(unique(edgelist$from))), rep(node_type[2], length(unique(edgelist$to)))),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  vertices <- dplyr::distinct(vertices, name, .keep_all = TRUE)
  if (!is.null(vertex_df)) {
    if (!"name" %in% colnames(vertex_df)) {
      vertex_df$name <- rownames(vertex_df)
      message("No 'name' in the colnames(vertex_df), use rownames(vertex_df) as the 'name'.")
    }
    vertex_df <- data.frame(vertex_df[, "name", drop = FALSE], vertex_df[, setdiff(colnames(vertex_df), "name"), drop = FALSE])
    vertices <- dplyr::left_join(vertices, vertex_df)
  }

  go <- igraph::graph_from_data_frame(edgelist, directed = direct, vertices = vertices)
  if (!is.null(e_type)) E(go)$e_type <- edgelist[, e_type]
  if (!is.null(e_class)) E(go)$e_class <- edgelist[, e_class]
  go <- c_net_update(go, initialize = TRUE)
  if (!"v_group" %in% colnames(vertices)) go <- c_net_set(go, vertex_group = "_type")
  if (!"v_class" %in% colnames(vertices)) go <- c_net_set(go, vertex_class = "_type")
  go
}

update_c_net_rda <- function() {
  otutab <- metadata <- taxonomy <- NULL
  data("otutab", package = "pcutils", envir = environment())
  t(otutab) -> totu
  metadata[, 3:10] -> env
  c_net_calculate(totu) -> corr
  c_net_build(corr, r_threshold = 0.65) -> co_net
  c_net_build(corr, r_threshold = 0.69) -> co_net_rmt

  c_net_calculate(totu, env) -> corr2
  c_net_build(corr2) -> co_net2

  co_net <- c_net_set(co_net, taxonomy, data.frame("Abundance" = colSums(totu)),
    vertex_class = "Phylum", vertex_size = "Abundance"
  )
  co_net2 <- c_net_set(co_net2, taxonomy, data.frame(name = colnames(env), env = colnames(env)),
    vertex_class = c("Phylum", "env")
  )
  co_net2 <- c_net_set(co_net2, data.frame("Abundance" = colSums(totu)), vertex_size = "Abundance")
  save(co_net, co_net2, co_net_rmt, file = "data/c_net.rda")
}
