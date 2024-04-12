# ==========2.1 manipulate========

#' Set basic attributes from totu table
#'
#' @param go metanet an igraph object
#' @param ... some data.frames to annotate go
#' @param vertex_group choose which column to be vertex_group (map to vertex_shape)
#' @param vertex_class choose which column to be vertex_class (map to vertex_color)
#' @param vertex_size choose which column to be vertex_size (map to vertex_size)
#' @param edge_type choose which column to be edge_type (map to edge_color)
#' @param edge_class choose which column to be edge_class (map to edge_linetype)
#' @param edge_width choose which column to be edge_width (map to edge_width)
#' @param node_break node_break if v_class is numeric, default: 5
#' @param edge_break edge_break if e_type is numeric, default: 5
#'
#' @return a metanet object
#' @export
#' @family build
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab) -> totu
#' metadata[, 3:10] -> env
#'
#' data("c_net")
#' co_net <- c_net_set(co_net, taxonomy, data.frame("Abundance" = colSums(totu)),
#'   vertex_class = "Phylum", vertex_size = "Abundance"
#' )
#' co_net2 <- c_net_set(co_net2, taxonomy, data.frame(name = colnames(env), env = colnames(env)),
#'   vertex_class = c("Phylum", "env")
#' )
#' co_net2 <- c_net_set(co_net2, data.frame("Abundance" = colSums(totu)), vertex_size = "Abundance")
c_net_set <- function(go, ..., vertex_group = "v_group", vertex_class = "v_class", vertex_size = "size",
                      edge_type = "e_type", edge_class = "e_class", edge_width = "width", node_break = 5, edge_break = 5) {
  size <- e_class <- width <- NULL
  c_net_update(go, verbose = FALSE) -> go
  name <- v_group <- v_class <- e_type <- color <- NULL

  # annotation vertex
  anno_dfs <- list(...)
  if (length(anno_dfs) > 0) {
    anno_dfs2 <- list()
    for (i in seq_len(length(anno_dfs))) {
      x <- anno_dfs[[i]]
      if ("name" %in% colnames(x)) {
        rownames(x) <- x$name
        x <- dplyr::select(x, -name)
      }
      anno_dfs2[[i]] <- x
    }

    if (any(duplicated(lapply(anno_dfs2, names) %>% unlist()))) stop("Duplicated column names in your annotation tables, please check!")

    Reduce(\(x, y)merge(x, y, by = "row.names", all = TRUE) %>%
      tibble::column_to_rownames("Row.names"), anno_dfs2) -> all_anno

    anno_vertex(go, all_anno) -> go
  }
  get_v(go) -> v_index
  get_e(go) -> e_index

  # set something
  # ！！！这里的set要改成跟c_net_update一样的逻辑
  if (!setequal(vertex_group, "v_group")) dplyr::select(v_index, v_group, !!vertex_group) %>% condance() -> v_index$v_group

  if (!setequal(vertex_class, "v_class")) {
    old_color <- twocol2vector(v_index[, c("v_class", "color")])
    new_color_name <- c()

    # 给每一个v_group加上v_class调整颜色
    # 可能某一个group用numeric做v_class，所以要分开上色
    for (i in unique(v_index$v_group)) {
      tmp_index <- v_index[v_index$v_group == i, ]
      tmp_v_class <- dplyr::select(tmp_index, v_class, !!vertex_class) %>% condance()
      if (identical(tmp_v_class, tmp_index$v_class)) {
        new_color_name <- c(new_color_name, unique(tmp_index$v_class))
        next
      }
      if (is.numeric(tmp_v_class)) {
        tmp_v_color <- color_generate(tmp_v_class, n_break = node_break, mode = "v")
        tmp_v_class <- color_generate(tmp_v_class, n_break = node_break, mode = "label")
        v_index[v_index$v_group == i, "v_class"] <- tmp_v_class
        v_index[v_index$v_group == i, "color"] <- tmp_v_color
      } else {
        new_color_name <- c(new_color_name, unique(tmp_index$v_class))
      }
    }
    # 总体分类颜色是否改变，没变的话就不该，变了的话全部重新赋
    new_color_name <- unique(new_color_name)
    if (!all(new_color_name %in% names(old_color))) {
      new_color <- setNames(pcutils::get_cols(length(new_color_name), pal = default_v_color), new_color_name)
      v_index$color <- condance(data.frame(
        v_index$color,
        pcutils::tidai(v_index$v_class, new_color)
      ))
    }
  }

  if (!setequal(vertex_size, "size")) dplyr::select(v_index, size, !!vertex_size) %>% condance() -> v_index$size

  if (!setequal(edge_type, "e_type")) {
    tmp_e_type <- dplyr::select(e_index, e_type, !!edge_type) %>% condance()
    if (!identical(tmp_e_type, e_index$e_type)) {
      tmp_e_color <- color_generate(tmp_e_type, edge_break, mode = "e")
      tmp_e_type <- color_generate(tmp_e_type, edge_break, mode = "label")
      e_index$e_type <- tmp_e_type
      e_index$color <- tmp_e_color
    }
  }
  if (!setequal(edge_class, "e_class")) dplyr::select(e_index, e_class, !!edge_class) %>% condance() -> e_index$e_class
  if (!setequal(edge_width, "width")) dplyr::select(e_index, width, !!edge_width) %>% condance() -> e_index$width

  as.list(v_index) -> igraph::vertex.attributes(go)
  as.list(e_index) -> igraph::edge.attributes(go)

  c_net_update(go, verbose = FALSE) -> go2
  return(go2)
}


#' Is this object a metanet object?
#'
#' @param go a test object
#'
#' @return logical
#' @export
#' @aliases is.metanet
#' @family manipulate
#' @examples
#' data(c_net)
#' is_metanet(co_net)
is_metanet <- function(go) {
  is.igraph(go) & inherits(go, "metanet")
}

#' Get vertex information
#'
#' @param go metanet object
#' @param name attribute name, default: NULL
#' @family manipulate
#' @return data.frame
#' @export
get_v <- function(go, name = NULL) {
  # 规定name只能为字符
  if (is.null(V(go)$name)) V(go)$name <- as.character(V(go))
  # df <- as.data.frame(igraph::vertex.attributes(go))
  igraph::as_data_frame(go, what = "v") -> df
  rownames(df) <- NULL
  if (!is.null(name)) {
    return(dplyr::select(df, !!name))
  } else {
    return(df)
  }
}

#' Get edge information
#' @param go metanet object
#' @param name attribute name, default: NULL
#' @return data.frame
#' @family manipulate
#' @export
get_e <- function(go, name = NULL) {
  id <- NULL
  tmp_e <- cbind_new(igraph::as_data_frame(go), data.frame(id = seq_len(igraph::ecount(go))))
  tmp_e <- dplyr::select(tmp_e, id, dplyr::everything())
  if (!is.null(name)) {
    return(dplyr::select(tmp_e, !!name))
  } else {
    return(tmp_e)
  }
}

#' Get network information
#'
#' @param go metanet object
#' @param name attribute name, default: NULL
#' @param simple logical, get simple index
#' @family manipulate
#' @return data.frame
#' @export
get_n <- function(go, name = NULL, simple = FALSE) {
  gls <- igraph::graph.attributes(go)
  if (simple) {
    gls <- lapply(gls, \(x){
      if (inherits(x, "data.frame")) {
        return(NULL)
      }
      if (is.array(x)) {
        return(NULL)
      }
      if (is.list(x)) {
        return(NULL)
      }
      if (length(x) > 1) {
        return(NULL)
      }
      return(x)
    })
  } else {
    gls <- lapply(gls, \(x){
      if (inherits(x, "data.frame")) {
        return(paste0(ncol(x), "-columns df"))
      }
      if (is.array(x)) {
        return(paste0(length(x), "-elements ", class(x)))
      }
      if (is.list(x)) {
        return(paste0(length(x), "-elements ", class(x)))
      }
      if (length(x) > 1) {
        return(paste0(length(x), "-elements vector"))
      }
      return(x)
    })
  }
  df <- as.data.frame(do.call(cbind, gls))
  if (!is.null(name)) {
    return(dplyr::select(df, !!name))
  } else {
    return(df)
  }
}

#' Filter a network according to some attributes
#'
#' @param go metanet object
#' @param ... some attributes of vertex and edge
#' @param mode "v" or "e"
#'
#' @return metanet
#' @export
#' @family manipulate
#' @examples
#' data("multi_net")
#' c_net_filter(multi1, v_group %in% c("omic1", "omic2"))
c_net_filter <- function(go, ..., mode = "v") {
  if (mode == "v") {
    go1 <- filter_v(go, ...)
  } else if (mode == "e") {
    go1 <- filter_e(go, ...)
  } else {
    stop("mode should be 'v' or 'e'")
  }
  if (length(V(go1)) == 0) {
    message("The network is empty.")
  }
  go1
}

filter_v <- function(go, ...) {
  get_v(go) -> tmp_v
  tmp_v <- dplyr::filter(tmp_v, ...)
  tmp_v$name -> vid
  igraph::subgraph(go, vid) -> go1
  class(go1) <- c("metanet", "igraph")
  go1
}

filter_e <- function(go, ...) {
  get_e(go) -> tmp_e
  tmp_e <- dplyr::filter(tmp_e, ...)
  tmp_e$id -> eid
  igraph::subgraph.edges(go, eid) -> go1
  class(go1) <- c("metanet", "igraph")
  go1
}


#' Union two networks
#'
#' @param go1 metanet object
#' @param go2 metanet object
#'
#' @return metanet
#' @export
#' @family manipulate
#' @examples
#' data("c_net")
#' co_net_union <- c_net_union(co_net, co_net2)
#' c_net_plot(co_net_union)
c_net_union <- function(go1, go2) {
  tmp_v1 <- get_v(go1)
  tmp_v2 <- get_v(go2)
  cols <- c("name", "label", "size", "v_group", "shape", "v_class", "color")
  tmp_v <- rbind(tmp_v1[cols], tmp_v2[cols])
  message("Duplicated vertexes: ", sum(duplicated(tmp_v$name)), "\nUse the attributes of the first network.")
  tmp_v <- tmp_v[!duplicated(tmp_v$name), ]

  tmp_e1 <- get_e(go1)
  tmp_e2 <- get_e(go2)
  cols <- c("from", "to", "e_type", "color", "e_class", "lty", "width")
  tmp_e <- rbind(tmp_e1[cols], tmp_e2[cols])
  message("Duplicated edges: ", sum(duplicated(tmp_e[, c("from", "to")])), "\nUse the attributes of the first network.")
  tmp_e <- tmp_e[!duplicated(tmp_e[, c("from", "to")]), ]

  go <- igraph::union(go1, go2)
  go <- clean_igraph(go, direct = FALSE)
  go <- c_net_annotate(go, tmp_v, mode = "v")
  go <- c_net_annotate(go, tmp_e, mode = "e")
  go <- c_net_annotate(go, list(n_type = "combine_net"), mode = "n")
  go <- c_net_update(go, initialize = TRUE)
  go
}


#' Annotate a metanet
#'
#' @param go metanet object
#' @param anno_tab a dataframe using to annotate (mode v, e), or a list (mode n)
#' @param mode "v" for vertex, "e" for edge, "n" for network
#' @param verbose logical
#'
#' @return a annotated metanet object
#' @export
#' @family manipulate
#' @examples
#' data("c_net")
#' anno <- data.frame("name" = "s__Pelomonas_puraquae", new_atr = "new")
#' co_net_new <- c_net_annotate(co_net, anno, mode = "v")
#' get_v(co_net_new, c("name", "new_atr"))
#'
#' anno <- data.frame("from" = "s__Pelomonas_puraquae", "to" = "s__un_g__Rhizobium", new_atr = "new")
#' co_net_new <- c_net_annotate(co_net, anno, mode = "e")
#' get_e(co_net_new, c("from", "to", "new_atr"))
#'
#' co_net_new <- c_net_annotate(co_net, list(new_atr = "new"), mode = "n")
#' get_n(co_net_new)
c_net_annotate <- function(go, anno_tab, mode = "v", verbose = TRUE) {
  mode <- match.arg(mode, c("v", "e", "n"))
  if (mode == "v") {
    anno_vertex(go, anno_tab, verbose = verbose) -> go
  } else if (mode == "e") {
    anno_edge(go, anno_tab, verbose = verbose) -> go
  } else if (mode == "n") {
    igraph::graph.attributes(go) <-
      pcutils::update_param(igraph::graph.attributes(go), anno_tab)
  }
  go
}


#' Use data.frame to annotate vertexes of metanet
#'
#' @param go metanet object
#' @param verbose logical
#' @param anno_tab a dataframe using to annotate (with rowname or a "name" column)
#'
#' @return a annotated metanet object
#' @aliases anno_node
#' @export
#' @family manipulate
#' @examples
#' data("c_net")
#' data("otutab", package = "pcutils")
#' anno_vertex(co_net, taxonomy)
anno_vertex <- function(go, anno_tab, verbose = TRUE) {
  if (is.null(anno_tab)) {
    return(go)
  }
  get_v(go) -> v_atr
  if (!"name" %in% colnames(anno_tab)) rownames(anno_tab) -> anno_tab$name
  if (any(duplicated(anno_tab$name))) {
    stop(
      "Duplicated name in annotation tables: ",
      paste0(anno_tab$name[duplicated(anno_tab$name)], collapse = ", ")
    )
  }
  v_atr <- dplyr::left_join(v_atr, anno_tab, by = "name", suffix = c(".x", ""))
  grep(".x", colnames(v_atr), value = TRUE) %>% gsub(".x", "", .) -> du
  if (length(du) > 0) message(length(du), (" attributes will be overwrited:\n"), paste0(du, collapse = ", "), "\n")
  v_atr %>% dplyr::select(!dplyr::ends_with(".x")) -> v_atr

  as.list(v_atr) -> igraph::vertex.attributes(go)
  return(go)
}

#' Use dataframe to annotate edges of an igraph
#'
#' @param go metanet an igraph object
#' @param verbose logical
#' @param anno_tab a dataframe using to annotate (with rowname or a name column)
#'
#' @return a annotated igraph object
#' @export
#' @family manipulate
#' @examples
#' data("c_net")
#' anno <- data.frame("from" = "s__Pelomonas_puraquae", "to" = "s__un_g__Rhizobium", new_atr = "new")
#' anno_edge(co_net, anno) -> anno_net
anno_edge <- function(go, anno_tab, verbose = TRUE) {
  name <- NULL
  if (is.null(anno_tab)) {
    return(go)
  }
  get_e(go) -> e_atr
  if (all(c("from", "to") %in% colnames(anno_tab))) {
    e_atr <- dplyr::left_join(e_atr, anno_tab, by = c("from", "to"), suffix = c(".x", ""))
    grep(".x", colnames(e_atr), value = TRUE) %>% gsub(".x", "", .) -> du
    if (length(du) > 0) {
      if (verbose) message(length(du), (" attributes will be overwrited:\n"), paste0(du, collapse = ","), "\n")
    }
    e_atr %>% dplyr::select(!dplyr::ends_with(".x")) -> e_atr
  } else {
    if (verbose) message("No 'from' and 'to' columns in annotation table, will use 'name_from' and 'name_to' instead.")
    if (!"name" %in% colnames(anno_tab)) rownames(anno_tab) -> anno_tab$name
    anno_tab %>% dplyr::select(name, dplyr::everything()) -> anno_tab
    # from
    tmp <- anno_tab
    colnames(tmp) <- paste0(colnames(anno_tab), "_from")
    e_atr <- dplyr::left_join(e_atr, tmp, by = c("from" = "name_from"), suffix = c(".x", ""))
    grep(".x", colnames(e_atr), value = TRUE) %>% gsub(".x", "", .) -> du
    if (length(du) > 0) {
      if (verbose) message(length(du), (" attributes will be overwrited:\n"), paste0(du, collapse = ","), "\n")
    }
    e_atr %>% dplyr::select(!dplyr::ends_with(".x")) -> e_atr
    # to
    tmp <- anno_tab
    colnames(tmp) <- paste0(colnames(anno_tab), "_to")
    e_atr <- dplyr::left_join(e_atr, tmp, by = c("to" = "name_to"), suffix = c(".x", ""))
    grep(".x", colnames(e_atr), value = TRUE) %>% gsub(".x", "", .) -> du
    if (length(du) > 0) {
      if (verbose) message(length(du), (" attributes will be overwrited:\n"), paste0(du, collapse = ","), "\n")
    }
    e_atr %>% dplyr::select(!dplyr::ends_with(".x")) -> e_atr
  }
  as.list(e_atr) -> igraph::edge.attributes(go)
  return(go)
}

#' Save network file
#'
#' @param go metanet network
#' @param filename filename
#' @param format "data.frame","graphml"
#' @return No value
#' @family manipulate
#' @export
c_net_save <- function(go, filename = "net", format = "data.frame") {
  if (format == "data.frame") {
    get_v(go) %>% write.csv(., paste0(filename, "_nodes.csv"), row.names = FALSE)
    get_e(go) %>%
      dplyr::select(-1) %>%
      write.csv(., paste0(filename, "_edges.csv"), row.names = FALSE)
  } else if (format == "graphml") {
    go <- igraph::delete_edge_attr(go, "id")
    if (!grepl("\\.graphml$", filename)) filename <- paste0(filename, ".graphml")
    igraph::write_graph(go, filename, format = "graphml")
  } else {
    if (!grepl(paste0("\\.", format), filename)) filename <- paste0(filename, ".", format)
    igraph::write_graph(go, filename, format = format)
  }
  message(paste0(filename, " saved sucessfully!"))
}

#' Load network file
#'
#' @inheritParams c_net_save
#'
#' @return metanet
#' @export
#' @family manipulate
c_net_load <- function(filename, format = "data.frame") {
  if (format == "data.frame") {
    nodes <- read.csv(paste0(filename, "_nodes.csv"), stringsAsFactors = FALSE)
    edges <- read.csv(paste0(filename, "_edges.csv"), stringsAsFactors = FALSE)
    c_net_from_edgelist(edges, nodes) -> go
  } else if (format == "graphml") {
    if (!grepl("\\.graphml$", filename)) filename <- paste0(filename, ".graphml")
    igraph::read_graph(paste0(filename, ".graphml"), format = "graphml") -> go
    go <- c_net_update(go, initialize = TRUE)
  } else {
    if (!grepl(paste0("\\.", format), filename)) filename <- paste0(filename, ".", format)
    igraph::read_graph(filename, format = format) -> go
    go <- c_net_update(go, initialize = TRUE)
  }
  go
}

#' Summaries two columns information
#' @param df data.frame
#' @param from first column name or index
#' @param to second column name or index
#' @param count (optional) weight column, if no, each equal to 1
#' @param direct consider direct? default: FALSE
#'
#' @return data.frame
#' @export
#' @examples
#' test <- data.frame(
#'   a = sample(letters[1:4], 10, replace = TRUE),
#'   b = sample(letters[1:4], 10, replace = TRUE)
#' )
#' summ_2col(test, direct = TRUE)
#' summ_2col(test, direct = FALSE)
#' if (requireNamespace("circlize")) {
#'   summ_2col(test, direct = TRUE) %>% pcutils::my_circo()
#' }
summ_2col <- function(df, from = 1, to = 2, count = 3, direct = FALSE) {
  if (ncol(df) < 2) stop("need at least two columns")
  if (ncol(df) == 2) {
    tmp <- cbind(df, count = 1)
  } else {
    tmp <- dplyr::select(df, !!from, !!to, !!count)
  }
  cols <- colnames(tmp)
  colnames(tmp) <- c("from", "to", "count")

  if (direct) {
    tmp <- (dplyr::group_by(tmp, from, to) %>% dplyr::summarise(count = sum(count)))
    colnames(tmp) <- cols
    return(as.data.frame(tmp))
  }

  com <- \(group1, group2, levels){
    factor(c(group1, group2), levels = levels) %>% sort()
  }

  group <- factor(c(tmp[, 1], tmp[, 2]))
  tmp1 <- apply(tmp, 1, function(x) com(x[1], x[2], levels(group))) %>%
    t() %>%
    as.data.frame()

  tmp1 <- cbind(tmp1, tmp$count)
  colnames(tmp1) <- c("from", "to", "count")
  tmp1 <- dplyr::group_by(tmp1, from, to) %>% dplyr::summarise(count = sum(count))
  colnames(tmp1) <- cols
  return(as.data.frame(tmp1))
}
