# ========3.layout========
#' Layout coordinates
#'
#' @param go igraph or metanet
#' @param method
#'  (1) as_line(), as_arc(), as_polygon(), as_polyarc(), as_polycircle(), as_circle_tree();
#'  (2) as_star(), as_tree(), in_circle(), nicely(), on_grid(), on_sphere(),randomly(), with_dh(), with_fr(), with_gem(), with_graphopt(), with_kk(),with_lgl(), with_mds(),. see \code{\link[igraph]{layout_}};
#'  (3) a character, "auto","backbone","centrality","circlepack","dendrogram",
#'  "eigen","focus","hive","igraph","linear","manual","matrix",
#'  "partition","pmds","stress","treemap","unrooted". see \code{\link[ggraph]{create_layout}}
#' @param order_by order nodes according to a node attribute
#' @param order_ls manual the discrete variable with a vector, or continuous variable with "desc" to decreasing
#' @param seed random seed
#' @param line_curved consider line curved, only for some layout methods like as_line(), as_polygon().default:0
#' @param ... add
#' @aliases c_net_lay
#' @family layout
#' @return coors object: coordinates for nodes, columns: name, X, Y; curved for edges, columns: from, to, curved;
#' @export
#' @examples
#' library(igraph)
#' c_net_layout(co_net) -> coors
#' c_net_plot(co_net, coors)
#' c_net_plot(co_net, c_net_layout(co_net, in_circle()), vertex.size = 2)
#' c_net_plot(co_net, c_net_layout(co_net, in_circle(), order_by = "v_class"), vertex.size = 2)
#' c_net_plot(co_net, c_net_layout(co_net, in_circle(), order_by = "size", order_ls = "desc"))
#' c_net_plot(co_net, c_net_layout(co_net, as_polygon(3)))
c_net_layout <- function(go, method = igraph::nicely(), order_by = NULL, order_ls = NULL,
                         seed = 1234, line_curved = 0.5, ...) {
  set.seed(seed)
  name <- x <- y <- NULL
  if ("igraph_layout_spec" %in% class(method)) {
    coors <- igraph::layout_(go, method)
  } else if ("poly" %in% class(method)) {
    coors <- method(go, group2 = order_by, group2_order = order_ls)
  } else if ("layout" %in% class(method)) {
    coors <- method(go)
  } else if (is.character(method)) {
    message("Use method from `ggraph`: ", method)
    lib_ps("ggraph", library = FALSE)
    data <- ggraph::create_layout(clean_igraph(go), layout = method, ...)
    coors <- data %>% dplyr::select(name, x, y)
    colnames(coors) <- c("name", "X", "Y")
  } else {
    stop("No valid method")
  }
  if (inherits(coors, "coors")) {
    return(coors)
  }

  # order
  if (is.matrix(coors)) {
    get_v(go) -> tmp_v
    coors <- order_tmp_v_name(tmp_v, coors, order_by = order_by, order_ls = order_ls)
  }

  curved <- NULL
  # if line type, need to consider edge.curved
  if ("line" %in% class(method)) {
    tmp_e <- data.frame(igraph::as_data_frame(go))[, c("from", "to")]
    if (nrow(tmp_e) > 0) {
      curved <- data.frame(tmp_e, curved = line_curved, row.names = NULL)
    }
  }
  coors <- data.frame(coors, row.names = NULL)
  coors <- structure(list(coors = coors, curved = curved), class = "coors")
  return(coors)
}

order_tmp_v_name <- function(tmp_v, coors_mat, order_by = NULL, order_ls = NULL) {
  if (is.null(order_by)) {
    coors <- data.frame(name = tmp_v$name, X = coors_mat[, 1], Y = coors_mat[, 2], row.names = NULL)
  } else {
    ordervec <- tmp_v[, order_by]
    if (is.numeric(ordervec)) {
      name <- tmp_v[order(ordervec, decreasing = is.null(order_ls)), "name"]
    } else {
      ordervec <- pcutils::change_fac_lev(ordervec, order_ls)
      name <- tmp_v[order(ordervec), "name"]
    }
    coors <- data.frame(name = name, X = coors_mat[, 1], Y = coors_mat[, 2], row.names = NULL)
  }
  return(coors)
}

is_layout <- \(x){
  any(class(x) %in% c("igraph_layout_spec", "layout"))
}

get_coors <- \(coors, go, ...){
  edge_curved <- NULL
  # 1.如果coors是NULL，去graph_attr找一下，没有的话就默认nicely计算
  if (is.null(coors)) {
    if (is.null(igraph::graph_attr(go, "coors"))) {
      coors <- c_net_layout(go, igraph::nicely(), ...)
    } else {
      coors <- igraph::graph_attr(go, "coors")
    }
  }
  # 2.如果是layout函数，那就计算
  if (is_layout(coors)) coors <- c_net_layout(go, coors, ...)
  # 3.如果是一个提供的coors对象（list），那就把curved传到edge里，coors导出为df
  if (inherits(coors, "coors")) {
    if (!is.null(coors$curved)) {
      edge_curved <- dplyr::left_join(get_e(go), coors$curved, by = c("from", "to"), suffix = c(".x", "")) %>%
        dplyr::select("from", "to", "curved")
      edge_curved[is.na(edge_curved)] <- 0
      edge_curved <- data.frame(edge_curved, row.names = NULL)
    }
    if (!is.null(coors$coors)) {
      coors <- coors$coors
    } else {
      coors <- c_net_layout(go, igraph::nicely(), ...)
    }
  }
  # 4.如果是matrix，变成df
  if (is.data.frame(coors)) {
    if (!"name" %in% colnames(coors)) coors <- as.matrix(coors)
  }
  if (is.matrix(coors)) coors <- data.frame(name = V(go)$name, X = coors[, 1], Y = coors[, 2], row.names = NULL)
  # 5.如果是df了，那就对齐name用于下一步的绘图，
  if (is.data.frame(coors)) {
    coors <- coors[match(V(go)$name, coors$name), ]
    return(structure(list(coors = coors, curved = edge_curved), class = "coors"))
  }
  stop("coors wrong!")
}

combine_coors <- function(..., list = NULL) {
  name <- from <- to <- NULL
  list <- c(list(...), list)
  if (!all(vapply(list, \(i)inherits(i, "coors"), logical(1)))) stop("some input are not coors object")
  coors <- lapply(list, \(i)i[["coors"]]) %>% do.call(rbind, .)
  curved <- lapply(list, \(i)i[["curved"]]) %>% do.call(rbind, .)
  if (any(duplicated(coors$name))) {
    warning("some duplicated name in coors$coors")
    coors <- dplyr::distinct(coors, name, .keep_all = TRUE)
  }
  if (!is.null(curved)) {
    curved2 <- dplyr::distinct(curved, from, to, .keep_all = TRUE)
    if (nrow(curved2) != nrow(curved)) warning("some duplicates in coors$curved")
    curved <- curved2
  }
  coors <- structure(list(coors = coors, curved = curved), class = "coors")
  return(coors)
}

#' Layout as a line
#'
#' @param angle anticlockwise rotation angle
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#' @family layout
#' @examples
#' as_line()(co_net)
#' c_net_plot(co_net, coors = as_line(pi / 2))
as_line <- function(angle = 0) {
  fun <- \(go){
    nv <- length(V(go))
    data.frame(
      x = seq(-cos(angle), cos(angle), len = nv),
      y = seq(-sin(angle), sin(angle), len = nv)
    ) %>%
      as.matrix() %>%
      round(., 4)
  }
  class(fun) <- c("line", "layout", "function")
  fun
}

#' Layout as a arc
#'
#' @param angle anticlockwise rotation angle
#' @param arc the radian of arc
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#' @family layout
#' @examples
#' as_arc()(co_net)
#' c_net_plot(co_net, coors = as_arc(pi / 2), rescale = FALSE)
as_arc <- function(angle = 0, arc = pi) {
  fun <- \(go){
    # (0,0) is the midpoint of circle
    nv <- length(V(go))
    theta <- seq(-arc / 2 + angle, arc / 2 + angle, len = nv)
    coor <- data.frame(x = cos(theta), y = sin(theta))
    as.matrix(coor) %>% round(., 4)
  }
  class(fun) <- c("layout", "function")
  fun
}

#' Layout as a polygon
#'
#' @param n how many edges of this polygon
#' @param line_curved line_curved 0~0.5
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#' @family layout
#' @examples
#' as_polygon()(co_net)
as_polygon <- function(n = 3, line_curved = 0.5) {
  fun <- \(go, group2 = NULL, group2_order = NULL){
    V(go)$poly_group <- rep(paste0("omic", seq_len(n)), len = length(go))
    if (n < 2) stop("n should bigger than 1")
    g_layout_polygon(go,
      group = "poly_group", group2 = group2, group2_order = group2_order,
      line_curved = line_curved
    ) -> oridata
    oridata
  }
  class(fun) <- c("poly", "layout", "function")
  fun
}

#' Layout as a polyarc
#'
#' @param n how many arcs of this poly_arc
#' @param space the space between each arc, default: pi/3
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#' @family layout
#' @examples
#' as_polyarc()(co_net)
as_polyarc <- \(n = 3, space = pi / 3){
  fun <- \(go, group2 = NULL, group2_order = NULL){
    V(go)$poly_group <- rep(paste0("omic", 1:n), len = length(go))
    if (n < 2) stop("n should bigger than 1")
    g_layout_polyarc(go, "poly_group", space = space, group2 = group2, group2_order = group2_order) -> oridata
    oridata
  }
  class(fun) <- c("poly", "layout", "function")
  fun
}

#' Layout as a polycircle
#'
#' @param n how many circles of this polycircle
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#' @family layout
#' @examples
#' as_polycircle()(co_net)
as_polycircle <- \(n = 2){
  fun <- \(go, group2 = NULL, group2_order = NULL){
    V(go)$poly_group <- rep(paste0("omic", 1:n), len = length(go))
    if (n < 2) stop("n should bigger than 1")
    g_layout_polycircle(go, "poly_group", group2 = group2, group2_order = group2_order) -> oridata
    oridata
  }
  class(fun) <- c("poly", "layout", "function")
  fun
}

#' Layout as a circle_tree
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#' @family layout
as_circle_tree <- \(){
  fun <- \(go){
    name <- x <- y <- NULL
    lib_ps("ggraph", library = FALSE)
    data <- ggraph::create_layout(clean_igraph(go), layout = "igraph", algorithm = "tree", circular = TRUE)
    coors <- data %>% dplyr::select(name, x, y)
    colnames(coors) <- c("name", "X", "Y")
    coors
  }
  class(fun) <- c("layout", "function")
  fun
}

big_layout <- \(zoom1, layout1, nodeGroup){
  c_net_update(igraph::make_ring(nlevels(nodeGroup$group))) -> tmp_da_net
  da <- get_coors(layout1, tmp_da_net)[["coors"]]
  da$name <- NULL

  if (!all(levels(nodeGroup$group) %in% rownames(da))) rownames(da) <- levels(nodeGroup$group)

  # center of each group
  scale_f <- (ceiling(max(da) - min(da)))
  scale_f <- ifelse(scale_f == 0, 1, scale_f / 2)
  da <- da / scale_f * zoom1
  colnames(da) <- c("X", "Y")
  return(da)
}

#' Layout with group
#'
#' @param go igraph or metanet object
#' @param group group name (default: module)
#' @param zoom1 big network layout size
#' @param zoom2 average sub_network layout size, or numeric vector, or "auto"
#' @param layout1 layout1 method, one of
#'      (1) a dataframe or matrix: rowname is group, two columns are X and Y
#'      (2) function: layout method for \code{\link{c_net_layout}} default: in_circle()
#' @param layout2 one of functions: layout method for \code{\link{c_net_layout}}, or a list of functions.
#' @param show_big_layout show the big layout to help you adjust.
#' @param ... add
#' @param group_order group_order
#'
#' @return coors
#' @export
#' @family g_layout
#' @examples
#' \donttest{
#' data("c_net")
#' module_detect(co_net, method = "cluster_fast_greedy") -> co_net_modu
#' g_layout(co_net_modu, group = "module", zoom1 = 30, zoom2 = "auto", layout2 = as_line()) -> oridata
#' plot(co_net_modu, coors = oridata)
#' }
g_layout <- function(go, group = "module", group_order = NULL, layout1 = in_circle(), zoom1 = 20, layout2 = in_circle(),
                     zoom2 = 3, show_big_layout = FALSE, ...) {
  name <- ID <- NULL

  stopifnot(is_igraph(go))
  if (!group %in% igraph::vertex_attr_names(go)) stop("no group named ", group, " !")
  get_v(go) %>% dplyr::select(name, !!group) -> nodeGroup
  colnames(nodeGroup) <- c("ID", "group")
  nodeGroup$group <- factor(nodeGroup$group)
  if (!is.null(group_order)) nodeGroup$group <- pcutils::change_fac_lev(nodeGroup$group, group_order)

  da <- big_layout(zoom1, layout1, nodeGroup)
  if (show_big_layout) {
    message("Big layout:")
    print(da)
    graphics::plot.new()
    plot(da, pch = 21, bty = "n", bg = get_cols(nrow(da), "col2"), main = "Big layout coordinates")
    graphics::text(da$X * 0.8, da$Y * 0.9, rownames(da))
    return(invisible())
  }

  # layout of vertexes in one group
  {
    layoutls <- list()
    if (is_layout(layout2)) {
      layoutls <- rep(list(layout2), nlevels(nodeGroup$group))
    } else if (all(class(layout2) == "list")) {
      if (is.null(names(layout2))) {
        layoutls <- rep(layout2, len = nlevels(nodeGroup$group))
      } else {
        for (i in levels(nodeGroup$group)) {
          if (i %in% names(layout2)) {
            layoutls[[i]] <- layout2[[i]]
          } else {
            layoutls[[i]] <- igraph::in_circle()
            warning("layout of ", i, " not set, use in_circle()")
          }
        }
      }
    }
    if (is.null(names(layoutls))) names(layoutls) <- levels(nodeGroup$group)

    zoom2 <- rep(zoom2, nlevels(nodeGroup$group))
    if (zoom2[1] == "auto") zoom2 <- ceiling((table(nodeGroup$group))^(1 / 3))
    names(zoom2) <- levels(nodeGroup$group)
  }

  # get coors
  all_coors <- setNames(vector("list", nlevels(nodeGroup$group)), levels(nodeGroup$group))
  for (i in levels(nodeGroup$group)) {
    nodeGroup[nodeGroup[, "group"] == i, "ID"] -> tmpid
    igraph::subgraph(go, tmpid) -> tmp_net

    get_coors(layoutls[[i]], tmp_net, ...) -> coors
    data <- coors$coors
    if ("igraph_layout_spec" %in% class(layoutls[[i]])) {
      data[, c("X", "Y")] <- igraph::norm_coords(as.matrix(data[, c("X", "Y")]))
    }
    data[, "X"] <- data[, "X"] * zoom2[i] + da[i, "X"]
    data[, "Y"] <- data[, "Y"] * zoom2[i] + da[i, "Y"]

    coors$coors <- data
    all_coors[[i]] <- coors
  }
  coors <- combine_coors(list = all_coors)
  return(coors)
}

#' Layout with group as a polygon
#'
#' @param go igraph
#' @param group group name (default:v_group)
#' @param group2 group2 name, will order nodes in each group according to group2_order
#' @param group2_order group2_order
#' @param line_curved line_curved 0~1
#' @param group_order group_order
#'
#' @return coors
#' @export
#' @family g_layout
#' @examples
#' g_layout_polygon(multi1) -> oridata
#' c_net_plot(multi1, oridata)
#' g_layout_polyarc(multi1, group2 = "v_class", group2_order = c(LETTERS[4:1])) -> oridata
#' c_net_plot(multi1, oridata)
#' g_layout_polycircle(co_net2, group2 = "v_class") -> oridata
#' c_net_plot(co_net2, oridata)
g_layout_polygon <- function(go, group = "v_group", group_order = NULL, group2 = NULL, group2_order = NULL, line_curved = 0.5) {
  n <- length(unique(igraph::vertex.attributes(go)[[group]]))

  if (n < 2) stop("n should bigger than 1")
  angle_ls <- -pi / 2 + (seq(0, n - 1, 1)) * 2 * pi / n
  fun_ls <- lapply(angle_ls, \(i)as_line(i))

  g_layout(go,
    group_order = group_order,
    group = group, zoom1 = 1, zoom2 = 0.9 * (ifelse(n > 2, tan(pi / n), 2)),
    layout2 = fun_ls, order_by = group2, order_ls = group2_order
  ) -> oridata

  if (is.data.frame(oridata$curved)) oridata$curved$curved <- line_curved
  oridata
}


#' Layout with group as a polyarc
#'
#' @param space the space between each arc, default: pi/4
#' @param scale_node_num scale with the node number in each group
#'
#' @rdname g_layout_polygon
#' @export
g_layout_polyarc <- function(go, group = "v_group", group_order = NULL,
                             group2 = NULL, group2_order = NULL, space = pi / 4, scale_node_num = TRUE) {
  get_v(go) -> tmp_v
  group1 <- as.factor(tmp_v[, group])
  n <- nlevels(group1)
  if (n < 2) stop("n should bigger than 1")

  if (!is.null(group_order)) group1 <- pcutils::change_fac_lev(group1, group_order)
  # consider each group numbers!!!
  g_num <- table(group1)
  sep <- space / n

  if (scale_node_num) {
    arc_r <- (2 * pi - space) * as.numeric(g_num) / length(group1)
  } else {
    arc_r <- rep((2 * pi - space) / n, n)
  }

  names(arc_r) <- levels(group1)

  # coordinate
  coors <- data.frame()
  theta1 <- 0
  for (i in names(arc_r)) {
    tmp_t <- seq(theta1, theta1 + arc_r[i], len = g_num[i])
    tmp_v1 <- tmp_v[tmp_v[, group] == i, ]
    tmp_coor <- order_tmp_v_name(tmp_v1, data.frame(X = cos(tmp_t), Y = sin(tmp_t)), group2, group2_order)

    coors <- rbind(coors, tmp_coor)
    theta1 <- theta1 + arc_r[i] + sep
  }
  coors
}


#' Layout with group as a polyarc
#'
#' @param space the space between each arc, default: pi/4
#' @param scale_node_num scale with the node number in each group
#'
#' @rdname g_layout_polygon
#' @export
g_layout_polycircle <- function(go, group = "v_group", group_order = NULL, group2 = NULL, group2_order = NULL) {
  name <- NULL
  n <- length(unique(igraph::vertex.attributes(go)[[group]]))
  if (n < 2) stop("n should bigger than 1")
  get_v(go) %>% dplyr::select(name, !!group) -> nodeGroup

  if (is.null(group_order)) {
    group_order <- table(nodeGroup[, group]) %>%
      sort() %>%
      names()
  }

  g_layout(go,
    group_order = group_order,
    group = group, layout1 = matrix(0, nrow = n, ncol = 2),
    zoom1 = 1, zoom2 = 1:n,
    layout2 = igraph::in_circle(), order_by = group2, order_ls = group2_order
  ) -> oridata
  oridata
}

#' Layout with group nicely
#'
#' @param go igraph or metanet
#' @param group group name (default: module)
#' @param mode circlepack, treemap, backbone, stress
#' @param ... add
#'
#' @export
#'
#' @rdname g_layout
#'
#' @examples
#' \donttest{
#' data("c_net")
#' module_detect(co_net, method = "cluster_fast_greedy") -> co_net_modu
#' if (requireNamespace("ggraph")) {
#'   plot(co_net_modu, coors = g_layout_nice(co_net_modu, group = "module"))
#'   plot(co_net_modu, coors = g_layout_nice(co_net_modu, group = "module", mode = "treemap"))
#' }
#' }
g_layout_nice <- function(go, group = "module", mode = "circlepack", ...) {
  name <- leaf <- x <- y <- NULL
  lib_ps("ggraph", library = FALSE)
  stopifnot(is_igraph(go))

  mode <- match.arg(mode, c("circlepack", "treemap", "backbone", "stress"))

  if (!group %in% vertex_attr_names(go)) stop("no group named ", group, " !")
  get_v(go) %>% dplyr::select(name, !!group) -> nodeGroup
  colnames(nodeGroup) <- c("ID", "group")
  nodeGroup$group <- as.factor(nodeGroup$group)

  edge <- data.frame(from = paste("group_", nodeGroup$group, sep = ""), to = nodeGroup$ID)

  directed <- TRUE
  if (mode %in% c("backbone")) directed <- FALSE

  mygraph <- igraph::graph_from_data_frame(edge, directed = directed)
  data <- ggraph::create_layout(mygraph, layout = mode, ...)
  coor <- data %>% dplyr::select(name, x, y)
  colnames(coor) <- c("name", "X", "Y")

  return(structure(list(coors = coor, curved = NULL), class = "coors"))
}

#' Transform a dataframe to a network edgelist.
#'
#' @param test df
#' @param fun default: sum
#'
#' @return metanet
#' @export
#'
#' @examples
#' data("otutab", package = "pcutils")
#' cbind(taxonomy, num = rowSums(otutab))[1:20, ] -> test
#' df2net_tree(test) -> ttt
#' plot(ttt)
#' plot(ttt, coors = as_circle_tree())
df2net_tree <- function(test, fun = sum) {
  flag <- FALSE
  if (!is.numeric(test[, ncol(test)])) {
    test$num <- 1
  } else {
    flag <- TRUE
    name <- colnames(test)[ncol(test)]
  }
  nc <- ncol(test)
  if (nc < 3) stop("as least 3-columns dataframe")

  link <- pcutils::df2link(test, fun = fun)

  nodes <- link$nodes
  links <- link$links
  if (flag) {
    colnames(links)[3] <- colnames(nodes)[3] <- name
  } else {
    name <- "weight"
  }

  # c_net_from_edgelist(as.data.frame(links),vertex = nodes)
  net <- igraph::graph_from_data_frame(as.data.frame(links), vertices = nodes)
  net <- c_net_update(net)
  net <- c_net_set(net, vertex_class = "level", vertex_size = name, edge_width = name)
  graph_attr(net, "coors") <- c_net_layout(net, as_tree())
  net
}

# ========4.plot========

#' Plot a metanet
#'
#' @param x metanet object
#' @param ... add
#'
#' @return plot
#' @exportS3Method
#' @method plot metanet
plot.metanet <- function(x, ...) {
  # 具有n_type的metanet的默认参数记录在这里，用于快速绘图。
  go <- x
  if (is.null(get_n(go)$n_type)) {
    c_net_plot(go, ...)
  } else if (get_n(go)$n_type == "skeleton") {
    skeleton_plot(go, ...)
  } else if (get_n(go)$n_type == "module") {
    default_arg <- list(
      labels_num = 0,
      group_legend_title = "Module"
    )
    do.call(c_net_plot, append(list(go = go), pcutils::update_param(default_arg, list(...))))
  } else if (get_n(go)$n_type == "venn") {
    nice_size <- ceiling(60 / sqrt(length(V(go)))) + 1
    default_arg <- list(
      labels_num = "all",
      vertex_size_range = list("Group" = c(1.5 * nice_size, 1.5 * nice_size), "elements" = c(0.5 * nice_size, 0.5 * nice_size)),
      vertex.shape = "circle",
      legend = FALSE, edge.curved = 0.3,
      edge.color = unique(V(go)$color)
    )
    do.call(c_net_plot, append(list(go = go), pcutils::update_param(default_arg, list(...))))
  } else if (get_n(go)$n_type == "twocol") {
    nice_size <- ceiling(60 / sqrt(length(V(go)))) + 1
    default_arg <- list(
      labels_num = 0,
      vertex.shape = "circle",
      edge_legend = FALSE,
      edge.color = "black"
    )
    do.call(c_net_plot, append(list(go = go), pcutils::update_param(default_arg, list(...))))
  } else if (get_n(go)$n_type == "ko_net") {
    nice_size <- ceiling(60 / sqrt(length(V(go)))) + 1
    default_arg <- list(
      labels_num = "all",
      vertex.shape = "circle",
      vertex_size_range = list("Pathway" = c(1.2 * nice_size, 1.2 * nice_size), "KOs" = c(0.6 * nice_size, 0.6 * nice_size)),
      edge_legend = FALSE,
      edge.color = "black",
      mark_alpha = 0.1
    )
    do.call(c_net_plot, append(list(go = go), pcutils::update_param(default_arg, list(...))))
  } else {
    c_net_plot(go, ...)
  }
}

get_net_main <- function(n_index) {
  main <- "Network"
  if (!is.null(n_index$n_type)) {
    switch(n_index$n_type,
      "single" = {
        main <- "Correlation network"
      },
      "bipartite" = {
        main <- "Bipartite network"
      },
      "multi_full" = {
        main <- "Multi-omics network"
      },
      "module" = {
        main <- paste0(n_index$n_modules, "-modules network")
      },
      "skeleton" = {
        main <- paste0(n_index$skeleton, " skeleton network")
      },
      "venn" = {
        main <- "Venn network"
      },
      default = {
        main <- "Network"
      }
    )
  }
  return(main)
}

scale_size_width <- function(tmp_v, tmp_e, vertex_size_range, edge_width_range) {
  {        v_groups <- unique(tmp_v$v_group)
    nice_size <- ceiling(60 / sqrt(nrow(tmp_v))) + 1

    vertex_size_range_default <- rep(list(c(max(nice_size * 0.4, 3), min(nice_size * 1.6, 12))), length(v_groups))
    names(vertex_size_range_default) <- v_groups

    if (!is.null(vertex_size_range)) {
      if (!is.list(vertex_size_range)) vertex_size_range <- list(vertex_size_range)
      if (is.null(names(vertex_size_range))) {
        vertex_size_range <- rep(vertex_size_range, length(v_groups))
        names(vertex_size_range) <- v_groups
      }
      vertex_size_range <- pcutils::update_param(vertex_size_range_default, vertex_size_range)
    } else {
      vertex_size_range <- vertex_size_range_default
    }

    node_size_text <- setNames(as.list(numeric(length(v_groups))), v_groups)
    for (i in v_groups) {
      node_size_text[[i]] <- c(
        min(tmp_v[tmp_v$v_group == i, "size"], na.rm = TRUE),
        max(tmp_v[tmp_v$v_group == i, "size"], na.rm = TRUE)
      )
      tmp_v[tmp_v$v_group == i, "size"] <- do.call(pcutils::mmscale, append(
        list(tmp_v[tmp_v$v_group == i, "size"]),
        as.list(vertex_size_range[[i]][1:2])
      ))
    }    }

  {
    edge_width_range_default <- vertex_size_range_default[[1]] / 6
    if (is.null(edge_width_range)) edge_width_range <- edge_width_range_default
    edge_width_text <- c(min(tmp_e$width, na.rm = TRUE), max(tmp_e$width, na.rm = TRUE))
    tmp_e$width <- pcutils::mmscale(tmp_e$width, edge_width_range[1], edge_width_range[2])
  }

  envir <- parent.frame()
  assign("node_size_text", node_size_text, envir)
  assign("edge_width_text", edge_width_text, envir)
  assign("tmp_e", tmp_e, envir)
  assign("tmp_v", tmp_v, envir)
}

some_custom_paras <- function(tmp_v, tmp_e, ...) {
  params <- list(...)
  params_name <- names(params)

  tmp_v$vertex.label.color <- "black"
  if ("vertex.size" %in% params_name) tmp_v$size <- params[["vertex.size"]]
  if ("vertex.color" %in% params_name) {
    tmp_v$color <- condance(data.frame(
      tmp_v$color,
      pcutils::tidai(tmp_v$v_class, params[["vertex.color"]])
    ))
  }
  if ("vertex.shape" %in% params_name) {
    tmp_v$shape <- condance(data.frame(
      tmp_v$shape,
      pcutils::tidai(tmp_v$v_group, params[["vertex.shape"]])
    ))
  }
  if ("vertex.label" %in% params_name) tmp_v$label <- params[["vertex.label"]]
  if ("vertex.label.color" %in% params_name) {
    tmp_v$vertex.label.color <- condance(data.frame(
      "black",
      pcutils::tidai(tmp_v$v_group, params[["vertex.label.color"]])
    ))
  }

  if ("edge.color" %in% params_name) {
    tmp_e$color <- condance(data.frame(
      tmp_e$color,
      pcutils::tidai(tmp_e$e_type, params[["edge.color"]])
    ))
  }
  if ("edge.lty" %in% params_name) {
    tmp_e$lty <- condance(data.frame(
      tmp_e$lty,
      pcutils::tidai(tmp_e$e_class, params[["edge.lty"]])
    ))
  }
  if ("edge.width" %in% params_name) tmp_e$width <- params[["edge.width"]]

  envir <- parent.frame()
  assign("tmp_e", tmp_e, envir)
  assign("tmp_v", tmp_v, envir)
}

get_show_labels <- function(tmp_v, labels_num) {
  name <- size <- color <- e_type <- lty <- e_class <- v_class <- shape <- NULL
  {
    if (labels_num == "all") {
      tmp_v %>% dplyr::pull(name) -> toplabel
    } else {
      if (labels_num >= 1) {
        tmp_v %>%
          dplyr::top_n(labels_num, size) %>%
          dplyr::pull(name) %>%
          head(labels_num) -> toplabel
      } else {
        tmp_v %>%
          dplyr::top_frac(labels_num, size) %>%
          dplyr::pull(name) %>%
          head(ceiling(labels_num * nrow(tmp_v))) -> toplabel
      }
    }
    tmp_v$label <- ifelse(tmp_v$name %in% toplabel, tmp_v$label, NA)
  }
  return(tmp_v)
}

module_set_for_plot <- function(tmp_v, mark_module, mark_color) {
  if (mark_module) {
    new_modu <- as_module(setNames(tmp_v$module, tmp_v$name))
    new_modu[["others"]] <- NULL

    module_color <- pcutils::get_cols(length(new_modu))
    if (!is.null(mark_color)) module_color <- condance(data.frame(module_color, pcutils::tidai(names(new_modu), mark_color)))

    module_color <- setNames(module_color, names(new_modu))
    module_color <- module_color[names(module_color) != "others"]
  } else {
    new_modu <- module_color <- NULL
  }
  envir <- parent.frame()
  assign("new_modu", new_modu, envir)
  assign("module_color", module_color, envir)
}

get_module_coors <- function(go = NULL, coors = NULL, tmp_v = NULL, ori_coors = NULL, module_label_just = c(0.5, 0.5), rescale_flag = TRUE) {
  name <- size <- color <- e_type <- lty <- e_class <- v_class <- shape <- NULL
  X <- Y <- module <- minx <- maxx <- miny <- maxy <- NULL
  if (is.null(go)) {
    if (is.null(tmp_v) & is.null(ori_coors)) message("input `tmp_v` and `ori_coors` when `go` is null.")
  } else {
    tmp_v <- get_v(go)
    ori_coors <- get_coors(coors, go)
    coors <- ori_coors$coors[, c("X", "Y")] %>% as.matrix()
  }
  module_coors <- dplyr::left_join(tmp_v[, c("name", "module")], ori_coors$coors, by = "name")
  if (rescale_flag) module_coors <- dplyr::mutate(module_coors, X = mmscale(X, -1, 1), Y = mmscale(Y, -1, 1))
  module_coors <- dplyr::group_by(module_coors, module) %>%
    dplyr::summarise(minx = min(X), maxx = max(X), miny = min(Y), maxy = max(Y))
  module_label_just <- rep(module_label_just, 2)
  module_coors <- mutate(module_coors,
    X = minx + module_label_just[1] * (maxx - minx),
    Y = miny + module_label_just[2] * (maxy - miny)
  )
  return(module_coors)
}

produce_c_net_legends <- function(tmp_v, tmp_e,
                                  legend_position, legend_number, legend_cex,
                                  node_size_text, edge_width_text,
                                  group_legend_title, group_legend_order,
                                  color_legend, color_legend_order,
                                  size_legend, size_legend_title,
                                  edge_legend, edge_legend_title, edge_legend_order,
                                  width_legend, width_legend_title,
                                  lty_legend, lty_legend_title, lty_legend_order, ...) {
  name <- size <- color <- e_type <- lty <- e_class <- v_class <- shape <- left_leg_x <- right_leg_x <- NULL

  legend_position_default <- c(left_leg_x = -2, left_leg_y = 1, right_leg_x = 1.2, right_leg_y = 1)

  if (is.null(legend_position)) legend_position <- legend_position_default
  if (is.null(names(legend_position))) {
    legend_position <- setNames(legend_position, names(legend_position_default)[seq_along(legend_position)])
  }
  legend_position <- pcutils::update_param(legend_position_default, legend_position)
  here_env <- environment()
  lapply(names(legend_position), \(i){
    assign(i, legend_position[i], here_env)
  })

  vgroups <- pcutils::change_fac_lev(tmp_v$v_group, group_legend_order)
  vgroups <- levels(vgroups)

  if (color_legend) {
    pchls <- c("circle" = 21, "square" = 22)

    if (is.null(group_legend_title)) {
      group_legend_title <- setNames(vgroups, vgroups)
    } else if (is.null(names(group_legend_title))) {
      group_legend_title <- setNames(rep(group_legend_title, len = length(vgroups)), vgroups)
    }

    for (g_i in vgroups) {
      tmp_v1 <- tmp_v[tmp_v$v_group == g_i, c("v_class", "color", "shape")]

      tmp_v1$v_class <- factor(tmp_v1$v_class, levels = stringr::str_sort(unique(tmp_v1$v_class), numeric = TRUE))

      vclass <- pcutils::change_fac_lev(tmp_v1$v_class, color_legend_order)
      vclass <- levels(vclass)

      node_cols <- dplyr::distinct(tmp_v1, color, v_class)
      node_cols <- setNames(node_cols$color, node_cols$v_class)
      node_shapes <- dplyr::distinct(tmp_v1, shape, v_class)
      node_shapes <- setNames(node_shapes$shape, node_shapes$v_class)

      if (legend_number) {
        eee <- table(tmp_v1$v_class)
        le_text <- paste(vclass, eee[vclass], sep = ": ")
      } else {
        le_text <- vclass
      }
      if (length(le_text) == 0) le_text <- ""
      legend(left_leg_x, left_leg_y,
        cex = 0.7 * legend_cex, adj = 0,
        legend = le_text, title.cex = 0.8 * legend_cex,
        title = group_legend_title[g_i], title.font = 2, title.adj = 0,
        col = "black", pt.bg = node_cols[vclass], bty = "n", pch = pchls[node_shapes[vclass]]
      )

      left_leg_y <- left_leg_y - (length(vclass) * 0.12 + 0.2) * legend_cex
    }
  }

  if (size_legend) {
    legend(
      x = right_leg_x, y = right_leg_y,
      cex = 0.7 * legend_cex, title.font = 2, title = size_legend_title, title.adj = 0,
      legend = c(
        paste(lapply(node_size_text[vgroups], \(i)round(i[1], 3)), collapse = "/ "),
        paste(lapply(node_size_text[vgroups], \(i)round(i[2], 3)), collapse = "/ ")
      ),
      adj = 0, title.cex = 0.8 * legend_cex,
      col = "black", bty = "n", pch = 21, pt.cex = c(min(tmp_v$size), max(tmp_v$size)) * legend_cex / 5
    )
    right_leg_y <- right_leg_y - 0.5 * legend_cex
  }

  if (edge_legend) {
    edges <- pcutils::change_fac_lev(tmp_e$e_type, edge_legend_order)
    edges <- levels(edges)
    edge_cols <- dplyr::distinct(tmp_e, color, e_type)
    edge_cols <- setNames(edge_cols$color, edge_cols$e_type)
    if (legend_number) {
      eee <- table(tmp_e$e_type)
      le_text <- paste(edges, eee[edges], sep = ": ")
    } else {
      le_text <- edges
    }
    legend(right_leg_x, right_leg_y,
      cex = 0.7 * legend_cex, title.font = 2, title = edge_legend_title, title.adj = 0,
      legend = le_text, adj = 0, title.cex = 0.8 * legend_cex,
      col = edge_cols[edges], bty = "n", lty = 1
    )
    right_leg_y <- right_leg_y - (length(unique(tmp_e$color)) * 0.12 + 0.2) * legend_cex
  }

  if (width_legend) {
    legend(right_leg_x, right_leg_y,
      cex = 0.7 * legend_cex, title.font = 2, title = width_legend_title, title.adj = 0,
      legend = edge_width_text %>% round(., 3), adj = 0, title.cex = 0.8 * legend_cex,
      col = "black", bty = "n", lty = 1, lwd = c(min(tmp_e$width), max(tmp_e$width))
    )
    right_leg_y <- right_leg_y - 0.5 * legend_cex
  }

  if (lty_legend) {
    edges <- pcutils::change_fac_lev(tmp_e$e_class, lty_legend_order)
    edges <- levels(edges)
    edge_ltys <- dplyr::distinct(tmp_e, lty, e_class)
    edge_ltys <- setNames(edge_ltys$lty, edge_ltys$e_class)

    if (legend_number) {
      eee <- table(tmp_e$e_class)
      le_text <- paste(edges, eee[edges], sep = ": ")
    } else {
      le_text <- edges
    }
    legend(right_leg_x, right_leg_y,
      cex = 0.7 * legend_cex, title.font = 2, title = lty_legend_title, title.adj = 0,
      legend = le_text, adj = 0, title.cex = 0.8 * legend_cex,
      col = "black", bty = "n", lty = edge_ltys[edges]
    )
  }
}

#' Plot a metanet
#'
#' @param go an igraph or metanet object
#' @param coors the coordinates you saved
#' @param ... additional parameters for \code{\link[igraph]{igraph.plotting}}
#' @param labels_num show how many labels,>1 indicates number, <1 indicates fraction, "all" indicates all, default:5
#' @param vertex_size_range the vertex size range, e.g. c(1,10)
#' @param edge_width_range the edge width range, e.g. c(1,10)
#'
#' @param plot_module logical, plot module?
#' @param mark_module logical, mark the modules?
#' @param mark_color mark colors
#' @param mark_alpha mark fill alpha, default 0.3
#' @param module_label module_label
#' @param module_label_cex module_label_cex
#' @param module_label_color module_label_color
#' @param module_label_just module_label_just
#'
#' @param legend all legends
#' @param legend_number legend with numbers
#' @param legend_cex character expansion factor relative to current par("cex"), default: 1
#' @param legend_position legend_position, default: c(left_leg_x=-1.9,left_leg_y=1,right_leg_x=1.2,right_leg_y=1)
#'
#' @param group_legend_title group_legend_title, length must same to the numbers of v_group
#' @param group_legend_order group_legend_order vector
#' @param color_legend logical
#' @param color_legend_order color_legend_order vector
#' @param size_legend logical
#' @param size_legend_title size_legend_title
#'
#' @param edge_legend logical
#' @param edge_legend_title edge_legend_title
#' @param edge_legend_order edge_legend_order vector, e.g. c("positive","negative")
#' @param width_legend logical
#' @param width_legend_title width_legend_title
#'
#' @param lty_legend logical
#' @param lty_legend_title lty_legend_title
#' @param lty_legend_order lty_legend_order
#'
#' @param seed random seed, default:1234, make sure each plot is the same.
#'
#' @family plot
#' @return a network plot
#' @export
#'
#' @examples
#' data("c_net")
#' c_net_plot(co_net)
#' c_net_plot(co_net2)
#' c_net_plot(multi1)
c_net_plot <- function(go, coors = NULL, ..., labels_num = 5,
                       vertex_size_range = NULL, edge_width_range = NULL,
                       plot_module = FALSE,
                       mark_module = FALSE, mark_color = NULL, mark_alpha = 0.3,
                       module_label = FALSE, module_label_cex = 2, module_label_color = "black",
                       module_label_just = c(0.5, 0.5),
                       legend = TRUE, legend_number = FALSE, legend_cex = 1,
                       legend_position = c(left_leg_x = -2, left_leg_y = 1, right_leg_x = 1.2, right_leg_y = 1),
                       group_legend_title = NULL, group_legend_order = NULL,
                       color_legend = TRUE, color_legend_order = NULL,
                       size_legend = FALSE, size_legend_title = "Node Size",
                       edge_legend = TRUE, edge_legend_title = "Edge type", edge_legend_order = NULL,
                       width_legend = FALSE, width_legend_title = "Edge width",
                       lty_legend = FALSE, lty_legend_title = "Edge class", lty_legend_order = NULL,
                       seed = 1234) {
  name <- size <- color <- e_type <- lty <- e_class <- v_class <- shape <- NULL
  new_modu <- module_color <- node_size_text <- edge_width_text <- NULL
  set.seed(seed)

  if (!"metanet" %in% class(go)) go <- c_net_update(go)
  # modules
  if (plot_module) {
    go <- to_module_net(go)
    if (is.null(group_legend_title)) group_legend_title <- "Module"
  }

  # get network type
  get_net_main(get_n(go)) -> main
  get_v(go) -> tmp_v
  get_e(go) -> tmp_e

  # get coordinates
  ori_coors <- get_coors(coors, go, seed = seed)
  coors <- ori_coors$coors[, c("X", "Y")] %>% as.matrix()
  if (is.null(ori_coors$curved)) {
    edge_curved <- 0
  } else {
    edge_curved <- ori_coors$curved$curved
  }

  # scale the size and width
  scale_size_width(tmp_v, tmp_e, vertex_size_range, edge_width_range)

  # some custom parameters
  some_custom_paras(tmp_v, tmp_e, ...)

  # show labels
  tmp_v <- get_show_labels(tmp_v, labels_num)

  # modules set
  module_set_for_plot(tmp_v, mark_module, mark_color)

  # main plot
  {
    old_xpd <- graphics::par(mar = c(4, 2, 2, 2), xpd = TRUE)
    on.exit(graphics::par(old_xpd), add = TRUE)

    igraph::plot.igraph(go,
      layout = coors,
      vertex.size = tmp_v$size,
      vertex.color = tmp_v$color,
      vertex.shape = tmp_v$shape,
      vertex.label.color = tmp_v$vertex.label.color,
      edge.color = tmp_e$color,
      edge.lty = tmp_e$lty,
      edge.width = tmp_e$width,
      mark.groups = new_modu,
      mark.col = pcutils::add_alpha(module_color[names(new_modu)], mark_alpha),
      mark.border = module_color[names(new_modu)],
      ...,
      main = main,
      vertex.label.font = 1,
      vertex.label.cex = 0.07 * tmp_v$size,
      vertex.label = tmp_v$label,
      edge.arrow.size = 0.3,
      edge.arrow.width = 0.5,
      edge.curved = edge_curved,
      margin = c(0, 0, 0, 0)
    )
  }

  # add module_label
  if (module_label) {
    rescale_flag <- TRUE
    params <- list(...)
    if ("rescale" %in% names(params)) {
      if (!params[["rescale"]]) rescale_flag <- FALSE
    }
    module_coors <- get_module_coors(
      tmp_v = tmp_v, ori_coors = ori_coors,
      module_label_just = module_label_just, rescale_flag = rescale_flag
    )

    n_module <- nrow(module_coors)
    module_label_cex <- rep(module_label_cex, n_module)
    module_label_color <- rep(module_label_color, n_module)
    for (i in seq_len(n_module)) {
      text(
        x = module_coors[i, "X"], y = module_coors[i, "Y"],
        labels = module_coors[i, "module"],
        cex = module_label_cex, col = module_label_color[i]
      )
    }
  }

  if (!legend) {
    return(invisible())
  }

  # produce legends
  produce_c_net_legends(
    tmp_v, tmp_e,
    legend_position, legend_number, legend_cex,
    node_size_text, edge_width_text,
    group_legend_title, group_legend_order,
    color_legend, color_legend_order,
    size_legend, size_legend_title,
    edge_legend, edge_legend_title, edge_legend_order,
    width_legend, width_legend_title,
    lty_legend, lty_legend_title, lty_legend_order, ...
  )
}

#' Transfer an igraph object to a ggig
#'
#' @param go igraph or meatnet
#' @param coors coordinates for nodes,columns: name, X, Y
#'
#' @return ggig object
#' @export
#' @family plot
#' @examples
#' as.ggig(co_net, coors = c_net_layout(co_net)) -> ggig
#' plot(ggig)
#' as.ggig(multi1, coors = c_net_layout(multi1)) -> ggig
#' plot(ggig, labels_num = 0.3)
as.ggig <- function(go, coors = NULL) {
  list(n_index = get_n(go), v_index = get_v(go), e_index = get_e(go)) -> net_par_res

  if (is.null(coors)) coors <- c_net_layout(go)
  coors <- get_coors(coors, go)
  coors <- coors$coors

  # add coors
  coors <- coors[, 1:3] %>% na.omit()
  net_par_res$v_index %<>% dplyr::left_join(., coors, by = "name", suffix = c("", ".1"))
  net_par_res$e_index %<>% dplyr::left_join(., coors, by = c("from" = "name")) %>%
    dplyr::rename(X1 = "X", Y1 = "Y") %>%
    dplyr::left_join(., coors, by = c("to" = "name")) %>%
    dplyr::rename(X2 = "X", Y2 = "Y")

  class(net_par_res) <- c("ggig", "list")
  return(net_par_res)
}

#' Plot a ggig
#'
#' @param x ggig object
#' @inheritParams c_net_plot
#'
#' @family plot
#' @return ggplot
#' @exportS3Method
#' @method plot ggig
plot.ggig <- function(x, coors = NULL, ..., labels_num = 5,
                      vertex_size_range = NULL, edge_width_range = NULL,
                      plot_module = FALSE,
                      mark_module = FALSE, mark_color = NULL, mark_alpha = 0.3,
                      module_label = FALSE, module_label_cex = 2, module_label_color = "black",
                      module_label_just = c(0.5, 0.5),
                      legend_number = FALSE, legend = TRUE, legend_cex = 1,
                      legend_position = c(left_leg_x = -2, left_leg_y = 1, right_leg_x = 1.2, right_leg_y = 1),
                      group_legend_title = NULL, group_legend_order = NULL,
                      color_legend = TRUE, color_legend_order = NULL,
                      size_legend = FALSE, size_legend_title = "Node Size",
                      edge_legend = TRUE, edge_legend_title = "Edge type", edge_legend_order = NULL,
                      width_legend = FALSE, width_legend_title = "Edge width",
                      lty_legend = FALSE, lty_legend_title = "Edge class", lty_legend_order = NULL,
                      seed = 1234) {
  rename <- size <- name <- color <- e_type <- lty <- e_class <- v_class <- shape <- X1 <- Y1 <- X2 <- Y2 <- width <- X <- Y <- label <- NULL
  edge_width_text <- NULL

  ggig <- x
  ggig$v_index -> tmp_v
  ggig$e_index -> tmp_e

  set.seed(seed)
  # get coordinates
  if (!is.null(coors)) {
    tmp_v$X <- tmp_v$Y <- NULL
    tmp_e$X1 <- tmp_e$X2 <- tmp_e$Y1 <- tmp_e$Y2 <- NULL
    # add coors
    tmp_v %<>% dplyr::left_join(., coors, by = "name", suffix = c("", ".1"))
    tmp_e %<>% dplyr::left_join(., coors, by = c("from" = "name")) %>%
      rename(X1 = "X", Y1 = "Y") %>%
      dplyr::left_join(., coors, by = c("to" = "name")) %>%
      rename(X2 = "X", Y2 = "Y")
  }

  # get network type
  main <- get_net_main(ggig$n_index)

  # scale the size and width
  scale_size_width(tmp_v, tmp_e, vertex_size_range, edge_width_range)

  # new shapes
  tmp_v$shape <- tidai(tmp_v$v_group, 21:26)

  # some custom parameters
  some_custom_paras(tmp_v, tmp_e, ...)

  # show labels
  tmp_v <- get_show_labels(tmp_v, labels_num)

  if (TRUE) {
    tmp_e$e_type <- pcutils::change_fac_lev(tmp_e$e_type, edge_legend_order)
    edges <- levels(tmp_e$e_type)
    edge_cols <- dplyr::distinct(tmp_e, color, e_type)
    edge_cols <- setNames(edge_cols$color, edge_cols$e_type)
    if (legend_number) {
      eee <- table(tmp_e$e_type)
      edge_text <- paste(edges, eee[edges], sep = ": ")
    } else {
      edge_text <- edges
    }
  }

  if (TRUE) {
    edges1 <- levels(factor(tmp_e$e_class))
    edge_ltys <- dplyr::distinct(tmp_e, lty, e_class)
    edge_ltys <- setNames(edge_ltys$lty, edge_ltys$e_class)

    if (legend_number) {
      eee <- table(tmp_e$e_class)
      lty_text <- paste(edges1, eee[edges1], sep = ": ")
    } else {
      lty_text <- edges1
    }
  }

  if (TRUE) {
    pchls <- c("circle" = 21, "square" = 22)

    vgroups <- pcutils::change_fac_lev(tmp_v$v_group, group_legend_order)

    node_size_text <- c(
      paste(lapply(node_size_text[levels(vgroups)], \(i)round(i[1], 3)), collapse = "/ "),
      paste(lapply(node_size_text[levels(vgroups)], \(i)round(i[2], 3)), collapse = "/ ")
    )

    new_f <- c()
    for (g_i in levels(vgroups)) {
      tmp_v1 <- tmp_v[tmp_v$v_group == g_i, c("v_class", "color", "shape")]
      tmp_f <- pcutils::change_fac_lev(tmp_v1$v_class, color_legend_order)
      new_f <- c(new_f, levels(tmp_f))
    }

    tmp_v$v_class <- pcutils::change_fac_lev(tmp_v$v_class, new_f)
    vclass <- levels(tmp_v$v_class)

    node_cols <- dplyr::distinct(tmp_v, color, v_class)
    node_cols <- setNames(node_cols$color, node_cols$v_class)
    node_shapes <- dplyr::distinct(tmp_v, shape, v_class)
    node_shapes <- setNames(node_shapes$shape, node_shapes$v_class)

    if (legend_number) {
      eee <- table(tmp_v$v_class)
      le_text <- paste(vclass, eee[vclass], sep = ": ")
    } else {
      le_text <- vclass
    }
  }

  p <- ggplot() +
    geom_segment(aes(
      x = X1, y = Y1, xend = X2, yend = Y2, color = e_type,
      linewidth = width, linetype = e_class
    ), data = tmp_e, alpha = 0.7) + # draw edges
    scale_color_manual(
      name = edge_legend_title, values = edge_cols,
      label = edge_text, guide = ifelse(edge_legend, "legend", "none")
    ) + # edge colors
    scale_linetype_manual(
      name = lty_legend_title, values = edge_ltys,
      label = lty_text, guide = ifelse(lty_legend, "legend", "none")
    ) + # edge linetype
    scale_linewidth(
      name = width_legend_title, breaks = c(min(tmp_e$width), max(tmp_e$width)), range = c(0.5, 1),
      labels = edge_width_text, guide = ifelse(width_legend, "legend", "none")
    )

  p1 <- p +
    geom_point(aes(X, Y, fill = v_class, size = size, shape = v_class), data = tmp_v) + # draw nodes
    # scale_shape_manual(values =setNames(pchls[node_shapes],vclass))+#node shape
    scale_shape_manual(values = node_shapes) +
    scale_fill_manual(
      name = group_legend_title, values = node_cols[vclass],
      labels = le_text, guide = ifelse(color_legend, "legend", "none")
    ) + # node color
    scale_size(
      name = size_legend_title, breaks = c(min(tmp_v$size), max(tmp_v$size)),
      labels = node_size_text, guide = ifelse(size_legend, "legend", "none")
    ) + # node size

    ggnewscale::new_scale("size") +
    geom_text(aes(X, Y, size = size, label = label), col = "black", data = tmp_v, show.legend = FALSE) +
    scale_size(range = c(1, 3), guide = "none") +
    guides(
      fill = guide_legend(override.aes = list(shape = node_shapes[vclass])),
      shape = "none"
    )

  p2 <- p1 + labs(title = main) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    coord_fixed(ratio = 1) +
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(
      legend.background = element_rect(colour = NA),
      legend.box.background = element_rect(colour = NA),
      legend.key = element_rect(fill = NA)
    ) +
    theme(panel.background = element_rect(fill = "white", colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

  if (!legend) {
    return(p2 + theme(legend.position = "none"))
  }

  p2
}


#' Input a graphml file exported by Gephi
#'
#' @param file graphml file exported by Gephi
#' @family plot
#' @return list contains the igraph object and coordinates
#'
#' @export
input_gephi <- function(file) {
  X <- Y <- code <- NULL
  igraph::read.graph(file, format = "graphml") -> gephi
  get_v(gephi) -> tmp_v
  # extract coors
  if (!all(c("x", "y", "r", "g", "b", "id") %in% colnames(tmp_v))) {
    stop("This file is not exported by Gephi, please use igraph::read.graph()")
  }
  coors <- tmp_v[, c("x", "y")]
  coors <- data.frame(name = tmp_v$label, X = coors[, 1], Y = coors[, 2])
  coors %>% dplyr::mutate(X = pcutils::mmscale(X, -40, 40), Y = pcutils::mmscale(Y, -40, 40)) -> coors

  # transform color
  pcutils::rgb2code(tmp_v[, c("r", "g", "b")]) %>% dplyr::pull(code) -> tmp_v$color
  E(gephi)$color <- ifelse(E(gephi)$cor > 0, "#48A4F0", "#E85D5D")
  # scale size
  tmp_v$size <- pcutils::mmscale(tmp_v$size, 1, 5)
  E(gephi)$width <- pcutils::mmscale(E(gephi)$width, 0.05, 0.2)
  # delete
  tmp_v %>%
    dplyr::select(-c("label", "x", "y", "r", "g", "b", "id")) %>%
    as.list() -> vertex.attributes(gephi)
  edge.attributes(gephi)["Edge Label"] <- edge.attributes(gephi)["id"] <- NULL

  gephi <- c_net_update(gephi)
  igraph::graph_attr(gephi, "coors") <- coors
  return(list(go = gephi, coors = coors))
}

#' plot use networkD3
#'
#' @param go metanet
#' @param v_class which attributes use to be v_class
#' @param ... see \code{\link[networkD3]{forceNetwork}}
#' @return D3 plot
#' @export
#' @family plot
#' @examples
#' data("c_net")
#' plot(co_net2)
#' if (requireNamespace("networkD3")) {
#'   netD3plot(co_net2)
#' }
netD3plot <- function(go, v_class = "v_class", ...) {
  flag <- "y"
  if (length(V(go)) > 200) {
    message("Too big network, recommend using Gephi to layout,still use networkD3?")
    flag <- readline("yes/no(y/n):")
  }
  if (tolower(flag) %in% c("yes", "y")) {
    lib_ps("networkD3", library = FALSE)
    go <- c_net_set(go, vertex_class = v_class)
    get_v(go) -> tmp_v
    nodes <- tmp_v[, c("name", "v_class", "size", "color")]
    colnames(nodes) <- c("name", "group", "size", "color")
    nodes$size <- pcutils::mmscale(nodes$size, 2, 40)

    colors <- unique(nodes$color)

    get_e(go) -> tmp_e
    links <- tmp_e[, c("from", "to", "width", "color")]
    links$width <- pcutils::mmscale(links$width, 0.5, 1.5)
    # give ids
    links$IDsource <- match(links$from, nodes$name) - 1
    links$IDtarget <- match(links$to, nodes$name) - 1
    # Create force directed network plot
    networkD3::forceNetwork(
      Links = links, Nodes = nodes,
      Source = "IDsource", Target = "IDtarget", linkColour = links$color, linkDistance = 20,
      linkWidth = networkD3::JS("function(d) { return (d.width); }"), charge = -5,
      NodeID = "name", Group = "group", Nodesize = "size",
      colourScale = networkD3::JS(paste0("d3.scaleOrdinal([`", paste(colors, collapse = "`,`"), "`])")), legend = TRUE, ...
    )
  }
}

MetaNet_theme <- {
  ggplot2::theme_classic(base_size = 13) + ggplot2::theme(
    axis.text = element_text(color = "black"),
    plot.margin = grid::unit(rep(0.5, 4), "lines"),
    strip.background = ggplot2::element_rect(fill = NA)
  )
}

#' Venn network
#'
#' @param tab data.frame (row is elements, column is group), or a list (names is group, value is elements)
#'
#' @return plot
#' @export
#' @family plot
#' @examples
#' data(otutab, package = "pcutils")
#' tab <- otutab[400:485, 1:3]
#' venn_net(tab) -> v_net
#' plot(v_net)
venn_net <- function(tab) {
  # pcutils:::venn_cal(tab)->vennlist
  tab[is.na(tab)] <- 0
  edgelist <- data.frame()
  if (is.data.frame(tab)) {
    groupss <- colnames(tab)
    for (i in groupss) {
      if (sum(tab[, i] > 0) > 0) edgelist <- rbind(edgelist, data.frame(Group = i, elements = rownames(tab)[tab[, i] > 0]))
    }
  } else if (all(class(tab) == "list")) {
    vennlist <- tab
    groupss <- names(vennlist)
    for (i in groupss) {
      if (length(vennlist[[i]] > 0)) edgelist <- rbind(edgelist, data.frame(Group = i, elements = vennlist[[i]]))
    }
  } else {
    stop("wrong input tab")
  }

  nodelist <- rbind(
    data.frame(name = groupss, v_group = "Group", v_class = paste0("Group: ", groupss)),
    data.frame(name = unique(edgelist$elements), v_group = "elements", v_class = "elements")
  )
  venn_net <- c_net_from_edgelist(edgelist, vertex = nodelist)
  graph.attributes(venn_net)$n_type <- "venn"
  all_group <- get_e(venn_net)[, c("from", "to")] %>%
    pcutils::squash("from") %>%
    dplyr::rename(name = "to", all_group = "from")
  venn_net <- c_net_set(venn_net, all_group, vertex_class = "all_group", edge_type = "from")
  venn_net
}


#' Quick build a metanet from two columns table
#'
#' @param edgelist two columns table (no elements exist in two columns at same time)
#'
#' @return metanet
#' @export
#' @family plot
#' @examples
#' twocol <- data.frame(
#'   "col1" = sample(letters, 30, replace = TRUE),
#'   "col2" = sample(c("A", "B"), 30, replace = TRUE)
#' )
#' twocol_net <- twocol_edgelist(twocol)
#' plot(twocol_net)
#' c_net_plot(twocol_net, g_layout_polygon(twocol_net))
twocol_edgelist <- function(edgelist) {
  if (any(edgelist[, 1] %in% edgelist[, 2])) stop("Must no elements exist in two columns at same time")
  nodelist <- rbind(
    data.frame(name = unique(edgelist[, 1]), v_group = names(edgelist)[1], v_class = names(edgelist)[1]),
    data.frame(name = unique(edgelist[, 2]), v_group = names(edgelist)[2], v_class = names(edgelist)[2])
  )
  venn_net <- c_net_from_edgelist(edgelist, vertex = nodelist)
  graph.attributes(venn_net)$n_type <- "twocol"
  # venn_net=c_net_set(venn_net,edge_type = "from")
  venn_net
}

#' Plot olympic rings using network
#'
#' @return network plot
#' @export
#' @family plot
#' @examples
#' olympic_rings_net()
olympic_rings_net <- function() {
  r <- 1
  pensize <- r / 6
  rings_data <- data.frame(
    x = c(-2 * (r + pensize), -(r + pensize), 0, (r + pensize), 2 * (r + pensize)),
    y = c(r, 0, r, 0, r),
    color = c("#0081C8", "#FCB131", "#000000", "#00A651", "#EE334E")
  )
  g1 <- module_net(module_number = 5, n_node_in_module = 30)
  plot(g1,
    coors = g_layout(g1, layout1 = rings_data[, 1:2], zoom1 = 1.2, zoom2 = 0.5),
    rescale = FALSE, legend = FALSE, main = "Olympic Rings", vertex.frame.color = NA,
    edge.width = 0, vertex.color = setNames(rings_data$color, 1:5), vertex.size = 9
  )
}
