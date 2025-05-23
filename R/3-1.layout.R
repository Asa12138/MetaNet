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
#' @param rescale logical, scale the X, Y to (-1,1)
#' @param ... add
#'
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
                         seed = 1234, line_curved = 0.5, rescale = TRUE, ...) {
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
    if (rescale) {
      coors <- rescale_coors(coors)
    }
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
  coors <- as_coors(coors, curved)
  if (rescale) {
    coors <- rescale_coors(coors)
  }
  return(coors)
}

order_circle_func <- function(coors_mat) {
  dis <- X <- Y <- NULL
  mid_x <- mean(coors_mat[, 1])
  mid_y <- mean(coors_mat[, 2])
  coors_mat <- data.frame(X = coors_mat[, 1], Y = coors_mat[, 2]) %>%
    dplyr::mutate(dis = sqrt((X - mid_x)^2 + (Y - mid_y)^2)) %>%
    dplyr::arrange(dis) %>%
    dplyr::select_at(1:2) %>%
    as.matrix()
  coors_mat
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


#' Transfer to a coors object
#'
#' @param coors data.frame
#' @param curved line curved
#'
#' @return coors object
#' @export
#'
as_coors <- function(coors, curved = NULL) {
  if (inherits(coors, "coors")) {
    if (!is.null(curved)) attributes(coors)$curved <- curved
    return(coors)
  }
  if (!identical(colnames(coors), c("name", "X", "Y"))) {
    stop("The input should be a data frame with columns: name, X, Y")
  }
  attributes(coors)$class <- c("coors", class(coors))
  attributes(coors)$curved <- curved
  return(coors)
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
    if (!is.null(attributes(coors)$curved)) {
      edge_curved <- dplyr::left_join(get_e(go), attributes(coors)$curved, by = c("from", "to"), suffix = c(".x", "")) %>%
        dplyr::select("from", "to", "curved")
      edge_curved[is.na(edge_curved)] <- 0
      edge_curved <- data.frame(edge_curved, row.names = NULL)
    }
  }
  # 4.如果是matrix，变成df
  if (is.data.frame(coors)) {
    if (!"name" %in% colnames(coors)) coors <- as.matrix(coors)
  }
  if (is.matrix(coors)) coors <- data.frame(name = V(go)$name, X = coors[, 1], Y = coors[, 2], row.names = NULL)
  # 5.如果是df了，那就对齐name用于下一步的绘图，
  if (is.data.frame(coors)) {
    coors <- data.frame(coors, check.names = FALSE)
    setdiff(V(go)$name, coors$name) -> tmp_name
    if (length(setdiff(V(go)$name, coors$name)) > 0) {
      warning("Some nodes are not in coors, will be randomly assigned")
      coors <- rbind(
        coors,
        data.frame(
          name = tmp_name,
          X = runif(length(tmp_name), range(coors$X)[1], range(coors$X)[2]),
          Y = runif(length(tmp_name), range(coors$Y)[1], range(coors$Y)[2])
        )
      )
    }
    coors <- coors[match(V(go)$name, coors$name), ]
    coors <- as_coors(coors, edge_curved)
    return(coors)
  }
  stop("coors wrong!")
}

rescale_coors <- function(coors, keep_asp = TRUE) {
  stopifnot(inherits(coors, "coors"))
  if (keep_asp) {
    diff_x <- diff(range(coors$X))
    diff_y <- diff(range(coors$Y))
    if (diff_x > diff_y) {
      coors$X <- mmscale(coors$X, -1, 1)
      coors$Y <- mmscale(coors$Y, -1 * diff_y / diff_x, 1 * diff_y / diff_x)
    } else {
      coors$X <- mmscale(coors$X, -1 * diff_x / diff_y, 1 * diff_x / diff_y)
      coors$Y <- mmscale(coors$Y, -1, 1)
    }
  } else {
    coors$X <- mmscale(coors$X, -1, 1)
    coors$Y <- mmscale(coors$Y, -1, 1)
  }
  coors
}

combine_coors <- function(..., list = NULL) {
  name <- from <- to <- NULL
  list <- c(list(...), list)
  if (!all(vapply(list, \(i)inherits(i, "coors"), logical(1)))) stop("some input are not coors object")
  coors <- lapply(list, \(i)i) %>% do.call(rbind, .)
  curved <- lapply(list, \(i)attributes(i)$curved) %>% do.call(rbind, .)
  if (any(duplicated(coors$name))) {
    warning("some duplicated name in coors")
    coors <- dplyr::distinct(coors, name, .keep_all = TRUE)
  }
  if (!is.null(curved)) {
    curved2 <- dplyr::distinct(curved, from, to, .keep_all = TRUE)
    if (nrow(curved2) != nrow(curved)) warning("some duplicates in attributes(i)$curved")
    curved <- curved2
  }
  coors <- as_coors(coors, curved)
  return(coors)
}

#' Transform the layout of a 'coors' object
#'
#' This function applies various transformations to a 'coors' object, including
#' scaling, aspect ratio adjustment, rotation, mirroring, and pseudo-3D perspective transformation.
#'
#' @param coors An object of class 'coors', containing node coordinates.
#' @param scale A numeric value to scale the layout (default = 1).
#' @param aspect_ratio A numeric value to adjust the Y-axis scaling (default = 1).
#' @param rotation A numeric value in degrees to rotate the layout (default = 0).
#' @param mirror_x A logical value indicating whether to mirror along the X-axis (default = FALSE).
#' @param mirror_y A logical value indicating whether to mirror along the Y-axis (default = FALSE).
#' @param shear_x A numeric value to apply a shear transformation in the X direction (default = 0).
#' @param shear_y A numeric value to apply a shear transformation in the Y direction (default = 0).
#'
#' @return A transformed 'coors' object with updated coordinates.
#' @export
transform_coors <- function(coors, scale = 1, aspect_ratio = 1,
                            rotation = 0, mirror_x = FALSE, mirror_y = FALSE,
                            shear_x = 0, shear_y = 0) {
  stopifnot(inherits(coors, "coors"))

  # 复制原始数据
  new_coor <- coors

  mid_x <- mean(new_coor$X)
  mid_y <- mean(new_coor$Y)

  new_coor$X <- new_coor$X - mid_x
  new_coor$Y <- new_coor$Y - mid_y

  # 放大/缩小
  new_coor$X <- new_coor$X * scale
  new_coor$Y <- new_coor$Y * scale * aspect_ratio

  # 旋转（角度转换为弧度）
  theta <- rotation * pi / 180
  new_x <- new_coor$X * cos(theta) - new_coor$Y * sin(theta)
  new_y <- new_coor$X * sin(theta) + new_coor$Y * cos(theta)
  new_coor$X <- new_x
  new_coor$Y <- new_y

  # 透视投影（shear变换）
  new_coor$X <- new_coor$X + shear_x * new_coor$Y
  new_coor$Y <- new_coor$Y + shear_y * new_coor$X

  # 镜像变换
  if (mirror_x) new_coor$X <- -new_coor$X
  if (mirror_y) new_coor$Y <- -new_coor$Y
  # 还原到原始位置
  new_coor$X <- new_coor$X + mid_x
  new_coor$Y <- new_coor$Y + mid_y

  coors <- new_coor
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
#' c_net_plot(co_net, coors = as_arc(pi / 2))
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
as_polycircle <- \(n = 5){
  fun <- \(go, group2 = NULL, group2_order = NULL){
    V(go)$poly_group <- rep(paste0("omic", 1:n), distribute_points_in_concentric_circles(length(go), n))
    if (n < 2) stop("n should bigger than 1")
    g_layout_polycircle(go, "poly_group", group2 = group2, group2_order = group2_order) -> oridata
    oridata
  }
  class(fun) <- c("poly", "layout", "function")
  fun
}


distribute_points_in_concentric_circles <- function(n, k = NULL, mode = c("linear", "perimeter")) {
  mode <- match.arg(mode)

  if (is.null(k)) {
    # 自动决定层数：sqrt(n) 是个常用经验
    k <- ceiling(sqrt(n))
  }

  layer_weights <- switch(mode,
    linear = 1:k, # 权重 1, 2, ..., k
    perimeter = 1:k # 半径与层号成比例，周长也是
  )

  # 归一化
  weight_sum <- sum(layer_weights)
  raw_counts <- n * layer_weights / weight_sum

  # 四舍五入并微调使总和为n
  counts <- round(raw_counts)
  diff_n <- n - sum(counts)
  if (diff_n != 0) {
    # 按最大余数或最小绝对值误差来调整
    adjust_order <- order(abs(raw_counts - counts), decreasing = TRUE)
    for (i in seq_len(abs(diff_n))) {
      idx <- adjust_order[i]
      counts[idx] <- counts[idx] + sign(diff_n)
    }
  }

  return(counts)
}

#' Layout as a multi_layer
#'
#' @param n how many arcs of this multi_layer
#' @param layout see method in \code{\link[MetaNet]{c_net_layout}}
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#' @family layout
#' @examples
#' as_multi_layer()(co_net)
as_multi_layer <- \(n = 3, layout = on_grid()){
  fun <- \(go, group2 = NULL, group2_order = NULL){
    V(go)$poly_group <- rep(paste0("omic", 1:n), len = length(go))
    if (n < 2) stop("n should bigger than 1")
    g_layout_multi_layer(go, layout = layout, group = "poly_group", group2 = group2, group2_order = group2_order) -> oridata
    oridata
  }
  class(fun) <- c("poly", "layout", "function")
  fun
}

#' Layout as a multi_layer
#'
#' @param n how many arcs of this multi_layer
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#' @family layout
as_poly_sector <- \(n = 3){
  fun <- \(go, group2 = NULL, group2_order = NULL){
    V(go)$poly_group <- rep(paste0("omic", 1:n), len = length(go))
    if (n < 2) stop("n should bigger than 1")
    g_layout_poly_sector(go, group = "poly_group", group2 = group2, group2_order = group2_order) -> oridata
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
  tmp_g <- igraph::make_ring(nlevels(nodeGroup$group))
  V(tmp_g)$name <- levels(nodeGroup$group)
  c_net_update(tmp_g) -> tmp_da_net
  da <- get_coors(layout1, tmp_da_net)
  rownames(da) <- da$name
  da <- da[, c("X", "Y")]

  # center of each group
  scale_f <- (ceiling(max(da) - min(da)))
  scale_f <- ifelse(scale_f == 0, 1, scale_f / 2)
  da <- da / scale_f * zoom1
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
#' @param rescale logical, scale the X, Y to (-1,1)
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
                     zoom2 = 3, show_big_layout = FALSE, rescale = TRUE, ...) {
  name <- NULL

  stopifnot(is_igraph(go))
  if (!group %in% igraph::vertex_attr_names(go)) stop("no group named ", group, " !")
  get_v(go) %>% dplyr::select(name, !!group) -> nodeGroup
  colnames(nodeGroup) <- c("ID", "group")
  nodeGroup$group <- factor(nodeGroup$group)
  stopifnot(nlevels(nodeGroup$group) > 1)

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
    # print(i)
    nodeGroup[nodeGroup[, "group"] == i, "ID"] -> tmpid
    igraph::subgraph(go, tmpid) -> tmp_net

    get_coors(layoutls[[i]], tmp_net, ...) -> coors
    data <- coors
    if (nrow(data) == 1) {
      data <- data.frame(name = data$name, X = 0, Y = 0) %>% as_coors()
    }
    if ("igraph_layout_spec" %in% class(layoutls[[i]])) {
      data[, c("X", "Y")] <- igraph::norm_coords(as.matrix(data[, c("X", "Y")]))
    }
    data[, "X"] <- data[, "X"] * zoom2[i] + da[i, "X"]
    data[, "Y"] <- data[, "Y"] * zoom2[i] + da[i, "Y"]

    coors <- data
    all_coors[[i]] <- coors
  }
  coors <- combine_coors(list = all_coors)
  if (rescale) {
    coors <- rescale_coors(coors)
  }
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
#' g_layout_multi_layer(co_net2, group2 = "v_class") -> oridata
#' c_net_plot(co_net2, oridata)
g_layout_polygon <- function(go, group = "v_group", group_order = NULL, group2 = NULL, group2_order = NULL, line_curved = 0.5) {
  n <- length(unique(igraph::vertex.attributes(go)[[group]]))

  if (n < 2) stop("n should bigger than 1")
  `_internal_group` <- NULL
  get_v(go) %>% dplyr::pull(!!group) -> group1
  group1_level <- table(group1)
  V(go)$`_internal_group` <- group1

  angle_ls <- -pi / 2 + (seq(0, n - 1, 1)) * 2 * pi / n
  names(angle_ls) <- pcutils::change_fac_lev(names(group1_level), group_order) %>% levels()

  layout2_ls <- list()
  for (i in names(group1_level)) {
    c_net_filter(go, `_internal_group` == i) %>%
      c_net_layout(
        method = as_line(angle_ls[i]), order_by = group2,
        order_ls = group2_order, rescale = FALSE
      ) -> layout2_ls[[i]]
  }

  g_layout(go,
    group_order = group_order,
    group = group, zoom1 = 1, zoom2 = 0.9 * (ifelse(n > 2, tan(pi / n), 2)),
    layout2 = layout2_ls, order_by = group2, order_ls = group2_order
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
                             group2 = NULL, group2_order = NULL, space = pi / 4,
                             scale_node_num = TRUE) {
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
  coors <- as_coors(coors)
  coors <- rescale_coors(coors)
  coors
}


#' Layout with group as a polycircle
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

#' Layout with group as a multi_layer
#'
#' @param layout see method in \code{\link[MetaNet]{c_net_layout}}
#'
#' @rdname g_layout_polygon
#' @export
g_layout_multi_layer <- function(go, layout = igraph::in_circle(), group = "v_group", group_order = NULL,
                                 group2 = NULL, group2_order = NULL, scale_node_num = TRUE) {
  n <- length(unique(igraph::vertex.attributes(go)[[group]]))
  if (n < 2) stop("n should bigger than 1")

  get_v(go) %>%
    dplyr::pull(!!group) %>%
    as.factor() -> group1
  if (!is.null(group_order)) group1 <- pcutils::change_fac_lev(group1, group_order)
  group1_level <- table(group1)
  `_internal_group` <- NULL
  V(go)$`_internal_group` <- group1

  # layout of vertexes in one group
  {
    layoutls <- list()
    if (is_layout(layout)) {
      layoutls <- rep(list(layout), nlevels(group1))
    } else if (all(class(layout) == "list")) {
      if (is.null(names(layout))) {
        layoutls <- rep(layout, len = nlevels(group1))
      } else {
        for (i in levels(group1)) {
          if (i %in% names(layout)) {
            layoutls[[i]] <- layout[[i]]
          } else {
            layoutls[[i]] <- igraph::in_circle()
            warning("layout of ", i, " not set, use in_circle()")
          }
        }
      }
    }
    if (is.null(names(layoutls))) names(layoutls) <- levels(group1)
  }

  layout2_ls <- list()
  for (i in names(group1_level)) {
    c_net_filter(go, `_internal_group` == i) %>%
      c_net_layout(method = layoutls[[i]], order_by = group2, order_ls = group2_order) %>%
      transform_coors(shear_x = 3, aspect_ratio = 0.2) -> layout2_ls[[i]]
  }

  if (scale_node_num) {
    zoom2 <- group1_level / mean(group1_level)
  } else {
    zoom2 <- rep(1, nlevels(group1))
  }

  g_layout(go,
    group_order = group_order,
    group = group, layout1 = data.frame(X = 0, Y = seq(-1, 1, length = length(group1_level))),
    zoom1 = 1, zoom2 = zoom2,
    layout2 = layout2_ls
  ) -> oridata
  oridata
}

#' Layout with group nicely
#'
#' @param go igraph or metanet
#' @param group group name (default: module)
#' @param mode circlepack, treemap, backbone, stress
#' @param ... add
#' @param group_zoom zoom for each group
#'
#' @export
#' @return coors
#' @family g_layout
#' @examples
#' \donttest{
#' data("c_net")
#' module_detect(co_net, method = "cluster_fast_greedy") -> co_net_modu
#' if (requireNamespace("ggraph")) {
#'   plot(co_net_modu, coors = g_layout_nice(co_net_modu, group = "module"))
#'   plot(co_net_modu, coors = g_layout_nice(co_net_modu, group = "module", mode = "treemap"))
#' }
#' }
g_layout_nice <- function(go, group = "module", mode = "circlepack", group_zoom = 1, ...) {
  name <- x <- y <- NULL
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
  coors <- data %>% dplyr::select(name, x, y)
  colnames(coors) <- c("name", "X", "Y")
  coors <- as_coors(coors)

  for (i in unique(nodeGroup$group)) {
    tmp_id <- nodeGroup[nodeGroup[, "group"] == i, "ID"]
    transform_coors(coors[coors$name %in% tmp_id, ], scale = group_zoom) -> coors[coors$name %in% tmp_id, ]
  }

  coors <- rescale_coors(coors)
  return(coors)
}

#' @rdname g_layout_nice
#' @export
get_big_lay_nice <- function(go, group = "module", mode = "circlepack", ...) {
  coors <- g_layout_nice(go, group, mode, ...)
  # 提取circlepack的分组坐标
  big_lay <- coors[grepl("group_", coors$name), ]
  big_lay$name <- gsub("group_", "", big_lay$name)
  return(big_lay)
}

#' @rdname g_layout_nice
#' @export
g_layout_circlepack <- \(go, group = "module", ...){
  g_layout_nice(go, group = group, mode = "circlepack", ...)
}

#' @rdname g_layout_nice
#' @export
g_layout_treemap <- \(go, group = "module", ...){
  g_layout_nice(go, group = group, mode = "treemap", ...)
}

#' @rdname g_layout_nice
#' @export
g_layout_backbone <- \(go, group = "module", ...){
  g_layout_nice(go, group = group, mode = "backbone", ...)
}

#' @rdname g_layout_nice
#' @export
g_layout_stress <- \(go, group = "module", ...){
  g_layout_nice(go, group = group, mode = "stress", ...)
}


#' Generate spatial layout using `spatstat`
#'
#' @param go igraph or metanet object
#' @param win A spatstat window object (owin), e.g. disc(), owin(poly=...); Or sf object.
#' @param type Type of distribution: "random", "regular"
#' @param mode "surface", "boundary"
#' @param jitter for surface-regular, defalut 0
#' @param curved Optional curved attribute for coors
#' @param order_by order nodes according to a node attribute
#' @param order_ls manual the discrete variable with a vector, or continuous variable with "desc" to decreasing
#' @param seed random seed
#' @param rescale rescale the coordinates to (0,1)
#' @param order_circle order nodes from the center of a circle
#'
#' @return A coors object (data.frame with class "coors" and attribute "curved")
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("spatstat.geom") && requireNamespace("spatstat.random")) {
#'   poly_x <- c(0, 2, 2, 0)
#'   poly_y <- c(0, 0, 1, 1)
#'   win_poly <- spatstat.geom::owin(poly = list(x = poly_x, y = poly_y))
#'   plot(win_poly)
#'   coors1 <- spatstat_layout(co_net, win_poly, type = "regular", mode = "surface")
#'   plot(co_net, coors = coors1)
#'   coors2 <- spatstat_layout(co_net2, win_poly, type = "random", mode = "boundary")
#'   plot(co_net2, coors = coors2)
#'
#'   if (requireNamespace("sf")) {
#'     nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'     poly <- nc[1, ]
#'     coors <- spatstat_layout(go = multi1, win = poly, type = "regular", mode = "surface")
#'     plot(multi1, coors = coors)
#'   }
#' }
#' }
spatstat_layout <- function(go, win,
                            type = c("random", "regular"),
                            mode = c("surface", "boundary"),
                            jitter = 0,
                            curved = NULL,
                            order_by = NULL, order_ls = NULL, order_circle = FALSE,
                            seed = 1234, rescale = TRUE) {
  lib_ps("spatstat.geom", library = FALSE)
  lib_ps("spatstat.random", library = FALSE)

  type <- match.arg(type)
  mode <- match.arg(mode)

  if (inherits(win, "sf")) {
    win <- sf_to_owin(win)
  }
  if (!inherits(win, "owin")) stop("Input 'win' must be a spatstat 'owin' object.")

  n_nodes <- length(go)
  set.seed(seed)
  pts <- switch(mode,
    surface = switch(type,
      random = generate_random_points(n_nodes, win),
      regular = generate_regular_points(n_nodes, win, jitter = jitter)
    ),
    boundary = switch(type,
      random = generate_boundary_points_random(n_nodes, win),
      regular = generate_boundary_points_regular(n_nodes, win)
    )
  )
  if (!is.null(attributes(pts)$rejects)) {
    pts$x <- c(pts$x, attributes(pts)$rejects$x)
    pts$y <- c(pts$y, attributes(pts)$rejects$y)
  }

  if (!order_circle) {
    coors <- data.frame(X = pts$x, Y = pts$y) %>%
      dplyr::arrange_all() %>%
      as.matrix()
  } else {
    coors <- order_circle_func(as.matrix(data.frame(X = pts$x, Y = pts$y)))
  }

  # order
  get_v(go) -> tmp_v
  coors <- order_tmp_v_name(tmp_v, coors, order_by = order_by, order_ls = order_ls)

  coors <- as_coors(coors, curved = curved)
  if (rescale) coors <- rescale_coors(coors)
  coors
}

generate_regular_points <- function(n_nodes, win, max_iter = 20, jitter = 0) {
  xlen <- diff(win$xrange)
  ylen <- diff(win$yrange)
  ratio <- xlen / ylen

  for (factor in seq(1.1, 5, length.out = max_iter)) {
    ny <- ceiling(sqrt(n_nodes / ratio) * factor)
    nx <- ceiling(ny * ratio)

    gx <- seq(win$xrange[1], win$xrange[2], length.out = nx)
    gy <- seq(win$yrange[1], win$yrange[2], length.out = ny)
    grid <- expand.grid(x = gx, y = gy)

    if (jitter > 0) {
      gx_step <- mean(diff(gx))
      gy_step <- mean(diff(gy))
      grid$x <- grid$x + runif(nrow(grid), -jitter * gx_step, jitter * gx_step)
      grid$y <- grid$y + runif(nrow(grid), -jitter * gy_step, jitter * gy_step)
    }

    inside <- spatstat.geom::inside.owin(grid$x, grid$y, win)
    valid_pts <- grid[inside, ]
    if (nrow(valid_pts) >= n_nodes) {
      selected <- valid_pts[sample(nrow(valid_pts), n_nodes), ]
      return(spatstat.geom::ppp(selected$x, selected$y, window = win))
    }
  }
  stop("Failed to generate enough regularly spaced points inside window after multiple attempts.")
}


generate_random_points <- function(n_nodes, win) {
  spatstat.random::runifpoint(n_nodes, win = win)
}

generate_boundary_points_random <- function(n_nodes, win) {
  bnd <- spatstat.geom::edges(win)
  pts <- spatstat.random::runifpointOnLines(n_nodes, bnd)
  spatstat.geom::ppp(pts$x, pts$y, window = win)
}

generate_boundary_points_regular <- function(n_nodes, win) {
  bnd <- spatstat.geom::edges(win)
  len <- spatstat.geom::lengths_psp(bnd)
  total_len <- sum(len)
  t_seq <- seq(0, total_len, length.out = n_nodes + 1)[-1]

  pts_x <- numeric(n_nodes)
  pts_y <- numeric(n_nodes)
  cum_len <- cumsum(len)

  for (i in seq_along(t_seq)) {
    t <- t_seq[i]
    seg_idx <- which(t <= cum_len)[1]
    t0 <- if (seg_idx == 1) 0 else cum_len[seg_idx - 1]
    frac <- (t - t0) / len[seg_idx]
    seg <- bnd[seg_idx]
    pts_x[i] <- seg$ends$x0 + frac * (seg$ends$x1 - seg$ends$x0)
    pts_y[i] <- seg$ends$y0 + frac * (seg$ends$y1 - seg$ends$y0)
  }

  spatstat.geom::ppp(pts_x, pts_y, window = win)
}

sf_to_owin <- function(sf_obj) {
  lib_ps("sf", library = FALSE)
  lib_ps("spatstat.geom", library = FALSE)

  # 确保输入是 sf 对象
  if (!inherits(sf_obj, "sf")) stop("Input must be an 'sf' object.")
  if (!sf::st_is(sf_obj, "POLYGON") && !sf::st_is(sf_obj, "MULTIPOLYGON")) {
    stop("sf object must be a POLYGON or MULTIPOLYGON")
  }

  # 修复方向问题
  sf_obj <- sf::st_make_valid(sf_obj) # 处理无效几何问题
  # 强制外环为逆时针，内环为顺时针
  sf_obj <- sf::st_set_geometry(sf_obj, sf::st_geometry(sf_obj))
  sf_obj <- sf::st_simplify(sf_obj) # 简化几何

  # 将 sf 对象转为 polylist 坐标
  poly_coords <- sf::st_coordinates(sf_obj)[, 1:2]
  polylist <- list(x = poly_coords[, 1], y = poly_coords[, 2])

  # 创建 spatstat 的 owin 对象
  win <- spatstat.geom::owin(poly = list(polylist))
  return(win)
}


create_sector_window <- function(r0 = 0, r1 = 1, theta_start = 0, theta_end = pi) {
  lib_ps("spatstat.geom", library = FALSE)
  # theta in radians
  n <- 100
  angles <- seq(theta_start, theta_end, length.out = n)
  x_outer <- r1 * cos(angles)
  y_outer <- r1 * sin(angles)
  x_inner <- r0 * cos(rev(angles))
  y_inner <- r0 * sin(rev(angles))
  x <- c(x_outer, x_inner)
  y <- c(y_outer, y_inner)
  spatstat.geom::owin(poly = list(x = x, y = y))
}



#' Layout with group
#'
#' @param go igraph
#' @param group group name (default:v_group)
#' @param group2 group2 name, will order nodes in each group according to group2_order
#' @param group2_order group2_order
#' @param group_order group_order
#' @param space the space between each arc, default: pi/4
#' @param scale_node_num scale with the node number in each group
#' @param type Type of distribution: "random", "regular"
#' @param mode "surface", "boundary"
#' @param jitter for surface-regular, defalut 0
#' @param curved Optional curved attribute for coors
#'
#' @return coors
#' @export
#' @family g_layout
g_layout_poly_sector <- function(go,
                                 group = "v_group", group_order = NULL,
                                 group2 = NULL, group2_order = NULL,
                                 space = pi / 4, jitter = 0,
                                 scale_node_num = TRUE,
                                 type = c("regular", "random"),
                                 mode = c("surface", "boundary"),
                                 curved = NULL) {
  type <- match.arg(type)
  mode <- match.arg(mode)

  get_v(go) -> tmp_v
  group1 <- as.factor(tmp_v[[group]])
  if (!is.null(group_order)) group1 <- pcutils::change_fac_lev(group1, group_order)

  n <- nlevels(group1)
  if (n < 2) stop("n should be greater than 1")

  g_num <- table(group1)
  sep <- space / n

  # 分配每组所占的弧度
  if (scale_node_num) {
    arc_r <- (2 * pi - space) * as.numeric(g_num) / sum(g_num)
  } else {
    arc_r <- rep((2 * pi - space) / n, n)
  }
  names(arc_r) <- levels(group1)

  coors_all <- list()
  theta1 <- 0
  for (i in names(arc_r)) {
    tmp_nodes <- tmp_v[tmp_v[[group]] == i, ]
    sub_go <- igraph::make_empty_graph(n = nrow(tmp_nodes)) %>%
      igraph::set_vertex_attr("name", value = tmp_nodes$name)

    win <- create_sector_window(r0 = 0, r1 = 1, theta_start = theta1, theta_end = theta1 + arc_r[i])
    tmp_layout <- spatstat_layout(sub_go, win = win, type = type, mode = mode, curved = curved, rescale = FALSE, jitter = jitter)

    tmp_layout <- order_tmp_v_name(tmp_nodes, tmp_layout[, c("X", "Y")], group2, group2_order)
    coors_all[[i]] <- tmp_layout
    theta1 <- theta1 + arc_r[i] + sep
  }

  # 合并所有 layout
  coors_df <- do.call(rbind, coors_all)
  as_coors(coors_df, curved = curved) %>% rescale_coors()
}
