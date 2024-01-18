# ========3.layout========
#' Layout coordinates
#'
#' @param go igraph or metanet
#' @param method as_star(), as_tree(), in_circle(), nicely(), on_grid(), on_sphere(),randomly(), with_dh(), with_fr(), with_gem(), with_graphopt(), with_kk(),with_lgl(), with_mds(),as_line(), as_arc(), as_polygon(), as_polyarc(). see \code{\link[igraph]{layout_}}.
#' @param order_by order nodes according to a node attribute
#' @param order_ls manual the discrete variable with a vector, or continuous variable with "desc" to decreasing
#' @param seed random seed
#' @param line_curved consider line curved, only for some layout methods like as_line(), as_polygon().default:0
#'
#' @return coordinates for nodes,columns: name, X, Y
#' @export
#' @examples
#' library(igraph)
#' c_net_lay(co_net) -> coors
#' c_net_plot(co_net, coors)
#' c_net_plot(co_net, c_net_lay(co_net, in_circle()), vertex.size = 2)
#' c_net_plot(co_net, c_net_lay(co_net, in_circle(), order_by = "v_class"), vertex.size = 2)
#' c_net_plot(co_net, c_net_lay(co_net, in_circle(), order_by = "size", order_ls = "desc"))
#' c_net_plot(co_net, c_net_lay(co_net, as_polygon(3)))
c_net_lay <- function(go, method = igraph::nicely(), order_by = NULL, order_ls = NULL, seed = 1234, line_curved = 1) {
    set.seed(seed)

    if ("igraph_layout_spec" %in% class(method)) {
        coors <- igraph::layout_(go, method)
    } else if ("poly" %in% class(method)) {
        coors <- method(go, group2 = order_by, group2_order = order_ls)
    } else if ("layout" %in% class(method)) {
        coors <- method(go)
    } else {
        stop("No valid method")
    }

    # order
    if (is.matrix(coors)) {
        if (is.null(order_by)) {
            coors <- data.frame(name = V(go)$name, X = coors[, 1], Y = coors[, 2], row.names = NULL)
        } else {
            get_v(go) -> tmp_v
            ordervec <- tmp_v[, order_by]
            if (is.numeric(ordervec)) {
                name <- tmp_v[order(ordervec, decreasing = is.null(order_ls)), "name"]
            } else {
                ordervec <- pcutils::change_fac_lev(ordervec, order_ls)
                name <- tmp_v[order(ordervec), "name"]
            }
            coors <- data.frame(name = name, X = coors[, 1], Y = coors[, 2], row.names = NULL)
        }
    }

    # if line type, need to consider edge.curved

    if (line_curved) {
        if ("line" %in% class(method)) {
            tmp_e <- data.frame(igraph::as_data_frame(go))[, c("from", "to")]
            if (nrow(tmp_e) > 0) {
                curved <- data.frame(tmp_e, curved = line_curved, row.names = NULL)
            } else {
                curved <- NULL
            }
            coors <- list(coors = data.frame(coors, row.names = NULL), curved = curved)
            class(coors) <- c("coors", class(coors))
        }
    }
    return(coors)
}

is_lay <- \(x){
    any(class(x) %in% c("igraph_layout_spec", "layout"))
}

get_coors <- \(coors, go, ...){
    edge_curved <- NULL
    # 1.如果coors是NULL，去graph_attr找一下，没有的话就默认nicely计算
    if (is.null(coors)) {
        if (is.null(igraph::graph_attr(go, "coors"))) {
            coors <- c_net_lay(go, igraph::nicely(), ...)
        } else {
            coors <- igraph::graph_attr(go, "coors")
        }
    }
    # 2.如果是layout函数，那就计算
    if (is_lay(coors)) coors <- c_net_lay(go, coors, ...)
    # 3.如果是一个提供的coors对象（list），那就把curved传到edge里，coors导出为df
    if (inherits(coors, "coors")) {
        if (!is.null(coors$curved)) {
            edge_curved <- dplyr::left_join(get_e(go), coors$curved, by = c("from", "to"), suffix = c(".x", "")) %>%
                dplyr::select("from", "to", "curved")
        }
        if (!is.null(coors$coors)) {
            coors <- coors$coors
        } else {
            coors <- c_net_lay(go, igraph::nicely(), ...)
        }
    }
    # 4.如果是matrix，变成df
    if (is.data.frame(coors)) {
        if (!"name" %in% colnames(coors)) coors <- as.matrix(coors)
    }
    if (is.matrix(coors)) coors <- data.frame(name = V(go)$name, X = coors[, 1], Y = coors[, 2])
    # 5.如果是df了，那就对齐name用于下一步的绘图，
    if (is.data.frame(coors)) {
        coors <- coors[match(V(go)$name, coors$name), ]
        return(list(coors = coors, curved = edge_curved))
    }
    stop("coors wrong!")
}


#' Layout as a line
#'
#' @param angle anticlockwise rotation angle
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
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


#' Layout with group
#'
#' @param go igraph or metanet object
#' @param group group name (default: module)
#' @param zoom1 big network layout size
#' @param zoom2 average sub_network layout size, or numeric vector, or "auto"
#' @param layout1 layout1 method, one of
#'      (1) a dataframe or matrix: rowname is group, two columns are X and Y
#'      (2) function: layout method for \code{\link{c_net_lay}} default: in_circle()
#' @param layout2 one of functions: layout method for \code{\link{c_net_lay}}, or a list of functions.
#' @param show_big_lay show the big layout to help you adjust.
#' @param ... add
#'
#' @return coors
#' @export
#'
#' @examples
#' \donttest{
#' data("c_net")
#' modu_dect(co_net, method = "cluster_fast_greedy") -> co_net_modu
#' g_lay(co_net_modu, group = "module", zoom1 = 30, zoom2 = 1:5, layout2 = as_line()) -> oridata
#' plot(co_net_modu, coors = oridata)
#' g_lay_nice(co_net_modu, group = "module") -> oridata
#' plot(co_net_modu, coors = oridata)
#' }
g_lay <- function(go, group = "module", layout1 = in_circle(), zoom1 = 20, layout2 = in_circle(),
                  zoom2 = 3, show_big_lay = FALSE, ...) {
    name <- ID <- NULL

    stopifnot(is_igraph(go))
    if (!group %in% igraph::vertex_attr_names(go)) stop("no group named ", group, " !")
    get_v(go) %>% dplyr::select(name, !!group) -> nodeGroup
    colnames(nodeGroup) <- c("ID", "group")
    nodeGroup$group <- as.factor(nodeGroup$group)

    big_lay <- \(zoom1, layout1){
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
    da <- big_lay(zoom1, layout1)
    if (show_big_lay) {
        message("Big layout:")
        print(da)
        graphics::plot.new()
        plot(da, pch = 21, bty = "n", bg = get_cols(nrow(da), "col2"), main = "Big layout coordinates")
        graphics::text(da$X * 0.8, da$Y * 0.9, rownames(da))
        return(invisible())
    }

    # layout of vertexes in one group
    layoutls <- list()
    if (is_lay(layout2)) {
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

    # get coors
    oridata <- list()
    oridata2 <- list()
    for (i in levels(nodeGroup$group)) {
        nodeGroup %>%
            dplyr::filter(group == i) %>%
            dplyr::pull(ID) -> tmpid
        igraph::subgraph(go, tmpid) -> tmp_net

        if (TRUE) {
            get_coors(layoutls[[i]], tmp_net, ...) -> coors
            data <- coors$coors
            if ("igraph_layout_spec" %in% class(layoutls[[i]])) {
                data[, c("X", "Y")] <- igraph::norm_coords(as.matrix(data[, c("X", "Y")]))
            }
            data[, "X"] <- data[, "X"] * zoom2[i] + da[i, "X"]
            data[, "Y"] <- data[, "Y"] * zoom2[i] + da[i, "Y"]
        }
        oridata[[i]] <- data
        oridata2[[i]] <- coors$curved
    }
    oridata <- do.call(rbind, oridata)
    oridata2 <- do.call(rbind, oridata2)

    coors <- list(
        coors = data.frame(oridata, row.names = NULL),
        curved = data.frame(oridata2, row.names = NULL)
    )
    if (nrow(coors$curved) == 0) coors$curved <- NULL
    class(coors) <- "coors"
    return(coors)
}


#' Layout with group as a polygon
#'
#' @param go igraph
#' @param group group name (default:v_group)
#' @param group2 group2 name, will order nodes in each group according to group2_order (default:v_class)
#' @param group2_order group2_order
#' @param line_curved line_curved 0~1
#'
#' @return coors
#' @export
#'
#' @examples
#' g_lay_polygon(multi1) -> oridata
#' c_net_plot(multi1, oridata)
#' g_lay_polyarc(multi1, group2_order = c(LETTERS[4:1])) -> oridata
#' c_net_plot(multi1, oridata)
g_lay_polygon <- function(go, group = "v_group", group2 = NULL, group2_order = NULL, line_curved = 0.5) {
    n <- length(unique(igraph::vertex.attributes(go)[[group]]))

    if (n < 2) stop("n should bigger than 1")
    angle_ls <- -pi / 2 + (seq(0, n - 1, 1)) * 2 * pi / n
    fun_ls <- lapply(angle_ls, \(i)as_line(i))

    g_lay(go,
        group = group, zoom1 = 1, zoom2 = 0.9 * (ifelse(n > 2, tan(pi / n), 2)),
        layout2 = fun_ls, order_by = group2, order_ls = group2_order
    ) -> oridata

    if (is.data.frame(oridata$curved)) oridata$curved$curved <- line_curved
    oridata
}

g_lay_polygon2 <- function(go, group = "v_group", group2 = NULL, group2_order = NULL, line_curved = 0.5) {
    n <- length(unique(igraph::vertex.attributes(go)[[group]]))

    if (n < 2) stop("n should bigger than 1")
    angle_ls <- -pi / 2 + (seq(0, n - 1, 1)) * 2 * pi / n
    fun_ls <- lapply(angle_ls, \(i)as_line(i))

    g_lay(go,
        group = group, zoom1 = 1, zoom2 = 0.9 * (ifelse(n > 2, tan(pi / n), 2)),
        layout2 = fun_ls, order_by = group2, order_ls = group2_order
    ) -> oridata

    if (is.data.frame(oridata$curved)) oridata$curved$curved <- line_curved
    oridata
}

#' Layout with group as a polyarc
#'
#' @param space the space between each arc, default: pi/4
#'
#' @rdname g_lay_polygon
#' @export
g_lay_polyarc <- function(go, group = "v_group", group2 = NULL, group2_order = NULL, space = pi / 4) {
    get_v(go) -> tmp_v
    group1 <- as.factor(tmp_v[, group])
    n <- nlevels(group1)
    if (n < 2) stop("n should bigger than 1")

    # consider each group numbers!!!
    g_num <- table(group1)
    sep <- space / n
    arc_r <- (2 * pi - space) * g_num / length(group1)

    # coordinate
    coors <- data.frame()
    theta1 <- 0
    for (i in names(arc_r)) {
        tmp_t <- seq(theta1, theta1 + arc_r[i], len = g_num[i])
        tmp_v1 <- tmp_v[tmp_v[, group] == i, ]
        {
            if (!is.null(group2)) {
                ordervec <- tmp_v1[, group2]
                if (is.numeric(ordervec)) {
                    name <- tmp_v1[order(ordervec, decreasing = is.null(group2_order)), "name"]
                } else {
                    ordervec <- pcutils::change_fac_lev(ordervec, group2_order)
                    name <- tmp_v1[order(ordervec), "name"]
                }
            } else {
                name <- tmp_v1$name
            }
        }
        tmp_coor <- data.frame(name = name, X = cos(tmp_t), Y = sin(tmp_t))
        coors <- rbind(coors, tmp_coor)
        theta1 <- theta1 + arc_r[i] + sep
    }
    coors
}


#' Layout as a polygon
#'
#' @param n how many edges of this polygon
#' @param line_curved line_curved 0~0.5
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#' @examples
#' as_polygon()(co_net)
as_polygon <- function(n = 3, line_curved = 0.5) {
    fun <- \(go, group2 = NULL, group2_order = NULL){
        V(go)$poly_group <- rep(paste0("omic", seq_len(n)), len = length(go))
        if (n < 2) stop("n should bigger than 1")
        g_lay_polygon(go,
            group = "poly_group", group2 = group2, group2_order = group2_order,
            line_curved = line_curved
        ) -> oridata
        oridata
    }
    class(fun) <- c("poly", "layout", "function")
    fun
}

#' Layout as a poly_arc
#'
#' @param n how many arcs of this poly_arc
#' @param space the space between each arc, default: pi/3
#'
#' @return A two-column matrix, each row giving the coordinates of a vertex, according to the ids of the vertex ids.
#' @export
#'
#' @examples
#' as_polyarc()(co_net)
as_polyarc <- \(n = 3, space = pi / 3){
    fun <- \(go, group2 = NULL, group2_order = NULL){
        V(go)$poly_group <- rep(paste0("omic", 1:n), len = length(go))
        if (n < 2) stop("n should bigger than 1")
        g_lay_polyarc(go, "poly_group", space = space, group2 = group2, group2_order = group2_order) -> oridata
        oridata
    }
    class(fun) <- c("poly", "layout", "function")
    fun
}


#' Layout with group nicely
#'
#' @param go igraph or metanet
#' @param group group name (default: module)
#'
#' @export
#'
#' @rdname g_lay
g_lay_nice <- function(go, group = "module") {
    name <- size <- leaf <- x <- y <- NULL
    lib_ps("ggraph", library = FALSE)
    stopifnot(is_igraph(go))
    if (!group %in% vertex_attr_names(go)) stop("no group named ", group, " !")
    get_v(go) %>% dplyr::select(name, !!group) -> nodeGroup
    colnames(nodeGroup) <- c("ID", "group")
    nodeGroup$group <- as.factor(nodeGroup$group)

    edge <- data.frame(from = paste("group_", nodeGroup$group, sep = ""), to = nodeGroup$ID)

    vertices_t <- data.frame(name = unique(c(
        as.character(edge$from),
        as.character(edge$to)
    )))
    vertices_t$size <- sample(1:10, nrow(vertices_t), replace = TRUE)

    mygraph <- igraph::graph_from_data_frame(edge, vertices = vertices_t)
    data <- ggraph::create_layout(mygraph, layout = "circlepack", weight = size)
    coor <- data %>%
        dplyr::filter(leaf == TRUE) %>%
        dplyr::select(name, x, y)
    colnames(coor) <- c("name", "X", "Y")
    return(coor)
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


#' Plot a metanet
#'
#' @param go an igraph or metanet object
#' @param coors the coordinates you saved
#' @param ... additional parameters for \code{\link[igraph]{igraph.plotting}}
#' @param labels_num show how many labels,>1 indicates number, <1 indicates fraction, "all" indicates all, default:5
#' @param vertex_size_range the vertex size range, e.g. c(1,10)
#' @param legend_number legend with numbers
#' @param legend all legends
#' @param legend_position legend_position, default: c(left_leg_x=-1.9,left_leg_y=1,right_leg_x=1.2,right_leg_y=1)
#' @param legend_cex 	character expansion factor relative to current par("cex"), default: 1
#' @param lty_legend logical
#' @param size_legend logical
#' @param edge_legend logical
#' @param color_legend logical
#' @param width_legend logical
#' @param lty_legend_title lty_legend_title
#' @param size_legend_title size_legend_title
#' @param edge_legend_title edge_legend_title
#' @param edge_legend_order edge_legend_order vector, e.g. c("positive","negative")
#' @param width_legend_title width_legend_title
#' @param color_legend_order color_legend_order vector,
#' @param group_legend_title group_legend_title, length must same to the numbers of v_group
#' @param group_legend_order group_legend_order vector
#' @param mark_module logical, mark the modules?
#' @param mark_color mark colors
#' @param mark_alpha mark fill alpha, default 0.3
#' @param seed random seed, default:1234, make sure each plot is the same.
#' @param plot_module logical, plot module?
#' @param module_label module_label
#' @param module_label_cex module_label_cex
#' @param module_label_color module_label_color
#' @param module_label_just module_label_just
#'
#' @return a network plot
#' @export
#'
#' @examples
#' data("c_net")
#' c_net_plot(co_net)
#' c_net_plot(co_net2)
#' c_net_plot(multi1)
c_net_plot <- function(go, coors = NULL, ..., labels_num = 5, vertex_size_range = NULL,
                       legend_number = FALSE, legend = TRUE, legend_cex = 1,
                       legend_position = c(left_leg_x = -2, left_leg_y = 1, right_leg_x = 1.2, right_leg_y = 1),
                       lty_legend = FALSE, lty_legend_title = "Edge class",
                       size_legend = FALSE, size_legend_title = "Node Size",
                       edge_legend = TRUE, edge_legend_title = "Edge type", edge_legend_order = NULL,
                       width_legend = FALSE, width_legend_title = "Edge width",
                       color_legend = TRUE, color_legend_order = NULL,
                       group_legend_title = NULL, group_legend_order = NULL,
                       plot_module = FALSE, mark_module = FALSE, mark_color = NULL, mark_alpha = 0.3,
                       module_label=TRUE, module_label_cex=2, module_label_color="black", module_label_just=c(0,0),
                       seed = 1234) {
    name <- size <- color <- e_type <- lty <- e_class <- v_class <- shape <- NULL
    lib_ps("igraph", library = FALSE)
    set.seed(seed)

    # modules
    if (plot_module) {
        go <- to_module_net(go)
        if (is.null(group_legend_title)) group_legend_title <- "Module"
    }

    # get network type
    main <- "Network"
    if (!is.null(get_n(go)$n_type)) {
        switch(get_n(go)$n_type,
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
                main <- paste0(get_n(go)$n_modules, "-modules network")
            },
            "skeleton" = {
                main <- paste0(get_n(go)$skeleton, " skeleton network")
            },
            "venn" = {
                main <- "Venn network"
            },
            default = {
                main <- "Network"
            }
        )
    }

    get_v(go) -> tmp_v
    get_e(go) -> tmp_e

    # get coordinates
    ori_coors <- get_coors(coors, go, seed = seed)
    if (is.null(ori_coors$curved)) {
        edge_curved <- 0
    } else {
        edge_curved <- ori_coors$curved$curved
    }
    edge_curved[is.na(edge_curved)] <- 0

    coors <- ori_coors$coors[, c("X", "Y")] %>% as.matrix()

    # scale the size and width
    {
        v_groups <- unique(tmp_v$v_group)
        nice_size <- ceiling(60 / sqrt(length(V(go)))) + 1
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

        node_size_text1 <- c()
        node_size_text2 <- c()
        for (i in v_groups) {
            node_size_text1 <- c(node_size_text1, min(tmp_v[tmp_v$v_group == i, "size"]))
            node_size_text2 <- c(node_size_text2, max(tmp_v[tmp_v$v_group == i, "size"]))
            tmp_v[tmp_v$v_group == i, "size"] <- do.call(pcutils::mmscale, append(
                list(tmp_v[tmp_v$v_group == i, "size"]),
                as.list(vertex_size_range[[i]][1:2])
            ))
        }
        node_size_text1 <- round(node_size_text1, 3)
        node_size_text2 <- round(node_size_text2, 3)

        tmp_e$width <- pcutils::mmscale(tmp_e$width, 0.5, 1)
    }
    # some custom parameters
    params <- list(...)
    params_name <- names(params)
    if ("vertex.size" %in% params_name) tmp_v$size <- params[["vertex.size"]]
    if ("vertex.color" %in% params_name) tmp_v$color <- condance(data.frame(tmp_v$color, pcutils::tidai(tmp_v$v_class, params[["vertex.color"]])))
    if ("vertex.shape" %in% params_name) tmp_v$shape <- condance(data.frame(tmp_v$shape, pcutils::tidai(tmp_v$v_group, params[["vertex.shape"]])))
    if ("vertex.label" %in% params_name) tmp_v$label <- params[["vertex.label"]]
    if ("vertex.label.color" %in% params_name) {
        vertex.label.color <- condance(data.frame("black", pcutils::tidai(tmp_v$v_group, params[["vertex.label.color"]])))
    } else {
        vertex.label.color <- "black"
    }

    if ("edge.color" %in% params_name) tmp_e$color <- condance(data.frame(tmp_e$color, pcutils::tidai(tmp_e$e_type, params[["edge.color"]])))
    if ("edge.lty" %in% params_name) tmp_e$lty <- condance(data.frame(tmp_e$lty, pcutils::tidai(tmp_e$e_class, params[["edge.lty"]])))
    if ("edge.width" %in% params_name) tmp_e$width <- params[["edge.width"]]

    # show labels
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
    # modules set
    {
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
        if (module_label){
            rescale_flag=TRUE
            module_coors=dplyr::left_join(tmp_v[,c("name","module")],ori_coors$coors,by = "name")
            if("rescale"%in%names(params)){
                if(!params[["rescale"]])rescale_flag=FALSE
            }
            if(rescale_flag) module_coors=dplyr::mutate(module_coors,X=mmscale(X,-1,1),Y=mmscale(Y,-1,1))
            module_coors=dplyr::group_by(module_coors,module)%>%
                dplyr::summarise(minx=min(X),maxx=max(X),miny=min(Y),maxy=max(Y))
            module_label_just=rep(module_label_just,2)
            module_coors=mutate(module_coors,
                                X=minx+module_label_just[1]*(maxx-minx),
                                Y=miny+module_label_just[2]*(maxy-miny))
        }
    }
    # main plot
    {
        old_xpd <- graphics::par(mar = c(4, 2, 2, 2), xpd = TRUE)
        on.exit(graphics::par(old_xpd), add = TRUE)
        # oldpar <- graphics::par(no.readonly = TRUE)
        # on.exit(graphics::par(oldpar),add = TRUE)
        igraph::plot.igraph(go,
            layout = coors,
            vertex.size = tmp_v$size,
            vertex.color = tmp_v$color,
            vertex.shape = tmp_v$shape,
            edge.color = tmp_e$color,
            edge.lty = tmp_e$lty,
            edge.width = tmp_e$width,
            mark.groups = new_modu,
            mark.col = pcutils::add_alpha(module_color[names(new_modu)], mark_alpha),
            mark.border = module_color[names(new_modu)],
            vertex.label.color = vertex.label.color,
            ...,
            main = main,
            vertex.label.font = 1,
            vertex.label.cex = 0.08 * tmp_v$size,
            vertex.label = tmp_v$label,
            edge.arrow.size = 0.3,
            edge.arrow.width = 0.5,
            edge.curved = edge_curved,
            margin = c(0, 0, 0, 0)
        )
    }

    # add module_label
    if (module_label){
        n_module=nrow(module_coors)
        module_label_cex=rep(module_label_cex,n_module)
        module_label_color=rep(module_label_color,n_module)
        for (i in seq_len(n_module)) {
            text(module_coors[i,"X"],module_coors[i,"Y"],labels = module_coors[i,"module"],
                 cex=module_label_cex,col=module_label_color[i])
        }
    }

    if (!legend) {
        return(invisible())
    }

    # produce legends
    legend_position_default <- c(left_leg_x = -2, left_leg_y = 1, right_leg_x = 1.2, right_leg_y = 1)
    if (is.null(legend_position)) legend_position <- legend_position_default
    if (is.null(names(legend_position))) {
        legend_position <- setNames(
            legend_position,
            names(legend_position_default)[seq_along(legend_position)]
        )
    }
    legend_position <- pcutils::update_param(legend_position_default, legend_position)
    left_leg_x <- legend_position["left_leg_x"]
    left_leg_y <- legend_position["left_leg_y"]
    right_leg_x <- legend_position["right_leg_x"]
    right_leg_y <- legend_position["right_leg_y"]

    if (size_legend) {
        legend(right_leg_x, right_leg_y,
            cex = 0.7 * legend_cex, title.font = 2, title = size_legend_title, title.adj = 0,
            legend = c(paste(node_size_text1, collapse = "/"), paste(node_size_text2, collapse = "/")), adj = 0, title.cex = 0.8 * legend_cex,
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
            legend = c(min(E(go)$width), max(E(go)$width)) %>% round(., 3), adj = 0, title.cex = 0.8 * legend_cex,
            col = "black", bty = "n", lty = 1, lwd = c(min(tmp_e$width), max(tmp_e$width))
        )
        right_leg_y <- right_leg_y - 0.5 * legend_cex
    }

    if (lty_legend) {
        edges <- levels(factor(tmp_e$e_class))
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

    if (color_legend) {
        pchls <- c("circle" = 21, "square" = 22)
        vgroups <- pcutils::change_fac_lev(tmp_v$v_group, group_legend_order)
        vgroups <- unique(vgroups)
        if (is.null(group_legend_title)) {
            group_legend_title <- setNames(vgroups, vgroups)
        } else if (is.null(names(group_legend_title))) group_legend_title <- setNames(rep(group_legend_title, len = length(vgroups)), vgroups)

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
}

#' Transfer an igraph object to a ggig
#'
#' @param go igraph or meatnet
#' @param coors coordinates for nodes,columns: name, X, Y
#'
#' @return ggig object
#' @export
#'
#' @examples
#' to.ggig(co_net, coors = c_net_lay(co_net)) -> ggig
#' plot(ggig)
#' to.ggig(multi1, coors = c_net_lay(multi1)) -> ggig
#' plot(ggig, labels_num = 0.3)
to.ggig <- function(go, coors = NULL) {
    list(n_index = get_n(go), v_index = get_v(go), e_index = get_e(go)) -> net_par_res

    if (is.null(coors)) coors <- c_net_lay(go)
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
#' @param coors the coordinates you saved
#' @param ... additional parameters
#' @param labels_num show how many labels,>1 indicates number, <1 indicates fraction ,default:5
#' @param legend_number legend with numbers
#' @param legend all legends
#' @param lty_legend logical
#' @param size_legend logical
#' @param edge_legend logical
#' @param color_legend logical
#' @param width_legend logical
#' @param lty_legend_title lty_legend_title
#' @param size_legend_title size_legend_title
#' @param edge_legend_title edge_legend_title
#' @param edge_legend_order edge_legend_order vector, e.g. c("positive","negative")
#' @param width_legend_title width_legend_title
#' @param color_legend_order color_legend_order vector,
#' @param group_legend_title group_legend_title, length must same to the numbers of v_group
#' @param group_legend_order group_legend_order vector
#'
#' @return ggplot
#' @exportS3Method
#' @method plot ggig
plot.ggig <- function(x, coors = NULL, ..., labels_num = 0,
                      legend_number = FALSE, legend = TRUE,
                      lty_legend = FALSE, lty_legend_title = "Edge class",
                      size_legend = FALSE, size_legend_title = "Node Size",
                      edge_legend = TRUE, edge_legend_title = "Edge type", edge_legend_order = NULL,
                      width_legend = FALSE, width_legend_title = "Edge width",
                      color_legend = TRUE, color_legend_order = NULL,
                      group_legend_title = "Node class", group_legend_order = NULL) {
    rename <- size <- name <- color <- e_type <- lty <- e_class <- v_class <- shape <- X1 <- Y1 <- X2 <- Y2 <- width <- X <- Y <- label <- NULL
    ggig <- x
    lib_ps("ggplot2", "ggnewscale", library = FALSE)
    ggig$v_index -> tmp_v
    ggig$e_index -> tmp_e
    set.seed(1234)
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
    main <- "Network"
    if (!is.null(ggig$n_index$n_type)) {
        switch(ggig$n_index$n_type,
            "single" = {
                main <- "Correlation network"
            },
            "bipartite" = {
                main <- "Bipartite network"
            },
            "multi_full" = {
                main <- "Multi-omics network"
            },
            "skeleton" = {
                main <- paste0(ggig$n_index$skeleton, " skeleton network")
            },
            default = {
                main <- "Network"
            }
        )
    }
    # scale the size and width
    {
        node_size_text1 <- c()
        node_size_text2 <- c()
        for (i in unique(tmp_v$v_group)) {
            node_size_text1 <- c(node_size_text1, min(tmp_v[tmp_v$v_group == i, "size"]))
            node_size_text2 <- c(node_size_text2, max(tmp_v[tmp_v$v_group == i, "size"]))
            tmp_v[tmp_v$v_group == i, "size"] %<>% mmscale(., 2, 10)
        }
        node_size_text <- c(paste(node_size_text1, collapse = "/"), paste(node_size_text2, collapse = "/"))
        edge_width_text <- c(min(tmp_e$width), max(tmp_e$width))
        tmp_e$width <- pcutils::mmscale(tmp_e$width, 0.5, 1)
        # new shapes
        tmp_v$shape <- tidai(tmp_v$v_group, 21:26)

        # some custom parameters
        params <- list(...)
        params_name <- names(params)
        if ("vertex.size" %in% params_name) tmp_v$size <- params[["vertex.size"]]
        if ("vertex.color" %in% params_name) tmp_v$color <- condance(data.frame(tmp_v$color, tidai(tmp_v$v_class, params[["vertex.color"]])))
        if ("vertex.shape" %in% params_name) tmp_v$shape <- condance(data.frame(tmp_v$shape, tidai(tmp_v$v_group, params[["vertex.shape"]])))
        if ("edge.color" %in% params_name) tmp_e$color <- condance(data.frame(tmp_e$color, tidai(tmp_e$e_type, params[["edge.color"]])))
        if ("edge.lty" %in% params_name) tmp_e$lty <- condance(data.frame(tmp_e$lty, tidai(tmp_e$e_class, params[["edge.lty"]])))
        if ("edge.width" %in% params_name) tmp_e$width <- params[["edge.width"]]
        if ("vertex.label" %in% params_name) tmp_e$label <- params[["vertex.label"]]
    }
    # show labels
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
    tmp_v$label <- ifelse(tmp_v$name %in% toplabel, tmp_v$label, NA)

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
#'
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
#'
#' @examples
#' plot(co_net2)
#' netD3plot(co_net2)
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
#'
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
#'
#' @examples
#' twocol <- data.frame(
#'     "col1" = sample(letters, 30, replace = TRUE),
#'     "col2" = sample(c("A", "B"), 30, replace = TRUE)
#' )
#' twocol_net <- twocol_edgelist(twocol)
#' plot(twocol_net)
#' c_net_plot(twocol_net, g_lay_polygon(twocol_net))
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


#' My Network plot
#'
#' @param test a dataframe with hierarchical structure
#' @param vertex_anno vertex annotation table
#' @param vertex_group vertex_group
#' @param vertex_class vertex_class
#' @param vertex_size vertex_size
#' @param ... look for parameters in \code{\link[MetaNet]{c_net_plot}}
#'
#' @return metanet
#' @export
#'
#' @examples
#' \donttest{
#' data(otutab, package = "pcutils")
#' cbind(taxonomy, num = rowSums(otutab))[1:10, ] -> test
#' my_network_tree(test)
#' }
my_network_tree <- function(test, vertex_anno = NULL,
                            vertex_group = "v_group", vertex_class = "v_class",
                            vertex_size = "size", ...) {
    lib_ps("MetaNet", library = FALSE)
    test <- as.data.frame(test)
    nc <- ncol(test)
    if (nc < 3) stop("as least 3-columns dataframe")
    if (!is.numeric(test[, nc])) stop("the last column must be numeric")
    ttt <- MetaNet::df2net_tree(test)
    ttt <- MetaNet::c_net_set(ttt, vertex_anno,
        vertex_group = vertex_group,
        vertex_class = vertex_class, vertex_size = vertex_size
    )
    message("For more details for network visualization, please refer to MetaNet ('https://github.com/Asa12138/MetaNet').")
    plot(ttt, ...)
}


#' Plot olympic rings using network
#'
#' @return network plot
#' @export
#'
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
    g1 <- modu_net(module_number = 5, n_node_in_module = 30)
    plot(g1,
        coors = g_lay(g1, layout1 = rings_data[, 1:2], zoom1 = 1.2, zoom2 = 0.5),
        rescale = FALSE, legend = FALSE, main = "Olympic Rings", vertex.frame.color = NA,
        edge.width = 0, vertex.color = setNames(rings_data$color, 1:5), vertex.size = 9
    )
}
