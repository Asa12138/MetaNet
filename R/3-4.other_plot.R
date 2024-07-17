# ========4.1.other_plot========
#
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
plot.ggig <- function(x, coors = NULL, ..., labels_num = NULL,
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
                      params_list = NULL,
                      seed = 1234) {
  if (!is.null(params_list)) {
    as.list(match.call()[-1]) -> set_params_list
    set_params_list[["params_list"]] <- NULL
    for (i in seq_along(params_list)) {
      if (names(params_list)[i] %in% names(set_params_list)) {
        message("The parameter `", names(params_list)[i], "` is duplicated, the format argument will be used.")
      }
    }
    pcutils::update_param(params_list, set_params_list) -> set_params_list
    do.call(c_net_plot, set_params_list)
    return(invisible())
  }

  rename <- size <- color <- e_type <- lty <- e_class <- v_class <- shape <- X1 <- Y1 <- X2 <- Y2 <- width <- X <- Y <- label <- NULL
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
    # scale_shape_manual(values =setNames(default_v_shape[node_shapes],vclass))+#node shape
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

  gephi <- c_net_update(gephi, initialize = TRUE)
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
  venn_net <- c_net_from_edgelist(edgelist, vertex_df = nodelist)
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
  venn_net <- c_net_from_edgelist(edgelist, vertex_df = nodelist)
  graph.attributes(venn_net)$n_type <- "twocol"
  # venn_net=c_net_set(venn_net,edge_type = "from")
  venn_net
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
#' if (requireNamespace("ggraph")) plot(ttt, coors = as_circle_tree())
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
  net <- c_net_set(net, vertex_class = "level", vertex_size = name, edge_width = name)
  net <- c_net_update(net, initialize = TRUE, verbose = FALSE)
  graph_attr(net, "coors") <- c_net_layout(net, as_tree())
  net
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
