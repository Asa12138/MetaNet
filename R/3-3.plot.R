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
    nice_size <- ceiling(100 / sqrt(nrow(tmp_v))) + 1

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
        as.list(sort(vertex_size_range[[i]][1:2]))
      ))
    }    }

  {
    edge_width_range_default <- vertex_size_range_default[[1]] / 6
    if (is.null(edge_width_range)) edge_width_range <- edge_width_range_default
    edge_width_range <- sort(edge_width_range)
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

  # 下列都是mapping好的颜色，形状等，无法单独修改某一个vertex或edge的颜色，形状等
  # ！！！考虑增加额外参数，用于单独修改某一个vertex或edge的颜色，形状等
  tmp_v$label.color <- "black"
  if ("vertex.size" %in% params_name) tmp_v$size <- params[["vertex.size"]]
  if ("vertex.color" %in% params_name) {
    tmp_v$color <- condance(data.frame(
      tmp_v$color,
      pcutils::tidai(tmp_v$v_class, params[["vertex.color"]], fac = TRUE)
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
    tmp_v$label.color <- condance(data.frame(
      "black",
      pcutils::tidai(tmp_v$v_group, params[["vertex.label.color"]], fac = TRUE)
    ))
  }

  if ("edge.color" %in% params_name) {
    tmp_e$color <- condance(data.frame(
      tmp_e$color,
      pcutils::tidai(tmp_e$e_type, params[["edge.color"]], fac = TRUE)
    ))
  }
  if ("edge.lty" %in% params_name) {
    tmp_e$lty <- condance(data.frame(
      tmp_e$lty,
      pcutils::tidai(tmp_e$e_class, params[["edge.lty"]], fac = TRUE)
    ))
  }
  if ("edge.width" %in% params_name) tmp_e$width <- params[["edge.width"]]

  envir <- parent.frame()
  assign("tmp_e", tmp_e, envir)
  assign("tmp_v", tmp_v, envir)
}

pie_set_for_plot <- function(tmp_v, pie_value, pie_color) {
  if ("pie" %in% tmp_v$shape) {
    if (!is.null(pie_value)) {
      if (!is.data.frame(pie_value)) stop("pie_value must be a data.frame.")
      if (!"name" %in% colnames(pie_value)) {
        pie_value$name <- rownames(pie_value)
      }
      if (any(duplicated(pie_value$name))) {
        stop(
          "Duplicated name in annotation tables: ",
          paste0(pie_value$name[duplicated(pie_value$name)], collapse = ", ")
        )
      }
      tmp_merge_df <- left_join(tmp_v["name"], pie_value, by = "name")
      pie_value <- tmp_merge_df[, -1]
      pie_parts <- colnames(pie_value)

      pie_value[is.na(pie_value)] <- 0
      pie_value$`__others` <- ifelse(rowSums(pie_value) > 0, 0, 1)
      pie_value_list <- as.list(pcutils::t2(pie_value))

      default_pie_color <- setNames(pcutils::get_cols(length(pie_parts), "col1"), pie_parts)
      if (is.null(pie_color)) pie_color <- default_pie_color
      if (is.null(names(pie_color))) {
        pie_color <- setNames(pie_color, pie_parts)
      } else {
        pie_color <- pcutils::update_param(default_pie_color, pie_color)
      }

      pie_color <- pie_color[pie_parts]
      pie_color <- c(pie_color, `__others` = NA)
    } else {
      pie_value_list <- pie_color <- NULL
    }
  } else {
    pie_value_list <- pie_color <- NULL
  }
  envir <- parent.frame()
  assign("pie_value_list", pie_value_list, envir)
  assign("pie_color", pie_color, envir)
}

get_show_labels <- function(tmp_v, labels_num) {
  name <- size <- NULL
  if (is.null(labels_num)) {
    if (nrow(tmp_v) < 20) {
      labels_num <- "all"
    } else {
      labels_num <- 0
    }
  }
  {
    if (labels_num == "all") {
      tmp_v %>% dplyr::pull(name) -> toplabel
    } else {
      if (labels_num >= 1) {
        tmp_v %>%
          dplyr::top_n(labels_num, size) %>%
          dplyr::arrange(-size) %>%
          dplyr::pull(name) %>%
          head(labels_num) -> toplabel
      } else {
        tmp_v %>%
          dplyr::top_frac(labels_num, size) %>%
          dplyr::arrange(-size) %>%
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
  X <- Y <- module <- minx <- maxx <- miny <- maxy <- NULL
  if (is.null(go)) {
    if (is.null(tmp_v) && is.null(ori_coors)) message("input `tmp_v` and `ori_coors` when `go` is null.")
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

produce_c_net_legends <- function(tmp_v, tmp_e, vertex_frame_width,
                                  legend_position, legend_number, legend_cex,
                                  node_size_text, edge_width_text,
                                  group_legend_title, group_legend_order,
                                  color_legend, color_legend_order,
                                  size_legend, size_legend_title,
                                  edge_legend, edge_legend_title, edge_legend_order,
                                  width_legend, width_legend_title,
                                  lty_legend, lty_legend_title, lty_legend_order,
                                  module_legend, module_legend_title, module_legend_order, module_color, mark_alpha,
                                  pie_legend, pie_legend_title, pie_legend_order, pie_color,
                                  ...) {
  color <- e_type <- lty <- e_class <- v_class <- shape <- left_leg_x <- right_leg_x <- NULL

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
    # ！！！考虑添加更多的形状，至少是21:25，形状为pie时添加额外legend

    if (is.null(group_legend_title)) {
      group_legend_title <- setNames(vgroups, vgroups)
    } else if (is.null(names(group_legend_title))) {
      group_legend_title <- setNames(rep(group_legend_title, len = length(vgroups)), vgroups)
    }

    this_shape <- unique(tmp_v$shape)
    if (!all(this_shape %in% names(default_v_shape))) {
      default_v_shape <- pcutils::update_param(
        default_v_shape,
        setNames(rep(21, length(this_shape)), this_shape)
      )
    }

    for (g_i in vgroups) {
      if ("count" %in% names(tmp_v)) {
        tmp_v1 <- tmp_v[tmp_v$v_group == g_i, c("v_class", "color", "shape", "count")]
      } else {
        tmp_v1 <- tmp_v[tmp_v$v_group == g_i, c("v_class", "color", "shape")]
      }

      tmp_v1$v_class <- factor(tmp_v1$v_class, levels = custom_sort(unique(tmp_v1$v_class)))

      vclass <- pcutils::change_fac_lev(tmp_v1$v_class, color_legend_order)
      vclass <- levels(vclass)

      node_cols <- dplyr::distinct(tmp_v1, color, v_class)
      node_cols <- setNames(node_cols$color, node_cols$v_class)
      node_shapes <- dplyr::distinct(tmp_v1, shape, v_class)
      node_shapes <- setNames(node_shapes$shape, node_shapes$v_class)

      if (legend_number) {
        eee <- table(tmp_v1$v_class)
        if (!is.null(attributes(tmp_v)$skeleton)) {
          eee <- setNames(tmp_v1$count, tmp_v1$v_class)
          legend_number <- FALSE
        }
        le_text <- paste(vclass, eee[vclass], sep = ": ")
      } else {
        le_text <- vclass
      }
      if (length(le_text) == 0) le_text <- ""
      legend(left_leg_x, left_leg_y,
        cex = 0.7 * legend_cex, adj = 0, pt.lwd = vertex_frame_width,
        legend = le_text, title.cex = 0.8 * legend_cex,
        title = group_legend_title[g_i], title.font = 2, title.adj = 0,
        col = "black", pt.bg = node_cols[vclass], bty = "n", pch = default_v_shape[node_shapes[vclass]]
      )

      left_leg_y <- left_leg_y - (length(vclass) * 0.12 + 0.2) * legend_cex
    }
  }

  if (module_legend) {
    tmp_v$module <- factor(tmp_v$module, levels = custom_sort(unique(tmp_v$module)))
    modules <- pcutils::change_fac_lev(tmp_v$module, module_legend_order)
    modules <- levels(modules)

    if (legend_number) {
      eee <- table(tmp_v$module)
      le_text <- paste(modules, eee[modules], sep = ": ")
    } else {
      le_text <- modules
    }
    legend(left_leg_x, left_leg_y,
      cex = 0.7 * legend_cex, adj = 0,
      legend = le_text, title.cex = 0.8 * legend_cex,
      title = module_legend_title, title.font = 2, title.adj = 0,
      fill = pcutils::add_alpha(module_color[modules], mark_alpha),
      border = module_color[modules], bty = "n"
    )

    left_leg_y <- left_leg_y - (length(modules) * 0.12 + 0.2) * legend_cex
  }

  if (pie_legend) {
    pie_color <- pie_color[names(pie_color) != "__others"]
    pies <- pcutils::change_fac_lev(names(pie_color), pie_legend_order)
    pies <- levels(pies)

    le_text <- pies
    legend(left_leg_x, left_leg_y,
      cex = 0.7 * legend_cex, adj = 0,
      legend = le_text, title.cex = 0.8 * legend_cex,
      title = pie_legend_title, title.font = 2, title.adj = 0,
      fill = pie_color[pies],
      border = "black", bty = "n"
    )
    left_leg_y <- left_leg_y - (length(pies) * 0.12 + 0.2) * legend_cex
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
    tmp_e$e_type <- factor(tmp_e$e_type, levels = custom_sort(unique(tmp_e$e_type)))
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
#' @param labels_num show how many labels, >1 indicates number, <1 indicates fraction, "all" indicates all.
#' @param vertex_size_range the vertex size range, e.g. c(1,10)
#' @param edge_width_range the edge width range, e.g. c(1,10)
#'
#' @param plot_module logical, plot module?
#' @param mark_module logical, mark the modules?
#' @param mark_color mark color
#' @param mark_alpha mark fill alpha, default 0.3
#' @param module_label show module label?
#' @param module_label_cex module label cex
#' @param module_label_color module label color
#' @param module_label_just module label just, default c(0.5,0.5)
#'
#' @param pie_value a dataframe using to plot pie (with rowname or a "name" column)
#' @param pie_color color vector
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
#' @param module_legend logical
#' @param module_legend_title module_legend_title
#' @param module_legend_order module_legend_order
#'
#' @param pie_legend logical
#' @param pie_legend_title pie_legend_title
#' @param pie_legend_order pie_legend_order
#'
#' @param seed random seed, default:1234, make sure each plot is the same.
#' @param params_list a list of parameters, e.g. list(edge_legend = TRUE, lty_legend = FALSE), when the parameter is duplicated, the format argument will be used rather than the argument in params_list.
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
c_net_plot <- function(go, coors = NULL, ..., labels_num = NULL,
                       vertex_size_range = NULL, edge_width_range = NULL,
                       plot_module = FALSE,
                       mark_module = FALSE, mark_color = NULL, mark_alpha = 0.3,
                       module_label = FALSE, module_label_cex = 2, module_label_color = "black",
                       module_label_just = c(0.5, 0.5),
                       pie_value = NULL, pie_color = NULL,
                       legend = TRUE, legend_number = FALSE, legend_cex = 1,
                       legend_position = c(left_leg_x = -2, left_leg_y = 1, right_leg_x = 1.2, right_leg_y = 1),
                       group_legend_title = NULL, group_legend_order = NULL,
                       color_legend = TRUE, color_legend_order = NULL,
                       size_legend = FALSE, size_legend_title = "Node Size",
                       edge_legend = TRUE, edge_legend_title = "Edge type", edge_legend_order = NULL,
                       width_legend = FALSE, width_legend_title = "Edge width",
                       lty_legend = FALSE, lty_legend_title = "Edge class", lty_legend_order = NULL,
                       module_legend = FALSE, module_legend_title = "Module", module_legend_order = NULL,
                       pie_legend = FALSE, pie_legend_title = "Pie part", pie_legend_order = NULL,
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

  if (length(V(go)) == 0) {
    message("The network is empty.")
    return(invisible())
  }
  new_modu <- module_color <- node_size_text <- edge_width_text <- NULL
  set.seed(seed)

  if (!is_metanet(go)) go <- c_net_update(go)
  # go <- c_net_update(go)

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
    edge_curved <- NULL
  } else {
    edge_curved <- ori_coors$curved$curved
  }

  # scale the size and width
  scale_size_width(tmp_v, tmp_e, vertex_size_range, edge_width_range)

  # some custom parameters
  some_custom_paras(tmp_v, tmp_e, ...)

  # set pie
  pie_set_for_plot(tmp_v, pie_value, pie_color)
  if (is.null(pie_value) || !"pie" %in% tmp_v$shape) pie_legend <- FALSE
  if (is.null(pie_color)) {
    pie_color_list <- NULL
  } else {
    pie_color_list <- list(pie_color)
  }

  # show labels
  tmp_v <- get_show_labels(tmp_v, labels_num)

  # modules set
  module_set_for_plot(tmp_v, mark_module, mark_color)
  if (!mark_module) module_legend <- FALSE

  if (any(igraph::is.loop(go))) go <- clean_multi_edge_metanet(go)
  # main plot
  {
    old_xpd <- graphics::par(mar = c(4, 2, 2, 2), xpd = TRUE)
    on.exit(graphics::par(old_xpd), add = TRUE)

    igraph::plot.igraph(go,
      layout = coors,
      vertex.size = tmp_v$size,
      vertex.color = tmp_v$color,
      vertex.shape = tmp_v$shape,
      vertex.label.color = tmp_v$label.color,
      edge.color = tmp_e$color,
      edge.lty = tmp_e$lty,
      edge.width = tmp_e$width,
      mark.groups = new_modu,
      mark.col = pcutils::add_alpha(module_color[names(new_modu)], mark_alpha),
      mark.border = module_color[names(new_modu)],
      vertex.pie = pie_value_list,
      vertex.pie.color = pie_color_list,
      ...,
      vertex.frame.width = 0.5,
      main = main,
      vertex.label.font = 1,
      vertex.label.cex = 0.07 * tmp_v$size,
      vertex.label = tmp_v$label,
      edge.arrow.size = 0.3 * tmp_e$width * 3,
      edge.arrow.width = 0.6 * tmp_e$width * 3,
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

  if ("vertex.frame.width" %in% names(list(...))) {
    vertex_frame_width <- list(...)[["vertex.frame.width"]]
  } else {
    vertex_frame_width <- 0.5
  }
  # produce legends
  if (grepl("skeleton", main)) attributes(tmp_v)$skeleton <- TRUE
  produce_c_net_legends(
    tmp_v, tmp_e, vertex_frame_width,
    legend_position, legend_number, legend_cex,
    node_size_text, edge_width_text,
    group_legend_title, group_legend_order,
    color_legend, color_legend_order,
    size_legend, size_legend_title,
    edge_legend, edge_legend_title, edge_legend_order,
    width_legend, width_legend_title,
    lty_legend, lty_legend_title, lty_legend_order,
    module_legend, module_legend_title, module_legend_order, module_color, mark_alpha,
    pie_legend, pie_legend_title, pie_legend_order, pie_color, ...
  )
}
