# =========2.build======

#' Print method for 'metanet' objects
#'
#' @param x 'metanet' object
#' @param ... Additional arguments
#' @exportS3Method
#' @method print metanet
#' @return No value
print.metanet <- function(x, ...) {
    pcutils::dabiao(" metanet ", print = TRUE)
    print.igraph(x)
}


#' Construct a network from correlation table
#'
#' @param corr result of `c_net_cal`
#' @param r_thres r_threshold (default: >0.6)
#' @param p_thres p_threshold (default: <0.05)
#' @param use_p_adj use the p.adjust instead of p-value (default: TRUE)
#' @param del_single should delete single vertexes?
#'
#' @return an metanet object
#' @export
#'
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab) -> totu
#' metadata[, 3:10] -> env
#' c_net_cal(totu) -> corr
#' c_net_build(corr, r_thres = 0.65) -> co_net
#'
#' c_net_cal(totu, env) -> corr2
#' c_net_build(corr2) -> co_net2
c_net_build <- function(corr, r_thres = 0.6, p_thres = 0.05, use_p_adj = TRUE, del_single = TRUE) {
    stopifnot(inherits(corr, "corr"))
    # set thresholds to construct
    if ("r" %in% names(corr)) {
        occor.r <- corr$r
    } else {
        stop("No r in the input object")
    }

    if (("p.adjust" %in% names(corr)) & use_p_adj) {
        occor.p <- corr$p.adjust
    } else {
        if ("p.value" %in% names(corr)) {
            occor.p <- corr$p.value
            message("Have not do p-value adjust! use the p.value to build network")
        } else {
            occor.p <- occor.r
            occor.p[occor.p != 0] <- 0
            message("No p.value given, just use r threshold to build network!")
        }
    }

    occor.r[occor.p > p_thres | abs(occor.r) < r_thres] <- 0

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
    if (del_single) go <- igraph::delete.vertices(go, V(go)[igraph::degree(go) == 0])

    # set vertex attributes
    # set vertices shape
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
            tmp <- reshape2::melt(corr$p.value, varnames = c("from", "to"), value.name = "p.value")
            tmp_e <- dplyr::left_join(tmp_e, tmp, by = c("from", "to"))
        }
        if ("p.adjust" %in% names(corr)) {
            # E(go)$p.adjust <-get_e(go)%>%dplyr::select(from,to)%>%apply(., 1, \(x)occor.r$p.adjust[x[1],x[2]])
            tmp <- reshape2::melt(corr$p.adjust, varnames = c("from", "to"), value.name = "p.adjust")
            tmp_e <- dplyr::left_join(tmp_e, tmp, by = c("from", "to"))
        }
        # 应该直接expand，再left_join快很多
        igraph::edge.attributes(go) <- as.list(tmp_e)
    }

    # set edges type
    E(go)$e_type <- ifelse(go.weight > 0, "positive", "negative")
    # set edges width
    E(go)$width <- E(go)$weight

    # set edges from_to
    anno_edge(go, get_v(go)[, c("name", "v_group")]) -> go

    # set edge intra-inter
    tmp_e <- igraph::edge.attributes(go)
    E(go)$e_class <- ifelse(tmp_e$v_group_from == tmp_e$v_group_to, "intra", "inter")
    c_net_update(go) -> go

    return(go)
}


#' Multi-omics network build
#'
#' @param ... some omics abundance tables
#' @param mode "full"
#' @param method spearman, pearson
#' @param filename the prefix of saved .corr file or FALSE
#' @param p.adjust.method see \code{\link[stats]{p.adjust}}
#' @param r_thres r_threshold (default: >0.6)
#' @param p_thres p_threshold (default: <0.05)
#' @param use_p_adj use the p.adjust instead of p-value (default: TRUE)
#' @param del_single should delete single vertexes?
#'
#' @return metanet
#' @export
#'
#' @examples
#' data("multi_test")
#' multi1 <- multi_net_build(list(Microbiome = micro, Metabolome = metab, Transcriptome = transc))
#' multi1 <- c_net_set(multi1, micro_g, metab_g, transc_g,
#'     vertex_class = c("Phylum", "kingdom", "type")
#' )
#' multi1 <- c_net_set(multi1, data.frame("Abundance1" = colSums(micro)),
#'     data.frame("Abundance2" = colSums(metab)), data.frame("Abundance3" = colSums(transc)),
#'     vertex_size = paste0("Abundance", 1:3)
#' )
#' c_net_plot(multi1)
multi_net_build <- function(..., mode = "full", method = "spearman",
                            filename = FALSE, p.adjust.method = NULL,
                            r_thres = 0.6, p_thres = 0.05, use_p_adj = TRUE, del_single = TRUE) {
    tables <- list(...)
    if (all(class(tables[[1]]) == "list")) tables <- tables[[1]]
    tabs_name <- names(tables)

    tables <- check_tabs(tables)
    if (mode == "full") {
        all_totu <- do.call(cbind, tables)
        message("Calculating ", nrow(all_totu), " samples and", ncol(all_totu), " objects of ", length(tables), " groups.\n")
        all_corr <- c_net_cal(all_totu, method = method, filename = filename, p.adjust.method = p.adjust.method)

        c_net_build(all_corr, r_thres = r_thres, p_thres = p_thres, use_p_adj = use_p_adj, del_single = del_single) -> multi_net

        igraph::graph.attributes(multi_net)$n_type <- "multi_full"

        get_v(multi_net) -> tmp_v
        if (is.null(tabs_name)) tabs_name <- paste0("omic", seq_len(length(tables)))
        position <- rep(tabs_name, vapply(tables, ncol, numeric(1)))
        names(position) <- lapply(tables, colnames) %>% do.call(c, .)

        tmp_v$v_class <- tmp_v$v_group <- vapply(tmp_v$name, \(x)position[x], character(1))
        tmp_v %>% as.list() -> igraph::vertex_attr(multi_net)

        # set edges from_to
        suppressMessages(anno_edge(multi_net, tmp_v[, c("name", "v_group")]) -> multi_net)
        # set edge intra-inter
        tmp_e <- igraph::edge.attributes(multi_net)
        E(multi_net)$e_class <- ifelse(tmp_e$v_group_from == tmp_e$v_group_to, "intra", "inter")
        c_net_update(multi_net) -> multi_net
    }
    multi_net
}

# 多列表是否一一对应
# test 对应的列需要跟剩下的列具有唯一映射
e_match <- function(df, test = NULL) {
    dis_df <- dplyr::distinct_all(df) %>% as.data.frame()
    if (is.null(test)) {
        return(all_same(c(nrow(dis_df), apply(dis_df, 2, \(i)length(unique(i))))))
    } else {
        return(all_same(c(nrow(dis_df), nrow(dplyr::distinct_all(dis_df[, test, drop = FALSE])))))
    }
}

#' Update a metanet
#'
#' @param go metanet metanet object
#' @param node_break node_break if v_class is numeric, default: 5
#' @param edge_break edge_break if e_type is numeric, default: 5
#' @aliases as.metanet
#' @export
#' @return metanet
c_net_update <- function(go, node_break = 5, edge_break = 5) {
    stopifnot(inherits(go, "igraph"))

    # name
    if (!"name" %in% igraph::vertex_attr_names(go)) V(go)$name <- paste0("n", seq_len(length(go)))
    if (!"label" %in% igraph::vertex_attr_names(go)) V(go)$label <- V(go)$name
    get_v(go) -> tmp_v
    # v_size
    if (!"size" %in% colnames(tmp_v)) {
        tmp_v$size <- ceiling(60 / sqrt(length(V(go)))) + 1
    }

    flag <- TRUE
    # v_shape
    if (!"v_group" %in% colnames(tmp_v)) {
        tmp_v$v_group <- "v_group1"
    }
    if ("shape" %in% colnames(tmp_v)) {
        if (e_match(tmp_v[, c("v_group", "shape")], test = 1)) flag <- FALSE
    }
    if (flag) tmp_v$shape <- pcutils::tidai(tmp_v$v_group, c("circle", "square"))

    flag <- TRUE
    # v_color
    if (!"v_class" %in% colnames(tmp_v)) {
        tmp_v$v_class <- "v_class1"
    }
    if ("color" %in% colnames(tmp_v)) {
        if (e_match(tmp_v[, c("v_group", "v_class", "color")], test = 1:2)) flag <- FALSE
    }
    if (flag) {
        if (is.numeric(tmp_v$v_class)) {
            tmp_v$color <- as.character(cut(tmp_v$v_class,
                breaks = node_break,
                labels = pcutils::get_cols(node_break, RColorBrewer::brewer.pal(n = 9, name = "Reds")[2:9])
            ))
            tmp_v$v_class <- as.character(cut(tmp_v$v_class,
                breaks = node_break,
                labels = seq(min(tmp_v$v_class, na.rm = T), max(tmp_v$v_class, na.rm = T), length.out = node_break)
            ))
        } else {
            tmp_col <- paste0(tmp_v$v_group, "-", tmp_v$v_class)
            tmp_v$color <- pcutils::tidai(tmp_col, pcutils::get_cols(nlevels(factor(tmp_col)), "col3"), fac = TRUE)
        }
    }

    as.list(tmp_v) -> igraph::vertex_attr(go)

    flag <- TRUE
    # e_color
    if (!"e_type" %in% igraph::edge_attr_names(go)) {
        E(go)$e_type <- "e_type1"
    }
    if ("color" %in% igraph::edge_attr_names(go)) {
        if (e_match(data.frame(a = E(go)$e_type, b = E(go)$color))) flag <- FALSE
    }
    if (flag) {
        if (is.numeric(E(go)$e_type)) {
            E(go)$color <- as.character(cut(E(go)$e_type,
                breaks = edge_break,
                labels = pcutils::get_cols(edge_break, RColorBrewer::brewer.pal(n = 9, name = "Greens")[3:9])
            ))
            E(go)$e_type <- as.character(cut(E(go)$e_type, breaks = edge_break))
        } else {
            edge.color <- droplevels(as.factor(E(go)$e_type))
            if (all(levels(edge.color) %in% c("negative", "positive"))) {
                ncols <- c(negative = "#E85D5D", positive = "#48A4F0")
            } else if (all(levels(edge.color) %in% c("inter-module", "intra-module"))) {
                ncols <- c("inter-module" = "#FA789A", "intra-module" = "#A6CEE3")
            } else {
                ncols <- pcutils::get_cols(nlevels(edge.color), "col2")
            }
            E(go)$color <- pcutils::tidai(E(go)$e_type, ncols, fac = TRUE)
        }
    }

    flag <- TRUE
    # e_lty
    if (!"e_class" %in% igraph::edge_attr_names(go)) {
        E(go)$e_class <- "e_class1"
    }
    if ("lty" %in% igraph::edge_attr_names(go)) {
        if (e_match(data.frame(a = E(go)$e_class, b = E(go)$lty), test = 1)) flag <- FALSE
    }
    if (flag) E(go)$lty <- pcutils::tidai(E(go)$e_class, 1:4)

    # e_width
    if (!"width" %in% igraph::edge_attr_names(go)) {
        E(go)$width <- 1
    }

    class(go) <- c("metanet", "igraph")
    return(go)
}

#' @export
assign("as.metanet", c_net_update, envir = asNamespace(packageName()))

#' Construct a network from edge_list dataframe
#'
#' @param edgelist first is source, second is target, others are annotation
#' @param vertex vertex metadata
#' @param direct logical
#' @param e_type set e_type
#' @param e_class set e_class
#' @return metanet
#' @export
#' @examples
#' data(edgelist)
#' edge_net <- c_net_from_edgelist(arc_count, vertex = arc_taxonomy)
#' c_net_plot(edge_net)
c_net_from_edgelist <- function(edgelist, vertex = NULL, direct = FALSE, e_type = NULL, e_class = NULL) {
    lib_ps("igraph", library = FALSE)
    go <- igraph::graph_from_data_frame(edgelist, directed = direct, vertices = vertex)
    if (!is.null(e_type)) E(go)$e_type <- edgelist[, e_type]
    if (!is.null(e_class)) E(go)$e_class <- edgelist[, e_class]
    go <- c_net_update(go)
    go
}


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
#'
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab) -> totu
#' metadata[, 3:10] -> env
#'
#' data("c_net")
#' co_net <- c_net_set(co_net, taxonomy, data.frame("Abundance" = colSums(totu)),
#'     vertex_class = "Phylum", vertex_size = "Abundance"
#' )
#' co_net2 <- c_net_set(co_net2, taxonomy, data.frame(name = colnames(env), env = colnames(env)),
#'     vertex_class = c("Phylum", "env")
#' )
#' co_net2 <- c_net_set(co_net2, data.frame("Abundance" = colSums(totu)), vertex_size = "Abundance")
c_net_set <- function(go, ..., vertex_group = "v_group", vertex_class = "v_class", vertex_size = "size",
                      edge_type = "e_type", edge_class = "e_class", edge_width = "width", node_break = 5, edge_break = 5) {
    size <- e_class <- width <- NULL
    lib_ps("igraph", library = FALSE)
    c_net_update(go) -> go
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
    if (!setequal(vertex_group, "v_group")) dplyr::select(v_index, v_group, !!vertex_group) %>% condance() -> v_index$v_group
    if (!setequal(vertex_class, "v_class")) {
        tmp_color <- list()
        tmp_ann <- dplyr::select(v_index, !!vertex_class)
        # charcter
        tmp_col <- dplyr::select_if(tmp_ann, \(i)!is.numeric(i)) %>%
            condance() %>%
            unique() %>%
            na.omit()
        if (length(tmp_col) > 0) tmp_col <- setNames(pcutils::get_cols(nlevels(factor(tmp_col)), "col3"), tmp_col)
        # tmp_col=c(tmp_col,`NA`=NA)

        # 可能某一个group用numeric做v_class，所以要分开上色
        for (i in 1:ncol(tmp_ann)) {
            if (is.numeric(tmp_ann[, i])) {
                tmp_color[[i]] <- as.character(cut(tmp_ann[, i],
                    breaks = node_break,
                    labels = pcutils::get_cols(node_break, RColorBrewer::brewer.pal(n = 9, name = "Reds")[2:9])
                ))
                tmp_ann[, i] <- as.character(cut(tmp_ann[, i],
                    breaks = node_break,
                    labels = seq(min(tmp_ann[, i], na.rm = TRUE), max(tmp_ann[, i], na.rm = TRUE),
                        length.out = node_break
                    )
                ))
            } else {
                tmp_color[[i]] <- tmp_col[tmp_ann[, i]]
            }
        }

        cbind(dplyr::select(v_index, v_class), tmp_ann) %>% condance() -> v_index$v_class
        cbind(dplyr::select(v_index, color), as.data.frame(tmp_color)) %>% condance() -> v_index$color
    }
    if (!setequal(vertex_size, "size")) dplyr::select(v_index, size, !!vertex_size) %>% condance() -> v_index$size
    if (!setequal(edge_type, "e_type")) {
        tmp_color <- list()
        tmp_ann <- dplyr::select(e_index, !!edge_type)
        # charcter
        tmp_col <- dplyr::select_if(tmp_ann, \(i)!is.numeric(i)) %>%
            condance() %>%
            unique() %>%
            na.omit()
        if (length(tmp_col) > 0) {
            tmp_col <- setNames(pcutils::get_cols(nlevels(factor(tmp_col)), "col2"), tmp_col)
            tmp_col <- pcutils::update_param(tmp_col, c(
                negative = "#E85D5D", positive = "#48A4F0",
                "inter-module" = "#FA789A", "intra-module" = "#A6CEE3"
            ))
        }
        # tmp_col=c(tmp_col,`NA`=NA)

        for (i in 1:ncol(tmp_ann)) {
            if (is.numeric(tmp_ann[, i])) {
                tmp_color[[i]] <- as.character(cut(tmp_ann[, i],
                    breaks = edge_break,
                    labels = pcutils::get_cols(edge_break, RColorBrewer::brewer.pal(n = 9, name = "Greens")[3:9])
                ))
                tmp_ann[, i] <- as.character(cut(tmp_ann[, i],
                    breaks = edge_break,
                    labels = seq(min(tmp_ann[, i], na.rm = T), max(tmp_ann[, i], na.rm = T),
                        length.out = edge_break
                    )
                ))
            } else {
                tmp_color[[i]] <- tmp_col[tmp_ann[, i]]
            }
        }

        cbind(dplyr::select(e_index, e_type), tmp_ann) %>% condance() -> e_index$e_type
        cbind(dplyr::select(e_index, color), as.data.frame(tmp_color)) %>% condance() -> e_index$color
    }
    if (!setequal(edge_class, "e_class")) dplyr::select(e_index, e_class, !!edge_class) %>% condance() -> e_index$e_class
    if (!setequal(edge_width, "width")) dplyr::select(e_index, width, !!edge_width) %>% condance() -> e_index$width

    as.list(v_index) -> igraph::vertex.attributes(go)
    as.list(e_index) -> igraph::edge.attributes(go)

    c_net_update(go) -> go
    return(go)
}


all_same <- \(x){
    return(all(x == x[1]))
}

# choose the last not na value
condance <- \(aa){
    if (any(is.na(aa[, length(aa)]))) {
        res <- apply(aa, 1, \(x){
            tmp <- x[!is.na(x)]
            if (length(tmp) == 0) {
                return(NA)
            }
            return(tmp[length(tmp)])
        })
    } else {
        res <- aa[, length(aa)]
    }
    res
}


#' Transfer a dataframe to a network edgelist.
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
    graph_attr(net, "coors") <- c_net_lay(net, as_tree())
    net
}

# ==========2.1 manipulate========

#' Is this object an metanet graph?
#'
#' @param go igraph object
#'
#' @return logical
#' @export
#' @aliases is.metanet
#' @examples
#' data(c_net)
#' is_metanet(co_net)
is_metanet <- function(go) {
    is.igraph(go) & inherits(go, "metanet")
}

#' Get vertex information
#' @param go metanet
#' @return data.frame
#' @export
get_v <- function(go) {
    as.data.frame(igraph::vertex.attributes(go))
}

# bind two df with same columns, the last df will replace first df
cbind_new <- \(df, df1){
    if (ncol(df) < 1) {
        return(df1)
    }
    if (ncol(df1) < 1) {
        return(df)
    }
    inter <- intersect(colnames(df1), colnames(df))
    la <- setdiff(colnames(df), inter)
    cbind(df[, la, drop = F], df1)
}

#' Get edge information
#' @param go metanet
#' @return data.frame
#' @export
get_e <- function(go) {
    id <- NULL
    tmp_e <- cbind_new(igraph::as_data_frame(go), data.frame(id = seq_len(igraph::ecount(go))))
    dplyr::select(tmp_e, id, dplyr::everything())
}


#' Get network information
#' @param go metanet
#' @param simple logical, get simple index
#' @return data.frame
#' @export
get_n <- function(go, simple = FALSE) {
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
    as.data.frame(do.call(cbind, gls))
}

#' Filter a network according to some attributes
#'
#' @param go metanet
#' @param ... some attributes of vertex and edge
#' @param mode "v" or "e"
#'
#' @return igraph
#' @export
#'
#' @examples
#' data("multi_net")
#' c_net_filter(multi1, v_group %in% c("omic1", "omic2"))
c_net_filter <- function(go, ..., mode = "v") {
    if (mode == "v") {
        go1 <- filter_v(go, ...)
    } else if (mode == "e") {
        go1 <- filter_e(go, ...)
    }
    class(go1) <- c("metanet", "igraph")
    go1
}

# c_net_filter1<-function(go,...,exclude=FALSE){
#   f_ls=list(...)
#   v_ls=list();e_ls=list()
#   for (i in names(f_ls)) {
#     if((i %in% igraph::vertex_attr_names(go))){
#       if((i %in% igraph::edge_attr_names(go))){warning("both vertex and edge has attributes ",i," use `filter_e()` to filter edge")}
#       v_ls[[i]]=f_ls[[i]]
#     }
#     else if(i %in% igraph::edge_attr_names(go)){
#       e_ls[[i]]=f_ls[[i]]
#     }
#   }
#   go1=filter_v(go,v_ls,exclude=exclude)
#   go1=filter_e(go1,e_ls,exclude=exclude)
#   class(go1)=c("metanet","igraph")
#   go1
# }


filter_v <- function(go, ...) {
    get_v(go) -> tmp_v
    tmp_v <- dplyr::filter(tmp_v, ...)
    tmp_v$name -> vid
    igraph::subgraph(go, vid) -> go1
    class(go1) <- c("metanet", "igraph")
    go1
}

# filter_v1<-function(go,...,exclude=FALSE){
#   f_ls=list(...)
#   if(is.list(f_ls[[1]]))f_ls=f_ls[[1]]
#   get_v(go)->tmp_v
#   for (i in names(f_ls)) {
#     if(!i%in%colnames(tmp_v)){warning(i," do not in vertex_attributes, skip");next}
#     if(exclude)tmp_v=tmp_v[!tmp_v[[i]]%in%f_ls[[i]],]
#     else tmp_v=tmp_v[tmp_v[[i]]%in%f_ls[[i]],]
#   }
#
#   tmp_v$name->vid
#   igraph::subgraph(go,vid)->go1
#   class(go1)=c("metanet","igraph")
#   go1
# }


filter_e <- function(go, ...) {
    get_e(go) -> tmp_e
    tmp_e <- dplyr::filter(tmp_e, ...)
    tmp_e$id -> eid
    igraph::subgraph.edges(go, eid) -> go1
    class(go1) <- c("metanet", "igraph")
    go1
}

# filter_e1<-function(go,...,exclude=FALSE){
#   f_ls=list(...)
#   if(is.list(f_ls[[1]]))f_ls=f_ls[[1]]
#   get_e(go)->tmp_e
#   for (i in names(f_ls)) {
#     if(!i%in%colnames(tmp_e)){warning(i," do not in edge_attributes, skip");next}
#     if(exclude)tmp_e=tmp_e[!tmp_e[[i]]%in%f_ls[[i]],]
#     else tmp_e=tmp_e[tmp_e[[i]]%in%f_ls[[i]],]
#   }
#
#   tmp_e$id->eid
#   igraph::subgraph.edges(go,eid)->go1
#   class(go1)=c("metanet","igraph")
#   go1
# }

#' Use data.frame to annotate vertexes of an igraph
#'
#' @param go metanet an igraph object
#' @param anno_tab a dataframe using to annotate (with rowname or a name column)
#' @return a annotated igraph object
#' @aliases anno_node
#' @export
#'
#' @examples
#' data("c_net")
#' data("otutab", package = "pcutils")
#' anno_vertex(co_net, taxonomy)
anno_vertex <- function(go, anno_tab) {
    lib_ps("igraph", library = FALSE)
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
#' @param anno_tab a dataframe using to annotate (with rowname or a name column)
#' @return a annotated igraph object
#' @export
#'
#' @examples
#' data("c_net")
#' data("otutab", package = "pcutils")
#' anno_edge(co_net, taxonomy) -> a
anno_edge <- function(go, anno_tab) {
    name <- NULL
    lib_ps("igraph", library = FALSE)
    if (is.null(anno_tab)) {
        return(go)
    }
    get_e(go) -> e_atr
    if (!"name" %in% colnames(anno_tab)) rownames(anno_tab) -> anno_tab$name
    anno_tab %>% dplyr::select(name, dplyr::everything()) -> anno_tab
    # from
    tmp <- anno_tab
    colnames(tmp) <- paste0(colnames(anno_tab), "_from")
    e_atr <- dplyr::left_join(e_atr, tmp, by = c("from" = "name_from"), suffix = c(".x", ""))
    grep(".x", colnames(e_atr), value = TRUE) %>% gsub(".x", "", .) -> du
    if (length(du) > 0) message(length(du), (" attributes will be overwrited:\n"), paste0(du, collapse = ","), "\n")
    e_atr %>% dplyr::select(!dplyr::ends_with(".x")) -> e_atr
    # to
    tmp <- anno_tab
    colnames(tmp) <- paste0(colnames(anno_tab), "_to")
    e_atr <- dplyr::left_join(e_atr, tmp, by = c("to" = "name_to"), suffix = c(".x", ""))
    grep(".x", colnames(e_atr), value = TRUE) %>% gsub(".x", "", .) -> du
    if (length(du) > 0) message(length(du), (" attributes will be overwrited:\n"), paste0(du, collapse = ","), "\n")
    e_atr %>% dplyr::select(!dplyr::ends_with(".x")) -> e_atr

    as.list(e_atr) -> igraph::edge.attributes(go)
    return(go)
}

#' Save network file
#'
#' @param go metanet network
#' @param filename filename
#' @param format "data.frame","graphml"
#' @return No value
#' @export
c_net_save <- function(go, filename = "net", format = "data.frame") {
    if (format == "data.frame") {
        get_v(go) %>% write.csv(., paste0(filename, "_nodes.csv"), row.names = FALSE)
        get_e(go) %>%
            dplyr::select(-1) %>%
            write.csv(., paste0(filename, "_edges.csv"), row.names = FALSE)
    }
    if (format == "graphml") {
        go <- igraph::delete_edge_attr(go, "id")
        igraph::write.graph(go, paste0(filename, ".graphml"), format = "graphml")
    }
    print(paste0(filename, " saved sucessfully!"))
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
#'     a = sample(letters[1:4], 10, replace = TRUE),
#'     b = sample(letters[1:4], 10, replace = TRUE)
#' )
#' summ_2col(test, direct = TRUE)
#' summ_2col(test, direct = FALSE)
#' summ_2col(test, direct = TRUE) %>% pcutils::my_circo()
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


# =========2.2RMT optimize=====

#' Get RMT threshold for a correlation matrix
#'
#' @param occor.r a corr object or a correlation matrix
#' @param min_threshold min_threshold
#' @param max_threshold max_threshold
#' @param step step
#' @param gif render a .gif file?
#' @param verbose verbose
#' @param out_dir output dir
#'
#' @return a r-threshold
#' @export
#' @references 1.  J. Zhou, Y. Deng, FALSE. Luo, Z. He, Q. Tu, X. Zhi, Functional Molecular Ecological Networks (2010), doi:10.1128/mBio.00169-10.
#' @references 2. \code{https://matstat.org/content_en/RMT/RMThreshold_Intro.pdf}
#' @examples
#' \donttest{
#' data(otutab, package = "pcutils")
#' t(otutab) -> totu
#' c_net_cal(totu) -> corr
#' rmt(corr)
#' # recommend: 0.69
#' c_net_build(corr, r_thres = 0.69) -> co_net_rmt
#' }
RMT_threshold <- function(occor.r, out_dir, min_threshold = 0.5, max_threshold = 0.8,
                          step = 0.02, gif = FALSE, verbose = FALSE) {
    nwd <- getwd()
    on.exit(setwd(nwd))

    setwd(out_dir)
    if (inherits(occor.r, "corr")) occor.r <- occor.r$r
    if (!dir.exists("./RMT_temp")) dir.create("./RMT_temp")
    diag(occor.r) <- 0

    if (max_threshold >= max(abs(occor.r))) max_threshold <- (max(abs(occor.r)) - step)
    if (min_threshold >= max_threshold) min_threshold <- max_threshold - 10 * step

    thres_seq <- seq(min_threshold, max_threshold, step)

    res <- data.frame()
    for (i in seq_len(length(thres_seq))) {
        threshold <- thres_seq[i]
        if (!verbose) pcutils::dabiao(paste0("Calculating", i, ":  threshold =", signif(threshold, 3)), print = TRUE)
        corr_r1 <- occor.r
        corr_r1[abs(corr_r1) < threshold] <- 0
        # calculate eigenvalues
        rand.mat <- corr_r1
        eigenvalues <- eigen(rand.mat, only.values = TRUE)$values
        eigenvalues <- eigenvalues[order(eigenvalues)] / max(abs(eigenvalues))

        eigenvalues <- pcutils::remove.outliers(unique(eigenvalues))

        # get the NNDS
        { # uf <- rm.unfold.gauss(eigenvalues,pop.up = TRUE)
            dens <- density(eigenvalues, kernel = "gaussian")
            midpoints <- \(x)(x[-length(x)] + 0.5 * diff(x))
            scale.function <- approx(dens$x, dens$y, xout = midpoints(eigenvalues))
            ev.spacing <- diff(eigenvalues)
            ev.spacing <- ev.spacing * scale.function$y
            ev.spacing <- ev.spacing / mean(ev.spacing)
        }

        # if (FALSE) { uf <- rm.unfold.spline(eigenvalues,pop.up = TRUE)
        #     cumulative.distrib <- ecdf(eigenvalues)
        #     support <- seq(min(eigenvalues), max(eigenvalues), len = nr.fit.points)
        #     D <- splinefun(support, cumulative.distrib(support), method = "hyman")
        #     uf.ev <- D(eigenvalues)
        #     ev.spacing <- diff(uf.ev)
        #     ev.spacing <- ev.spacing / mean(ev.spacing)
        # }
        # if (FALSE) {
        #     ssp <- smooth.spline(eigenvalues, control.spar = list(low = 0, high = 3))
        #     ev.spacing <- NNSD(ssp$y) / mean(NNSD(ssp$y))
        # }

        ev.spacing <- ev.spacing[ev.spacing <= 3]
        # test whether fit possion?
        p_ks_test <- ks.test(unique(ev.spacing), "pexp", 1)$p.value
        # get sse
        # sse = rm.sse(ev.spacing)
        sse <- get_sse(ev.spacing)
        log_sse <- log(sse)

        # if (FALSE) {
        #     nnsdw <- density(ev.spacing)
        #     poisson_d <- exp(seq(0, 3, len = 1000))
        #     nnsdpois <- density(NNSD(poisson_d))
        #     chival <- sum((nnsdw$y - nnsdpois$y)^2 / nnsdpois$y / 1e3)
        # }

        # maximum likelihood
        evs <- ev.spacing[ev.spacing != 0]
        N <- length(evs)
        log_LE <- -sum(evs) / N
        log_LW <- log(pi / 2) + sum(log(evs)) / N - 0.25 * pi * sum(evs^2) / N

        # save png
        {
            histo <- hist(ev.spacing, breaks = seq(min(ev.spacing), max(ev.spacing), len = 51), plot = FALSE)
            grDevices::png(paste0("RMT_temp/rmt_nnsd", i, ".png"), height = 600, width = 700, res = 130)
            nnsd_plot(
                histo = histo, title = "Eigenvalue spacing distribution (NNSD)", threshold = threshold,
                dis_GOE = log_LW, dis_possion = log_LE, p_ks_test = p_ks_test
            )
            grDevices::dev.off()
        }
        res <- rbind(res, data.frame(threshold, p_ks_test, log_sse, log_LW, log_LE))
    }
    message("The Intermediate files saved in ", out_dir, "/RMT_temp/.")
    # transfer to gif
    if (gif) {
        lib_ps("gifski", library = FALSE)
        gifski::gifski(paste0("RMT_temp/rmt_nnsd", seq_len(length(thres_seq)), ".png"),
            gif_file = "RMT_temp/rmt_nnsd.gif"
        )
    }
    r_threshold <- (res[which(res$log_LW == min(res$log_LW)), "threshold"] +
        res[which(res$log_LE == max(res$log_LE)), "threshold"]) / 2
    res <- list(res = res, r_threshold = r_threshold)
    class(res) <- c("rmt_res", class(res))
    return(res)
}

#' Plot a rmt_res
#'
#' @param x rmt_res
#' @param ... Additional arguments
#'
#' @return ggplot
#' @exportS3Method
#' @method plot rmt_res
plot.rmt_res <- function(x, ...) {
    threshold <- value <- variable <- xi <- y <- NULL
    res <- x$res
    linedf <- data.frame(
        variable = c("p_ks_test", "log_sse", "log_LW", "log_LE"),
        xi = c(
            res[which(res$p_ks_test == max(res$p_ks_test))[1], "threshold"],
            res[which(res$log_sse == min(res$log_sse))[1], "threshold"],
            res[which(res$log_LW == min(res$log_LW))[1], "threshold"],
            res[which(res$log_LE == max(res$log_LE))[1], "threshold"]
        ),
        x = max(res$threshold) - min(res$threshold),
        y = apply(res[, -1], 2, max)
    )

    reshape2::melt(res, "threshold") -> md

    # filter(threshold<0.77)%>%
    p <- ggplot(md, aes(threshold, value)) +
        geom_point(aes(col = variable)) +
        geom_line(aes(col = variable)) +
        scale_color_manual(values = get_cols(4, "col1")) +
        facet_wrap(. ~ variable, scales = "free_y") +
        theme_bw() +
        xlab(NULL) +
        geom_text(data = linedf, aes(x = xi - 0.02 * x, y = 0.5 * y, label = xi)) +
        geom_vline(data = linedf, aes(xintercept = xi), linetype = 2, col = "red") +
        theme(legend.position = "none")

    message(paste("recommend r_threshold: ", mean(linedf$xi)))
    return(p)
}

nnsd_plot <- \(histo = histo, title = title, threshold = threshold,
    dis_GOE = dis_GOE, dis_possion = dis_possion, p_ks_test = p_ks_test) {
    plot(histo, freq = FALSE, col = "#F4FCA1", main = title, font.main = 1, xlab = "eigenvalue spacing", ylab = "PDF of eigenvalue spacing")
    {
        actual.ymax <- par("yaxp")[2]
        x0 <- -log(actual.ymax * 0.98)
        possion_dis <- \(x)exp(-x)
        graphics::curve(possion_dis,
            from = max(x0, min(histo$breaks)),
            to = max(histo$breaks), n = 1001, add = TRUE, type = "l", lty = 1, col = "#EB34FF", lwd = 2
        )
    }
    {
        GOE <- function(x) pi / 2 * x * exp(-pi / 4 * x^2)
        graphics::curve(GOE,
            from = min(histo$breaks),
            to = max(histo$breaks), n = 1001, add = TRUE, type = "l",
            lty = 1, col = "blue", lwd = 2
        )
    }

    if ((!is.na(dis_GOE)) && (!is.na(dis_possion))) {
        graphics::mtext(side = 3, paste(
            "Distance to GOE =", signif(dis_GOE, 3),
            "\nDistance to Possion =", signif(dis_possion, 3), "; ks_test p.value for possion =", signif(p_ks_test, 3)
        ), col = "#878787", cex = 0.6)
    }

    if (!is.na(threshold)) graphics::mtext(side = 4, paste("threshold =", signif(threshold, 4)))

    graphics::legend("topright", inset = 0.05, c("Possion", "GOE"), col = c("#EB34FF", "blue"), lty = 1, lwd = 2, cex = 0.8)
}

trapez <- \(x, y){
    ind <- 2:length(x)
    as.double((x[ind] - x[ind - 1]) %*% (y[ind] + y[ind - 1])) / 2
}

get_sse <- \(ev.spacing){
    dens <- density(ev.spacing)
    N <- 20
    x <- seq(min(ev.spacing), max(ev.spacing), len = 1000)
    observed <- approx(dens$x, dens$y, xout = x)$y
    expected <- exp(-x)
    A <- exp(-min(ev.spacing)) - exp(-max(ev.spacing))
    xs <- numeric(N + 1)
    xs[1] <- min(ev.spacing)
    for (i in 1:N) xs[i + 1] <- -log(exp(-xs[i]) - A / N)
    area <- numeric(N)
    for (i in 1:N) {
        xsec <- x[(x > xs[i]) & (x < xs[i + 1])]
        xsec <- c(xs[i], xsec, xs[i + 1])
        ysec <- approx(dens$x, dens$y, xout = xsec)$y
        area[i] <- trapez(xsec, ysec)
    }
    sse <- sum((area[i] - A / N)^2)
    sse
}

#' Get RMT threshold for a correlation matrix roughly
#'
#' @export
#' @return recommend threshold
#' @rdname RMT_threshold
rmt <- function(occor.r, min_threshold = 0.5, max_threshold = 0.85, step = 0.01) {
    if (inherits(occor.r, "corr")) occor.r <- occor.r$r
    lib_ps("ggplot2", library = F)
    NNSD <- \(x)abs(diff(x))

    s <- seq(0, 3, 0.1)
    poisson_d <- exp(-s)
    nnsdpois <- density(NNSD(poisson_d))

    ps <- c()
    threshold <- c()

    for (i in seq(min_threshold, max_threshold, step)) {
        corr_r1 <- occor.r
        corr_r1[abs(corr_r1) < i] <- 0

        {
            eigen_res <- sort(eigen(corr_r1)$value)
            # spline to eigen_res
            check <- tryCatch(ssp <- smooth.spline(eigen_res, control.spar = list(low = 0, high = 3)),
                error = \(e) {
                    TRUE
                }
            )
            if (rlang::is_true(check)) next
            nnsdw <- density(NNSD(ssp$y))
            chival <- sum((nnsdw$y - nnsdpois$y)^2 / nnsdpois$y / 1e3)
        }

        ps <- c(ps, chival)
        threshold <- c(threshold, i)
        if (((i * 100) %% 5 == 0)) {
            message(paste0("Calculating: ", i))
        }
    }

    res <- data.frame(threshold, ps)
    recommend_thres <- res[which.min(res[, 2]), 1]
    p <- ggplot(res, aes(threshold, ps)) +
        geom_point() +
        geom_vline(xintercept = recommend_thres, linetype = 2, col = "red") +
        geom_text(x = recommend_thres + 0.01, y = 0.5 * max(res$ps), label = recommend_thres) +
        theme_bw(base_size = 15)
    print(p)

    res1 <- res[(res$threshold < (recommend_thres + 0.05)) & (res$threshold > (recommend_thres - 0.05)), ]
    p <- ggplot(res1, aes(threshold, ps)) +
        geom_point() +
        geom_vline(xintercept = recommend_thres, linetype = 2, col = "red") +
        geom_text(x = recommend_thres + 0.01, y = 0.5 * max(res1$ps), label = recommend_thres) +
        theme_bw(base_size = 15)
    print(p)

    message("We recommend r-threshold: ", recommend_thres, ", you can calculate again in a smaller region")
    recommend_thres
}
