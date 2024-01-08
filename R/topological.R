# ========5.topological=======

#' Extract sub-network from the whole network
#'
#' @param a_net the whole network
#' @param otutab otutab, these columns will be extract
#' @param threads threads, default: 1
#' @param save_net should save these sub_nets? FALSE or a filename
#' @param fast less indexes for faster calculate ?
#' @param verbose verbose
#'
#' @return a dataframe contains all sub_net parameters
#' @export
#'
#' @examples
#' data(otutab, package = "pcutils")
#' extract_sub_net(co_net, otutab) -> sub_net_pars
extract_sub_net <- function(a_net, otutab, threads = 1, save_net = FALSE, fast = TRUE, verbose = TRUE) {
    lib_ps("igraph", library = FALSE)
    V(a_net)$name -> v_name
    reps <- ncol(otutab)

    message("extracting")
    sub_nets <- lapply(1:reps, \(i){
        rownames(otutab)[otutab[, i] > 0] -> exist_sp
        subgraph(a_net, which(v_name %in% exist_sp)) -> spe_sub
        class(spe_sub) <- c("metanet", "igraph")
        return(spe_sub)
    })
    names(sub_nets) <- colnames(otutab)
    message("calculating topological indexes")

    # parallel
    # main function
    loop <- function(i) {
        spe_sub <- sub_nets[[i]]
        indexs <- net_par(spe_sub, mode = "n", fast = fast)[["n_index"]]
        wc <- igraph::cluster_fast_greedy(spe_sub, weights = abs(igraph::E(spe_sub)$weight))
        indexs$modularity <- igraph::modularity(wc)
        indexs
    }
    {
        if (threads > 1) {
            pcutils::lib_ps("foreach", "doSNOW", "snow")
            if (verbose) {
                pb <- utils::txtProgressBar(max = reps, style = 3)
                opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
            } else {
                opts <- NULL
            }
            cl <- snow::makeCluster(threads)
            doSNOW::registerDoSNOW(cl)
            res <- foreach::foreach(i = 1:reps, .options.snow = opts) %dopar% {
                loop(i)
            }
            snow::stopCluster(cl)
            gc()
            pcutils::del_ps("doSNOW", "snow", "foreach")
        } else {
            res <- lapply(1:reps, loop)
        }
    }
    # simplify method
    sub_net_pars <- do.call(rbind, res)

    rownames(sub_net_pars) <- colnames(otutab)

    if (is.logical(save_net)) {
        if (save_net) save_net <- paste0("sub_net_", date())
    }
    if (is.character(save_net)) {
        saveRDS(sub_nets, file = paste0(save_net, ".RDS"))
    }
    sub_net_pars
}

#' Calculate natural_connectivity
#'
#' @param p an igraph or metanet object
#' @return natural_connectivity (numeric)
#' @export
#' @references \code{`nc` in `ggClusterNet`}
#' @examples
#' igraph::make_ring(10) %>% nc()
nc <- function(p) {
    adj_matrix <- as.matrix(igraph::as_adj(p, sparse = FALSE))
    adj_matrix[abs(adj_matrix) != 0] <- 1

    lambda <- eigen(adj_matrix, only.values = TRUE)$values
    lambda <- sort(lambda, decreasing = TRUE)

    lambda_sum <- 0
    N <- length(lambda)
    for (i in 1:N) lambda_sum <- lambda_sum + exp(lambda[i])
    lambda_average <- log(lambda_sum / N, base = exp(1))
    lambda_average
}


#' Calculate all indexes of a network
#'
#' @param go an igraph or metanet object
#' @param fast less indexes for faster calculate ?
#' @param mode calculate what? c("v", "e", "n", "all")
#'
#' @return a 3-elements list
#' \item{n_index}{indexs of the whole network}
#' \item{v_index}{indexs of each vertex}
#' \item{e_index}{indexs of each edge}
#' @export
#'
#' @examples
#' igraph::make_graph("Walther") %>% net_par()
#' c_net_index(co_net) -> co_net_with_par
net_par <- function(go, mode = c("v", "e", "n", "all"), fast = TRUE) {
    lib_ps("igraph", "dplyr", library = FALSE)
    stopifnot(is_igraph(go))
    if ("all" %in% mode) mode <- c("v", "e", "n")

    n_index <- NULL
    v_index <- NULL
    e_index <- NULL
    # non-weighted network
    up <- go
    if (!is.null(igraph::edge_attr(up)[["weight"]])) up <- igraph::delete_edge_attr(up, "weight")
    if ("n" %in% mode) {
        # Calculate Network Parameters
        n_index <- data.frame(
            check.names = F,
            `Node_number` = length(igraph::V(go)), # number of nodes
            `Edge_number` = length(igraph::E(go)), # number of edges
            `Edge_density` = igraph::edge_density(go), # density of network, connectance
            `Negative_percentage` = ifelse(!is.null(E(go)$cor), sum(igraph::E(go)$cor < 0) / length(igraph::E(go)), NA), # negative edges percentage
            `Average_path_length` = igraph::average.path.length(up), # Average path length

            `Global_efficiency` = igraph::global_efficiency(up),
            `Average_degree` = mean(igraph::degree(go)), # Average degree
            `Average_weighted_degree` = ifelse(is.null(igraph::E(go)$weight), mean(igraph::degree(go)), sum(igraph::E(go)$weight) / length(igraph::V(go))), # weighted degree
            Diameter = igraph::diameter(up), # network diameter
            `Clustering_coefficent` = igraph::transitivity(go), # Clustering coefficient
            `Centralized_betweenness` = igraph::centralization.betweenness(go)$centralization, # Betweenness centralization
            `Natural_connectivity` = nc(go) # natural
        )

        if (!fast) {
            # mean_dist=mean_distance(go)#
            # w_mean_dist=ifelse(is.null(E(go)$weight),mean_dist,mean_distance(go))
            # v_conn= vertex.connectivity(go) #
            # e_conn= edge.connectivity(go) #
            # components= count_components(go) #
            modularity <- igraph::modularity(igraph::cluster_fast_greedy(go)) #
            rand.g <- igraph::erdos.renyi.game(length(V(go)), length(E(go)), type = "gnm")
            rand_m <- igraph::modularity(igraph::cluster_fast_greedy(rand.g))
            relative_modularity <- (modularity - rand_m) / rand_m #

            n_index <- data.frame(
                check.names = F,
                n_index,
                Modularity = modularity,
                `Relative_modularity` = relative_modularity,
                `Centralized_closeness` = igraph::centralization.closeness(go)$centralization, # Closeness centralization
                `Centralized_degree` = igraph::centralization.degree(go)$centralization, # Degree centralization
                `Centralized_eigenvector` = igraph::centralization.evcent(go)$centralization # eigenvector centralization
            )
        }
        n_index <- apply(n_index, 1, FUN = \(x)replace(x, is.nan(x), 0)) %>%
            t() %>%
            as.data.frame()
        n_index <- cbind_new(get_n(go, simple = TRUE), n_index)
    }
    if ("v" %in% mode) {
        # Calculate Vertices Parameters
        v_index <- data.frame(
            check.names = F,
            Degree = igraph::degree(go),
            `Clustering_coefficent` = igraph::transitivity(go, type = "local"), # local clustering coefficient
            Betweenness = igraph::betweenness(go), # betweenness
            Eccentricity = igraph::eccentricity(go),
            Closeness = igraph::closeness(go),
            `Hub_score` = igraph::hub_score(go)[["vector"]]
            # page_rank = page.rank(go)$vector
            # igraph::evcent(go)[["vector"]]
            # igraph::local_efficiency(go)
        )

        # weighted degree
        if (!is.null(E(go)$cor)) {
            get_e(go) -> edge_list
            edge_list %>%
                dplyr::select(from, cor) %>%
                rbind(., dplyr::select(edge_list, to, cor) %>% dplyr::rename(from = to)) %>%
                dplyr::group_by(from) %>%
                dplyr::summarise(w_degree = sum(cor)) -> w_degree
            v_index$`Average_weighted_degree` <- w_degree[match(rownames(v_index), w_degree$from), "w_degree"] %>% unlist()
        }

        v_index <- apply(v_index, 1, FUN = \(x)replace(x, is.nan(x), 0)) %>%
            t() %>%
            as.data.frame()
        v_index <- cbind_new(get_v(go), v_index)
    }
    if ("e" %in% mode) {
        # Calculate Edges Parameters
        e_index <- get_e(go)
        # if(!(edge_attr(go)%>%unlist()%>%is.null()))e_index=data.frame(edge_attr(go),e_index)
    }
    return(list(n_index = n_index, v_index = v_index, e_index = e_index))
}


#' Add topological indexes for a network
#' @param go igraph or metanet
#' @param force replace exist net_par
#'
#' @export
#' @rdname net_par
c_net_index <- function(go, force = FALSE) {
    if (!force) {
        if (!is.null(graph_attr(go)[["net_par"]])) stop("Already calculate net_pars, set `force = TRUE to replace existed net_par")
    }
    net_par(go, fast = FALSE) -> res
    graph_attr(go) <- as.list(res$n_index)
    graph_attr(go)[["net_par"]] <- TRUE
    vertex_attr(go) <- as.list(res$v_index)
    edge_attr(go) <- as.list(res$e_index)
    go
}


#' Get skeleton network according to a group
#'
#' @param go network
#' @param Group vertex column name
#' @param count take which column count, default: NULL
#' @param top_N top_N
#'
#' @return skeleton network
#' @export
#'
#' @examples
#' get_group_skeleton(co_net) -> ske_net
#' skeleton_plot(ske_net)
get_group_skeleton <- function(go, Group = "v_class", count = NULL, top_N = 8) {
    lib_ps("igraph")
    stopifnot(is_igraph(go))
    direct <- is_directed(go)

    if (!Group %in% vertex_attr_names(go)) stop("no Group named ", Group, " !")
    get_v(go) -> tmp_v
    tmp_v %>% dplyr::select(name, !!Group) -> nodeGroup
    colnames(nodeGroup) <- c("name", "Group")
    nodeGroup$Group <- as.factor(nodeGroup$Group)
    # summary edges counts in each e_type
    suppressMessages(anno_edge(go, nodeGroup) %>% get_e() -> edge)
    {
        if (is.null(count)) {
            edge$count <- 1
        } else {
            edge$count <- edge[, count]
        }
    }
    bb <- data.frame()
    for (i in unique(edge$e_type)) {
        tmp <- edge[edge$e_type == i, c("Group_from", "Group_to", "count")]
        tmp <- dplyr::mutate_if(tmp, is.factor, as.character)
        # tmp=pcutils:::gettop(tmp,top_N)
        bb <- rbind(bb, data.frame(summ_2col(tmp,
            direct = direct
        ), e_type = i))
    }
    tmp_go <- igraph::graph_from_data_frame(bb, directed = direct)
    nodeGroup <- cbind_new(nodeGroup, data.frame(v_group = tmp_v$v_group))

    # nodeGroup=mutate_all(nodeGroup,as.character)
    # nodeGroup=rbind(nodeGroup,c("others","others","others"))

    dplyr::distinct(nodeGroup, Group, v_group) %>% tibble::column_to_rownames("Group") -> v_group_tab

    V(tmp_go)$v_group <- v_group_tab[V(tmp_go)$name, "v_group"]
    V(tmp_go)$v_class <- V(tmp_go)$name
    V(tmp_go)$size <- stats::aggregate(tmp_v$size, by = list(tmp_v[, Group]), sum)[["x"]]
    suppressWarnings({
        V(tmp_go)$count <- tmp_v %>%
            dplyr::group_by_(Group) %>%
            dplyr::count() %>%
            dplyr::pull(n)
    })

    tmp_go <- c_net_update(tmp_go)
    get_e(tmp_go) -> tmp_e

    E(tmp_go)$width <- E(tmp_go)$label <- tmp_e$count

    graph.attributes(tmp_go)$n_type <- "skeleton"
    graph.attributes(tmp_go)$skeleton <- Group
    tmp_go
}

#' Skeleton plot
#'
#' @param ske_net skeleton
#' @param ... additional parameters for \code{\link[igraph]{igraph.plotting}}
#' @export
#' @rdname get_group_skeleton
skeleton_plot <- function(ske_net, ...) {
    flag <- TRUE
    tmp_go <- ske_net
    lib_ps("igraph", library = FALSE)
    if (get_n(tmp_go)$n_type != "skeleton") stop("Not a skeleton network")
    get_e(tmp_go) -> tmp_e
    get_v(tmp_go) -> tmp_v
    # c_net_lay(tmp_go)->coors

    # some
    params <- list(...)
    params_name <- names(params)
    legend_position <- params[["legend_position"]]
    legend_position_default <- c(left_leg_x = -1.9, left_leg_y = 1, right_leg_x = 1.2, right_leg_y = 1)
    if (!"legend_position" %in% params_name) {
        legend_position <- legend_position_default
    } else {
        if (is.null(names(legend_position))) legend_position <- setNames(legend_position, names(legend_position_default)[seq_along(legend_position)])
        legend_position <- pcutils::update_param(legend_position_default, legend_position)
    }
    if ("vertex.color" %in% params_name) tmp_v$color <- condance(data.frame(tmp_v$color, tidai(tmp_v$v_class, params[["vertex.color"]])))
    if ("legend_cex" %in% params_name) {
        legend_cex <- parms[["legend_cex"]]
    } else {
        legend_cex <- 1
    }

    for (i in unique(tmp_e$e_type)) {
        left_leg_x <- legend_position["left_leg_x"]
        left_leg_y <- legend_position["left_leg_y"]

        # main plot
        tmp_go1 <- c_net_filter(tmp_go, e_type == i, mode = "e")
        c_net_plot(tmp_go1,
            vertex.size = pcutils::mmscale(V(tmp_go1)$size, 10, 16),
            edge.width = pcutils::mmscale(E(tmp_go1)$width, 0.5, 8), col_legend = FALSE, ...,
            lty_legend = FALSE, legend_number = FALSE, size_legend = FALSE
        )

        if ("legend" %in% params_name) {
            if (!params[["legend"]]) flag <- FALSE
        }
        if (flag) {
            # color legend
            pchls <- c("circle" = 21, "square" = 22)
            for (i in 1:length(unique(tmp_v$v_group))) {
                g_i <- unique(tmp_v$v_group)[i]
                tmp_v1 <- tmp_v[tmp_v$v_group == g_i, c("v_class", "count", "color", "shape")]
                if (TRUE) {
                    le_text <- paste(tmp_v1$v_class, tmp_v1$count, sep = ": ")
                } else {
                    le_text <- unique(tmp_v1$v_class)
                }
                legend(left_leg_x, left_leg_y,
                    cex = 0.7 * legend_cex, adj = 0,
                    legend = le_text, title.cex = 0.8 * legend_cex,
                    title = g_i, title.font = 2, title.adj = 0,
                    col = "black", pt.bg = unique(tmp_v1$color), bty = "n", pch = pchls[unique(tmp_v1$shape)]
                )
                left_leg_y <- left_leg_y - (length(unique(tmp_v1$v_class)) * 0.12 + 0.2) * legend_cex
            }
        }
    }
}

#' Link summary of the network
#'
#' @param go igraph
#' @param group summary which group of vertex attribution in names(vertex_attr(go))
#' @param e_type "positive", "negative", "all"
#' @param topN topN of group, default:5
#' @param colors colors
#' @param legend_number legend with numbers
#' @param legend all legends
#' @param legend_position legend_position, default: c(left_leg_x=-1.9,left_leg_y=1,right_leg_x=1.2,right_leg_y=1)
#' @param legend_cex 	character expansion factor relative to current par("cex"), default: 1
#' @param col_legend_order col_legend_order vector,
#' @param group_legend_title group_legend_title, length must same to the numbers of v_group
#' @param group_legend_order group_legend_order vector
#'
#' @return plot
#' @export
#'
#' @examples
#' links_stat(co_net, topN = 10)
#' modu_dect(co_net) -> co_net_modu
#' links_stat(co_net_modu, group = "module")
links_stat <- function(go, group = "v_class", e_type = "all", topN = 6, colors = NULL,
                       legend_number = FALSE, legend = TRUE, legend_cex = 1,
                       legend_position = c(left_leg_x = -1.6, left_leg_y = 1, right_leg_x = 1.2, right_leg_y = 1),
                       col_legend_order = NULL,
                       group_legend_title = NULL, group_legend_order = NULL) {
    pcutils::lib_ps("igraph", "dplyr", library = FALSE)
    direct <- is_directed(go)
    go <- c_net_set(go, vertex_class = group)

    get_v(go) -> v_index
    v_index %>% dplyr::select("name", "v_class") -> map

    suppressMessages(anno_edge(go, map) %>% get_e() -> edge)
    # statistics
    if (e_type != "all") edge %>% dplyr::filter(e_type == !!e_type) -> edge
    summ_2col(edge[, paste0("v_class", c("_from", "_to"))], direct = direct) -> bb
    colnames(bb) <- c("from", "to", "count")

    dplyr::group_by(bb, from) %>%
        dplyr::summarise(n = sum(count)) %>%
        dplyr::arrange(-n) %>%
        dplyr::top_n(topN, n) %>%
        dplyr::pull(from) -> nnn

    get_v(go) -> tmp_v
    if (is.null(colors)) {
        colors <- setNames(unique(tmp_v$color), unique(tmp_v$v_class))
    }
    # plot
    bb2 <- mutate(bb, from = ifelse(from %in% nnn, from, "Others"), to = ifelse(to %in% nnn, to, "Others")) %>% summ_2col(direct = direct)
    pcutils::my_circo(bb2, reorder = FALSE, pal = colors)
}

# 每个分组可以构建一个网络，每个网络都可以用link_stat得到一些互作的数量（互作强度），可以再看这些数量和分组间某些指标的相关性。


#' Fit power-law distribution for an igraph
#'
#' @param go igraph
#' @param p.value calculate p.value
#'
#' @return ggplot
#' @export
#'
#' @examples
#' fit_power(co_net)
fit_power <- function(go, p.value = FALSE) {
    # igraph::degree distribution
    degree_dist <- table(igraph::degree(go))
    dat <- data.frame(degree = as.numeric(names(degree_dist)), count = as.numeric(degree_dist))
    # fit, set the original a & b
    mod <- stats::nls(count ~ a * degree^b, data = dat, start = list(a = 2, b = 1.5))
    summary(mod)
    # extract the coefficient
    a <- round(coef(mod)[1], 3)
    b <- round(coef(mod)[2], 3)
    fit <- fitted(mod)
    SSre <- sum((dat$count - fit)^2)
    SStot <- sum((dat$count - mean(dat$count))^2)
    R2 <- round(1 - SSre / SStot, 3)

    # bootstrap t get p.value
    if (p.value) {
        dat_rand <- dat
        p_num <- lapply(seq_len(999), \(i){
            dat_rand$count <- sample(dat_rand$count)
            SSre_rand <- sum((dat_rand$count - fit)^2)
            SStot_rand <- sum((dat_rand$count - mean(dat_rand$count))^2)
            R2_rand <- 1 - SSre_rand / SStot_rand
            R2_rand > R2
        })
        p_value <- (sum(p_num) + 1) / (999 + 1)
    }

    p <- ggplot(dat, aes(x = degree, y = count)) +
        geom_point() +
        theme_bw() +
        stat_smooth(method = "nls", formula = y ~ a * x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
        labs(x = "Degree", y = "Count")

    if (p.value) {
        label <- data.frame(
            x = 0.8 * max(dat$degree),
            y = c(0.9, 0.85, 0.8) * max(dat$count),
            formula = c(
                sprintf("italic(Y) == %.3f*italic(X)^%.3f", a, b),
                sprintf("italic(R^2) == %.3f", R2),
                sprintf("italic(P) < %.3f", p_value)
            )
        )
    } else {
        label <- data.frame(
            x = 0.8 * max(dat$degree),
            y = c(0.9, 0.85) * max(dat$count),
            formula = c(
                sprintf("italic(Y) == %.3f*italic(X)^%.3f", a, b),
                sprintf("italic(R^2) == %.3f", R2)
            )
        )
    }
    p + geom_text(aes(x = x, y = y, label = formula), data = label, parse = TRUE)
}


#' Degree distribution comparison with random network
#'
#' @param go igraph object
#'
#' @return ggplot
#' @export
#'
#' @examples
#' rand_net(co_net)
rand_net <- function(go = go) {
    lib_ps("igraph")
    # generate a random network
    rand.g <- igraph::erdos.renyi.game(length(V(go)), length(E(go)), type = "gnm")

    data1 <- data.frame(
        freq = igraph::degree_distribution(go), net = "Network",
        degree = 0:(length(degree_distribution(go)) - 1)
    )

    data2 <- data.frame(
        freq = igraph::degree_distribution(rand.g), net = "Random E-R",
        degree = 0:(length(degree_distribution(rand.g)) - 1)
    )

    # if data1[1,1]=0, it'is delete single vertex
    if (data1[1, 1] == 0) data1 <- data1[-1, ]

    data <- rbind(data1, data2)
    p1 <- ggplot(data) +
        geom_point(aes(x = degree, y = freq, group = net, fill = net), pch = 21, size = 2) +
        geom_smooth(aes(x = degree, y = freq, group = net, color = net), se = FALSE, method = "loess", formula = "y ~ x") +
        labs(x = "Degree", y = "Proportion") +
        scale_color_manual(values = c("#F58B8B", "#7AADF0")) +
        scale_fill_manual(values = c("#F58B8B", "#7AADF0")) +
        MetaNet_theme +
        theme(legend.position = c(0.8, 0.9), legend.title = element_blank())
    print(p1)
    return(rand.g)
}

#' Net_pars of many random network
#'
#' @param go igraph
#' @param reps simulation time
#' @param threads threads
#' @param verbose verbose
#'
#' @export
#' @rdname compare_rand
rand_net_par <- function(go, reps = 99, threads = 1, verbose = TRUE) {
    # parallel
    # main function
    loop <- function(i) {
        # generate a random network
        rand.g <- igraph::erdos.renyi.game(length(igraph::V(go)),
            length(igraph::E(go)),
            type = "gnm"
        )
        indexs <- net_par(rand.g, mode = "n")[["n_index"]]
        wc <- igraph::cluster_fast_greedy(rand.g)
        indexs$modularity <- igraph::modularity(wc)
        indexs
    }
    {
        if (threads > 1) {
            pcutils::lib_ps("foreach", "doSNOW", "snow")
            if (verbose) {
                pb <- utils::txtProgressBar(max = reps, style = 3)
                opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
            } else {
                opts <- NULL
            }
            cl <- snow::makeCluster(threads)
            doSNOW::registerDoSNOW(cl)
            res <- foreach::foreach(i = 1:reps, .options.snow = opts) %dopar% {
                loop(i)
            }
            snow::stopCluster(cl)
            gc()
            pcutils::del_ps("doSNOW", "snow", "foreach")
        } else {
            res <- lapply(1:reps, loop)
        }
    }
    # simplify method
    rand_net_pars <- do.call(rbind, res)

    rand_net_pars
}

#' Compare some indexes between your net with random networks
#'
#' @param pars your net pars resulted by net_pars()
#' @param randp random networks pars resulted by rand_net_par()
#' @param index compared indexes: "Average_path_length","Clustering_coefficent" or else
#'
#' @return ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' data("c_net")
#' rand_net_par(co_net_rmt, reps = 30) -> randp
#' net_par(co_net_rmt, fast = FALSE) -> pars
#' compare_rand(pars, randp)
#' }
compare_rand <- function(pars, randp, index = c("Average_path_length", "Clustering_coefficent")) {
    labss <- t(pars$n_index[, index, drop = FALSE]) %>% as.data.frame()
    rownames(labss) -> labss$indexes
    pcutils::group_box(randp[, index, drop = FALSE]) +
        geom_hline(data = labss, aes(yintercept = V1), linetype = 2) +
        geom_text(data = labss, aes(x = 1, y = V1 * 1.05, label = paste0("Net: ", round(V1, 3)))) +
        MetaNet_theme +
        theme(legend.position = "none", axis.text.x = element_blank())
}


#' Calculate small-world coefficient
#'
#' @param go igraph or metanet
#' @param reps simulation time
#' @param threads threads
#' @param verbose verbose
#'
#' @export
#' @examples
#' \dontrun{
#' smallworldness(co_net)
#' }
smallworldness <- function(go, reps = 99, threads = 1, verbose = TRUE) {
    rand_net_par(go, reps = reps, threads = threads, verbose = verbose) -> rands
    small_world_coefficient <- (igraph::transitivity(go) / mean(rands$Clustering_coefficent)) /
        (igraph::average.path.length(go) / mean(rands$`Average_path_length`))
    small_world_coefficient
}
