# ========5.topological=======

#' Extract each sample network from the whole network
#'
#' @param whole_net the whole network
#' @param otutab otutab, columns are samples, these columns will be extract
#' @param threads threads, default: 1
#' @param save_net should save these sub_nets? FALSE or a filename
#' @param fast less indexes for faster calculate ?
#' @param verbose verbose
#' @param remove_negative remove negative edge or not? default: FALSE
#'
#' @return a dataframe contains all sub_net parameters
#' @export
#' @family topological
#' @examples
#' data(otutab, package = "pcutils")
#' extract_sample_net(co_net, otutab) -> sub_net_pars
extract_sample_net <- function(whole_net, otutab, threads = 1, save_net = FALSE, fast = TRUE, remove_negative = FALSE, verbose = TRUE) {
  i <- NULL
  V(whole_net)$name -> v_name
  reps <- ncol(otutab)

  if (verbose) message("extracting")
  sub_nets <- lapply(1:reps, \(i){
    rownames(otutab)[otutab[, i] > 0] -> exist_sp
    subgraph(whole_net, which(v_name %in% exist_sp)) -> spe_sub
    class(spe_sub) <- c("metanet", "igraph")
    return(spe_sub)
  })
  names(sub_nets) <- colnames(otutab)
  if (verbose) message("calculating topological indexes")

  # parallel
  # main function
  loop <- function(i) {
    spe_sub <- sub_nets[[i]]
    indexs <- net_par(spe_sub, mode = "n", fast = fast, remove_negative = remove_negative)[["n_index"]]
    wc <- igraph::cluster_fast_greedy(spe_sub, weights = abs(igraph::E(spe_sub)$weight))
    indexs$modularity <- igraph::modularity(wc)
    indexs
  }
  {
    if (threads > 1) {
      pcutils::lib_ps("foreach", "doSNOW", "snow", library = FALSE)
      if (verbose) {
        pb <- utils::txtProgressBar(max = reps, style = 3)
        opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
      } else {
        opts <- NULL
      }
      cl <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cl)
      res <- foreach::`%dopar%`(
        foreach::foreach(i = 1:reps, .options.snow = opts),
        loop(i)
      )
      snow::stopCluster(cl)
      gc()
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
#' @family topological
#' @examples
#' igraph::make_ring(10) %>% nc()
nc <- function(p) {
  if (length(p) == 0) {
    return(NaN)
  }
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


#' Calculate all topological indexes of a network
#'
#' @param go an igraph or metanet object
#' @param fast less indexes for faster calculate ?
#' @param mode calculate what? c("v", "e", "n", "all")
#' @param remove_negative remove negative edge or not? default: FALSE
#'
#' @return a 3-elements list
#' \item{n_index}{indexs of the whole network}
#' \item{v_index}{indexs of each vertex}
#' \item{e_index}{indexs of each edge}
#' @export
#' @family topological
#' @examples
#' igraph::make_graph("Walther") %>% net_par()
#' c_net_index(co_net) -> co_net_with_par
net_par <- function(go, mode = c("v", "e", "n", "all"), fast = TRUE, remove_negative = FALSE) {
  from <- to <- NULL
  stopifnot(is_igraph(go))
  if ("all" %in% mode) mode <- c("v", "e", "n")

  n_index <- NULL
  v_index <- NULL
  e_index <- NULL

  Negative_percentage <- ifelse(!is.null(E(go)$cor), sum(igraph::E(go)$cor < 0) / length(igraph::E(go)), NA)
  # remove negative weight
  if (remove_negative) {
    if (!is.null(E(go)$cor)) {
      # message("Remove negative correlation edges")
      c_net_filter(go, cor > 0, mode = "e") -> go
    }
  }

  # non-weighted network
  up <- go
  if (!is.null(igraph::edge_attr(up)[["weight"]])) {
    up <- igraph::delete_edge_attr(up, "weight")
  }

  if ("n" %in% mode) {
    # Calculate Network Parameters
    n_index <- data.frame(
      check.names = F,
      `Node_number` = length(igraph::V(go)), # number of nodes
      `Edge_number` = length(igraph::E(go)), # number of edges
      `Edge_density` = igraph::edge_density(go), # density of network, connectance
      `Negative_percentage` = Negative_percentage, # negative edges percentage
      `Average_path_length` = igraph::average.path.length(up), # Average path length

      `Global_efficiency` = igraph::global_efficiency(up),
      `Average_degree` = mean(igraph::degree(go)), # Average degree
      `Average_weighted_degree` = ifelse(is.null(igraph::E(go)$weight), mean(igraph::degree(go)), sum(igraph::E(go)$weight) / length(igraph::V(go))), # weighted degree
      Diameter = igraph::diameter(up), # network diameter
      `Clustering_coefficient` = igraph::transitivity(go), # Clustering coefficient
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
    # n_index <- cbind_new(get_n(go, simple = TRUE), n_index)
  }
  if ("v" %in% mode) {
    negative_weight <- FALSE
    if (!is.null(igraph::edge_attr(go)[["weight"]])) {
      if (any(igraph::edge_attr(go)[["weight"]] < 0)) {
        negative_weight <- TRUE
        warning("Weight vector must be positive, drop the weight.")
      }
    }
    # Calculate Vertices Parameters
    v_index <- data.frame(
      check.names = F,
      name = igraph::vertex_attr(go, "name"),
      Degree = igraph::degree(go),
      `Clustering_coefficient` = igraph::transitivity(go, type = "local"), # local clustering coefficient
      Betweenness = ifelse(negative_weight, igraph::betweenness(up), igraph::betweenness(go)), # betweenness
      Eccentricity = igraph::eccentricity(go),
      Closeness = ifelse(negative_weight, igraph::closeness(up), igraph::closeness(go)),
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
    # v_index <- cbind_new(get_v(go), v_index)
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
#' @param force replace existed net_par
#'
#' @export
#' @rdname net_par
c_net_index <- function(go, force = FALSE) {
  if (!force) {
    if (!is.null(graph_attr(go)[["net_par"]])) stop("Already calculated net_pars, set `force = TRUE to replace existed net_par")
  }
  net_par(go, fast = FALSE) -> res
  graph_attr(go) <- as.list(res$n_index)
  graph_attr(go)[["net_par"]] <- TRUE
  vertex_attr(go) <- as.list(res$v_index)
  edge_attr(go) <- as.list(res$e_index)
  go
}


#' Fit power-law distribution for an igraph
#'
#' @param go igraph
#' @param p.value calculate p.value
#'
#' @return ggplot
#' @export
#' @family topological
#' @examples
#' fit_power(co_net)
fit_power <- function(go, p.value = FALSE) {
  x <- y <- formula <- NULL
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
    p_value <- (sum(unlist(p_num)) + 1) / (999 + 1)
  }

  p <- ggplot(dat, aes(x = degree, y = count)) +
    geom_point() +
    theme_bw() +
    stat_smooth(method = "nls", formula = y ~ a * x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
    labs(x = "Degree", y = "Count")

  if (p.value) {
    label <- data.frame(
      x = 0.8 * max(dat$degree),
      y = c(0.9, 0.8, 0.7) * max(dat$count),
      formula = c(
        sprintf("italic(Y) == %.3f*italic(X)^%.3f", a, b),
        sprintf("italic(R^2) == %.3f", R2),
        sprintf("italic(P) < %.3f", p_value)
      )
    )
  } else {
    label <- data.frame(
      x = 0.8 * max(dat$degree),
      y = c(0.9, 0.8) * max(dat$count),
      formula = c(
        sprintf("italic(Y) == %.3f*italic(X)^%.3f", a, b),
        sprintf("italic(R^2) == %.3f", R2)
      )
    )
  }
  p + geom_text(aes(x = x, y = y, label = formula), data = label, parse = TRUE)
}


#' Plot degree distribution of networks
#'
#' @param gols list of metanet
#' @param net_names names of networks
#'
#' @returns ggplot
#' @export
plot_net_degree <- function(gols, net_names = NULL) {
  if (is.null(names(gols))) {
    names(gols) <- paste0("Network", seq_along(gols))
  }
  if (!is.null(net_names)) {
    names(gols) <- net_names
  }
  lapply(seq_along(gols), \(i){
    data.frame(
      freq = igraph::degree_distribution(gols[[i]]), net = names(gols)[[i]],
      degree = 0:(length(degree_distribution(gols[[i]])) - 1)
    )
  }) -> data_list
  data <- do.call(rbind, data_list)

  # if data1[1,1]=0, it'is delete single vertex
  if (data[1, 1] == 0) data <- data[-1, ]

  p1 <- ggplot(data) +
    geom_point(aes(x = degree, y = freq, group = net, fill = net), pch = 21, size = 2) +
    geom_smooth(aes(x = degree, y = freq, group = net, color = net), se = FALSE, method = "loess", formula = "y ~ x") +
    labs(x = "Degree", y = "Proportion") +
    MetaNet_theme +
    theme(legend.position = c(0.8, 0.9), legend.title = element_blank())
  p1
}

#' Degree distribution comparison with random network
#'
#' @param go igraph object
#' @param plot plot or not
#'
#' @return ggplot
#' @export
#' @family topological
#' @examples
#' rand_net(co_net)
rand_net <- function(go = go, plot = TRUE) {
  freq <- net <- NULL
  # generate a random network
  rand.g <- igraph::erdos.renyi.game(length(V(go)), length(E(go)), type = "gnm")

  if (!plot) {
    return(rand.g)
  }
  plot_net_degree(list(go, rand.g), net_names = c("Network", "Random E-R")) -> p1
  print(p1 + scale_color_manual(values = c("#F58B8B", "#7AADF0")) +
    scale_fill_manual(values = c("#F58B8B", "#7AADF0")))
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
  i <- NULL
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
      pcutils::lib_ps("foreach", "doSNOW", "snow", library = FALSE)
      if (verbose) {
        pb <- utils::txtProgressBar(max = reps, style = 3)
        opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
      } else {
        opts <- NULL
      }
      cl <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cl)
      res <- foreach::`%dopar%`(
        foreach::foreach(i = 1:reps, .options.snow = opts),
        loop(i)
      )
      snow::stopCluster(cl)
      gc()
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
#' @param index compared indexes: "Average_path_length","Clustering_coefficient" or else
#'
#' @return ggplot
#' @export
#' @family topological
#' @examples
#' data("c_net")
#' rand_net_par(co_net_rmt, reps = 30) -> randp
#' net_par(co_net_rmt, fast = FALSE) -> pars
#' compare_rand(pars, randp)
compare_rand <- function(pars, randp, index = c("Average_path_length", "Clustering_coefficient")) {
  V1 <- NULL
  labss <- t(pars$n_index[, index, drop = FALSE]) %>% as.data.frame()
  rownames(labss) -> labss$indexes

  p <- pcutils::group_box(randp[, index, drop = FALSE])

  p <- p +
    geom_hline(data = labss, aes(yintercept = V1), linetype = 2, color = "blue3") +
    geom_text(
      data = labss, aes(x = 1, y = V1 * 1.05, label = paste0("Network: ", round(V1, 3))),
      color = "blue3"
    ) +
    MetaNet_theme +
    theme(legend.position = "none", axis.text.x = element_blank())
  p
}


#' Calculate small-world coefficient
#'
#' @param go igraph or metanet
#' @param reps simulation time
#' @param threads threads
#' @param verbose verbose
#'
#' @return number
#' @export
#' @family topological
#' @examples
#' \donttest{
#' # set reps at least 99 when you run.
#' smallworldness(co_net, reps = 9)
#' }
smallworldness <- function(go, reps = 99, threads = 1, verbose = TRUE) {
  rand_net_par(go, reps = reps, threads = threads, verbose = verbose) -> rands
  small_world_coefficient <- (igraph::transitivity(go) / mean(rands$Clustering_coefficient)) /
    (igraph::average.path.length(go) / mean(rands$`Average_path_length`))
  small_world_coefficient
}
