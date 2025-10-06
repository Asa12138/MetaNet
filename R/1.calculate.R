# =======1.calculate========
#' Calculate correlation for one or two t(otutab), or distance for one t(otutab).
#'
#' @param totu t(otutab), row are samples, column are features.
#' @param totu2 t(otutab2) or NULL, row are samples, column are features.
#' @param method "spearman" (default), "pearson", "sparcc", or distance index from \code{\link[vegan]{vegdist}}.
#' @param filename the prefix of saved .corr file or FALSE.
#' @param p.adjust.method see \code{\link[stats]{p.adjust}}
#' @param p.adjust.mode see \code{\link{p.adjust.table}}
#' @param threads threads, default: 1.
#' @param verbose verbose, default: TRUE.
#'
#' @return a corr object with 3 elements:
#' \item{r}{default: spearman correlation}
#' \item{p.value}{default: p-value of spearman correlation}
#' \item{p.adjust}{default p.adjust.method = NULL}
#' @family calculate
#' @export
#' @aliases c_net_cal
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab) -> totu
#' c_net_calculate(totu) -> corr
#' metadata[, 3:10] -> env
#' c_net_calculate(totu, env) -> corr2
c_net_calculate <- function(totu, totu2 = NULL, method = "spearman", filename = FALSE,
                            p.adjust.method = NULL, p.adjust.mode = "all", threads = 1, verbose = TRUE) {
  if (!is.null(totu2)) {
    tls <- check_tabs(totu, totu2)
    totu <- tls[[1]]
    totu2 <- tls[[2]]
  } else {
    totu <- check_tabs(totu)
  }

  if (method %in% c("spearman", "pearson")) {
    corr <- fast_cor(totu, totu2, method = method)
  } else if (method == "sparcc") {
    if (!is.null(totu2)) warning("sparcc only take the totu, ignore the totu2.")
    message("sparcc is not supported on CRAN version. Please install `zdk123/SpiecEasi` from github and use the function `sparcc`.")
    # corr <- par_sparcc(totu, threads = threads, verbose = verbose)
    corr <- NULL
  } else if (method %in% c(
    "manhattan", "euclidean", "canberra", "bray",
    "kulczynski", "gower", "morisita", "horn", "mountford",
    "jaccard", "raup", "binomial", "chao", "altGower", "cao",
    "mahalanobis", "clark", "chisq", "chord", "hellinger",
    "aitchison", "robust.aitchison"
  )) {
    corr <- cal_sim(totu, totu2, method = method)
  } else {
    stop("method should be one of 'spearman', 'pearson', 'sparcc', or distance index from vegan.")
  }

  if (!is.null(p.adjust.method)) {
    p.adjust <- p.adjust.table(corr$p.value, p.adjust.method)
    res <- list(r = corr$r, p.value = corr$p.value, p.adjust = p.adjust)
  } else {
    res <- list(r = corr$r, p.value = corr$p.value)
  }

  class(res) <- "corr"
  attributes(res)$method <- method

  # save the correlation result
  if (is.logical(filename)) {
    if (filename) filename <- paste0("c_net_", date())
  }
  if (is.character(filename)) {
    save_corr(res, filename)
  }
  return(res)
}

#' Check tables and extract common samples
#'
#' @param ... tables
#' @return formatted tables
#' @export
#' @examples
#' data("otutab", package = "pcutils")
#' check_tabs(otutab)
check_tabs <- function(...) {
  tables <- list(...)
  if (all(class(tables[[1]]) == "list")) tables <- tables[[1]]

  names(tables) <- NULL

  if (length(tables) > 1) {
    comm <- Reduce(intersect, lapply(tables, rownames))
    if (length(comm) < 2) stop("There are ", length(comm), " common sample! Can not calculate correlation.")
    if (all(lapply(tables, \(i)identical(rownames(i), comm)) %>% unlist())) {
      message("All samples matched.")
    } else {
      message("Extract ", length(comm), " commmon samples.")
    }
    dup <- lapply(tables, colnames) %>% do.call(c, .)
    dup <- dup[duplicated(dup)]
    if (length(dup) > 0) {
      stop("Duplicated colnames found: ", paste0(dup, collapse = ", "), "\nPlease check colnames of input tables.")
    } else {
      message("All features are OK.")
    }
    tables <- lapply(tables, \(i)i[comm, ])
  }

  if (length(tables) == 1) tables <- tables[[1]]
  return(tables)
}


#' Save a corr object
#'
#' @param corr a corr object
#' @param filename filename without extension, default: "corr"
#'
#' @return a .corr file
#' @export
#'
save_corr <- function(corr, filename = "corr") {
  stopifnot(inherits(corr, "corr"))
  if (!grepl("\\.corr$", filename)) filename <- paste0(filename, ".corr")
  if (t_flag(corr$r)) {
    # 节约一半的储存空间
    corr$r[upper.tri(corr$r)] -> r
    corr$p.value[upper.tri(corr$p.value)] -> p.value
    if (!is.null(corr$p.adjust)) {
      corr$p.adjust[upper.tri(corr$p.adjust)] -> p.adj
      if (all(p.value == p.adj)) {
        p.adj <- FALSE
      }
    } else {
      p.adj <- NULL
    }

    corr_names <- rownames(corr$r)
    saveRDS(
      list(
        r = r, p.value = p.value, p.adjust = p.adj,
        corr_names = corr_names, corr_attr = attributes(corr)
      ),
      file = filename
    )
  } else {
    saveRDS(corr, file = filename)
  }
}

#' Read a corr object
#'
#' @param filename filename of .corr
#'
#' @return a corr object
#' @export
#' @family calculate
read_corr <- function(filename) {
  # r <- read.csv(paste0(filename, "_r.csv"), row.names = 1, check.names = FALSE)
  # p.value <- read.csv(paste0(filename, "_p.csv"), row.names = 1, check.names = FALSE)
  # p.adjust <- read.csv(paste0(filename, "_p_adj.csv"), row.names = 1, check.names = FALSE)
  if (!grepl("\\.corr$", filename)) filename <- paste0(filename, ".corr")
  in_corr <- readRDS(filename)
  if ("corr_names" %in% names(in_corr)) {
    corr_names <- in_corr$corr_names
    r <- in_corr$r
    p.value <- in_corr$p.value
    p.adj <- in_corr$p.adj
    corr_attr <- in_corr$corr_attr

    new_p.adjust <- new_p <- new_r <- matrix(0, nrow = length(corr_names), ncol = length(corr_names), dimnames = list(corr_names, corr_names))
    new_r[upper.tri(new_r)] <- r
    new_r <- new_r + t(new_r)
    diag(new_r) <- 1
    new_p[upper.tri(new_p)] <- p.value
    new_p <- new_p + t(new_p)
    if (is.null(p.adj)) {
      new_p.adjust <- NULL
    } else if (is.logical(p.adj)) {
      new_p.adjust <- new_p
    } else {
      new_p.adjust[upper.tri(new_p.adjust)] <- p.adj
      new_p.adjust <- new_p.adjust + t(new_p.adjust)
      diag(new_p.adjust) <- 1
    }
    in_corr <- list()
    in_corr$r <- new_r
    in_corr$p.value <- new_p
    if (!is.null(new_p.adjust)) in_corr$p.adjust <- new_p.adjust
    attributes(in_corr) <- corr_attr
  }
  return(in_corr)
}

#' Plot a correlation heatmap
#'
#' @param corr a corr object, must contain 'r' and 'p.value' matrices.
#' @param show.p whether to show p-values as significance stars on the heatmap.
#' @param cluster_rows whether to cluster rows.
#' @param cluster_cols whether to cluster columns.
#' @param ... additional arguments passed to `pheatmap::pheatmap`.
#' @return plot of the correlation heatmap.
#' @export
plot_corr_heatmap <- function(corr,
                              show.p = FALSE,
                              cluster_rows = TRUE,
                              cluster_cols = TRUE,
                              ...) {
  # 检查必要包
  lib_ps(pheatmap, library = FALSE)

  # 检查输入结构
  if (!inherits(corr, "corr") ||
    !all(c("r", "p.value") %in% names(corr))) {
    stop("`corr` must be a corr object with 'r' and 'p.value' matrices.")
  }
  p.mat <- corr$p.value

  # 准备显著性标记（可选）
  if (show.p) {
    # 定义显著性符号
    signif_codes <- cut(
      p.mat,
      breaks = c(0, 0.001, 0.01, 0.05, 1),
      labels = c("***", "**", "*", ""),
      include.lowest = TRUE
    )

    # 将相关系数矩阵转换为带标记的字符矩阵
    r_with_stars <- matrix(
      # sprintf("%.2f%s", corr$r, signif_codes),
      signif_codes,
      nrow = nrow(corr$r),
      dimnames = dimnames(corr$r)
    )

    # 移除多余的+号（如0.10* -> 0.1*）
    r_with_stars <- gsub("0\\.(\\d)", ".\\1", r_with_stars)
  }

  # 绘制热图
  pheatmap::pheatmap(
    mat = corr$r,
    display_numbers = if (show.p) r_with_stars else FALSE,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    number_color = "black",
    border_color = if (nrow(corr$r) < 20) "gray60" else NA,
    fontsize_number = 10,
    ...
  )
}

#' Fast correlation calculation
#'
#' @param totu t(otutab), row are samples, column are features.
#' @param totu2 t(otutab) or NULL, row are samples, column are features.
#' @param method "spearman" or "pearson"
#'
#' @export
#'
#' @return a list with 2 elements:
#' \item{r}{default: spearman correlation}
#' \item{p.value}{default: p-value of spearman correlation}
#' @family calculate
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab[1:100, ]) -> totu
#' fast_cor(totu, method = "spearman") -> corr
fast_cor <- function(totu, totu2 = NULL, method = c("pearson", "spearman")) {
  method <- match.arg(method, c("pearson", "spearman"))
  totu <- as.matrix(totu)
  if (!is.null(totu2)) {
    totu2 <- as.matrix(totu2)
    r <- stats::cor(totu, totu2, method = method)
    df <- dim(totu)[1] - 2
    t <- r * sqrt(df / (1 - r^2))
    p <- -2 * expm1(stats::pt(abs(t), df, log.p = TRUE))
  } else {
    r <- stats::cor(totu, method = method)
    p <- r
    p[p != 0] <- 0
    r_tri <- r[upper.tri(r)]
    df <- dim(totu)[1] - 2
    t <- r_tri * sqrt(df / (1 - r_tri^2))
    p[upper.tri(p)] <- -2 * expm1(stats::pt(abs(t), df, log.p = TRUE))
    p <- p + t(p)
  }
  return(list(r = r, p.value = p))
}

Hmisc_cor <- function(totu, totu2 = NULL, method = c("spearman", "pearson")[1]) {
  lib_ps("Hmisc", library = FALSE)
  totu <- as.matrix(totu)
  if (!is.null(totu2)) {
    totu2 <- as.matrix(totu2)
    tmp <- Hmisc::rcorr(totu, totu2, type = method)
    r <- tmp$r[1:ncol(totu), (ncol(totu) + 1):ncol(tmp$r)]
    p <- tmp$P[1:ncol(totu), (ncol(totu) + 1):ncol(tmp$r)]
    return(list(r = r, p.value = p))
  }
  tmp <- Hmisc::rcorr(totu, type = method)
  p <- tmp$P
  p[is.na(p)] <- 0
  return(list(r = tmp$r, p.value = p))
}

#' Calculate similarity for one t(otutab)
#'
#' @param totu t(otutab), row are samples, column are features.
#' @param method Dissimilarity index, see \code{\link[vegan]{vegdist}}.
#' @param totu2 t(otutab) or NULL, row are samples, column are features.
#'
#' @family calculate
#' @return similarity = 1-distance
#' @export
#' @seealso \code{\link[vegan]{vegdist}}
#' @examples
#' if (requireNamespace("vegan")) {
#'   data("otutab", package = "pcutils")
#'   t(otutab) -> totu
#'   cal_sim(totu) -> sim_corr
#' }
cal_sim <- function(totu, totu2 = NULL, method = "bray") {
  lib_ps("vegan", library = FALSE)
  if (is.null(totu2)) {
    vegan::vegdist(t(totu), method = method) %>% as.matrix() -> dist
  } else {
    n1 <- ncol(totu)
    n2 <- ncol(totu2)
    dist <- matrix(NA, nrow = n1, ncol = n2, dimnames = list(colnames(totu), colnames(totu2)))
    for (i in seq_len(n1)) {
      for (j in seq_len(n2)) {
        dist[i, j] <- vegan::vegdist(rbind(totu[, i], totu2[, j]), method = method)
      }
    }
  }

  sim <- 1 - dist
  p <- dist
  message("p-value is not supported for distance index, all set as 0.")
  p[p != 0] <- 0

  return(list(r = sim, p.value = p))
}

cal_KLD <- function(totu, totu2 = NULL, method = "KLD") {
  # Kullback-Leibler divergence
  # KLD = function(p, q) {
  #   sum(p * log(p / q))
  # }
  # p = c(0.1, 0.2, 0.3, 0.4)
  # q = c(0.2, 0.3, 0.2, 0.3)
  # KLD(p, q)
  lib_ps("philentropy", library = FALSE)
  dat <- t(totu) / rowSums(t(totu))
  philentropy::KL(dat, unit = "log")
}

#' p.adjust apply on a correlation table (matrix or data.frame)
#'
#' @param pp table of p-values
#' @param method see \code{\link[stats]{p.adjust}}, default: "BH".
#' @param mode "all" for all values; "rows" adjust each row one by one; "columns" adjust each column one by one. Default: "all".
#'
#' @return a table of adjusted p-values
#' @export
#'
#' @family calculate
#' @examples
#' matrix(abs(rnorm(100, 0.01, 0.1)), 10, 10) -> pp
#' p.adjust.table(pp, method = "BH", mode = "all") -> pp_adj
p.adjust.table <- \(pp, method = "BH", mode = "all"){
  mode <- match.arg(mode, c("all", "rows", "columns"))
  pp <- as.matrix(pp)
  if (mode == "all") {
    if (t_flag(pp)) {
      lp <- lower.tri(pp)
      pa <- pp[lp]
      pa <- p.adjust(pa, method)
      pp[lower.tri(pp, diag = FALSE)] <- pa
      pp[upper.tri(pp, diag = FALSE)] <- 0
      pp <- pp + t(pp)
    } else {
      pp1 <- p.adjust(pp, method)
      pp <- matrix(pp1, nrow(pp), ncol(pp), dimnames = list(rownames(pp), colnames(pp)))
    }
  } else if (mode == "rows") {
    pp <- t(apply(pp, 1, p.adjust, method = method))
  } else if (mode == "columns") {
    pp <- apply(t(pp), 1, p.adjust, method = method)
  }
  return(pp)
}
