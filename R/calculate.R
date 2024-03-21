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
  } else if (method %in% c(
    "manhattan", "euclidean", "canberra", "bray",
    "kulczynski", "gower", "morisita", "horn", "mountford",
    "jaccard", "raup", "binomial", "chao", "altGower", "cao",
    "mahalanobis", "clark", "chisq", "chord", "hellinger",
    "aitchison", "robust.aitchison"
  )) {
    if (!is.null(totu2)) warning("distance only take the totu, ignore the totu2.")
    corr <- cal_sim(totu, method = method)
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

  # save the correlation result
  if (is.logical(filename)) {
    if (filename) filename <- paste0("c_net_", date())
  }
  if (is.character(filename)) {
    saveRDS(res, file = paste0(filename, ".corr"))
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
      message("Extract, ", length(comm), " commmon samples.")
    }
    dup <- lapply(tables, colnames) %>% do.call(c, .)
    dup <- dup[duplicated(dup)]
    if (length(dup) > 0) {
      stop("Duplicated colnames found: ", paste0(dup, collapse = ", "), "\nPlease check colnames of input tables.")
    } else {
      message("All objects are OK.")
    }
    tables <- lapply(tables, \(i)i[comm, ])
  }

  if (length(tables) == 1) tables <- tables[[1]]
  return(tables)
}

#' Import corr from .csv file
#'
#' @param filename filename of .corr
#'
#' @return a corr object
#' @export
#' @family calculate
input_corr <- function(filename) {
  # r <- read.csv(paste0(filename, "_r.csv"), row.names = 1, check.names = FALSE)
  # p.value <- read.csv(paste0(filename, "_p.csv"), row.names = 1, check.names = FALSE)
  # p.adjust <- read.csv(paste0(filename, "_p_adj.csv"), row.names = 1, check.names = FALSE)
  if (!grepl(".corr", filename)) filename <- paste0(filename, ".corr")
  return(readRDS(filename))
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
#' @param norm hellinger normalization in features (default: FALSE).
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
cal_sim <- function(totu, method = "bray", norm = FALSE) {
  otu <- t(totu)
  lib_ps("vegan", library = FALSE)

  if (norm) otu <- vegan::decostand(otu, "hellinger", 1) %>% as.data.frame()
  vegan::vegdist(otu, method = method) %>% as.matrix() -> dist

  sim <- 1 - dist
  p <- dist
  p[p != 0] <- 0

  return(list(r = sim, p.value = p))
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
