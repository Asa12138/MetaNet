# =======1.calculate========
#' Calculate spearman correlation for one or two t(otutab)
#'
#' @param totu t(otutab)
#' @param totu2 t(otutab2) or NULL
#' @param method spearman, pearson, sparcc,
#' @param filename the prefix of saved .corr file or FALSE
#' @param p.adjust.method see \code{\link[stats]{p.adjust}}
#' @param threads threads, default: 1.
#' @param p.adjust.mode see \code{\link{p.adjust.table}}
#' @param verbose verbose
#'
#' @return a corr list with 3 elements:
#' \item{r}{default: spearman correlation}
#' \item{p.value}{default: p-value of spearman correlation}
#' \item{p.adjust}{default p.adjust.method = NULL}
#' @export
#'
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab) -> totu
#' c_net_cal(totu) -> corr
#' metadata[, 3:10] -> env
#' c_net_cal(totu, env) -> corr2
c_net_cal <- function(totu, totu2 = NULL, method = "spearman", filename = FALSE,
                      p.adjust.method = NULL, p.adjust.mode = "all", threads = 1, verbose = TRUE) {
    if (!is.null(totu2)) {
        tls <- check_tabs(totu, totu2)
        totu <- tls[[1]]
        totu2 <- tls[[2]]
    } else {
        totu <- check_tabs(totu)
    }

    if (method %in% c("spearman", "pearson")) {
        corr <- f_cor(totu, totu2, method = method)
    } else if (method == "sparcc") {
        if (!is.null(totu2)) warning("sparcc only take the totu, ignore the totu2.")
        # corr <- par_sparcc(totu, threads = threads, verbose = verbose)
    } else if (method %in% c(
        "manhattan", "euclidean", "canberra", "bray",
        "kulczynski", "gower", "morisita", "horn", "mountford",
        "jaccard", "raup", "binomial", "chao", "altGower", "cao",
        "mahalanobis", "clark", "chisq", "chord", "hellinger",
        "aitchison", "robust.aitchison"
    )) {
        if (!is.null(totu2)) warning("distance only take the totu, ignore the totu2.")
        corr <- par_sim(totu, method = method)
    }

    # if(threads>1){
    # corr<-par_cor(totu,totu2,threads = threads,method = "spearman")
    # }

    # if (is.null(totu2)) diag(corr$r) <- 0

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

#' Print method for 'corr' objects
#'
#' @param x 'corr' object
#' @param ... additional arguments
#'
#' @exportS3Method
#' @method print corr
#' @return No value
print.corr <- function(x, ...) {
    cat("Correlation table:\n")
    cat("Table dimensions:", nrow(x$r), "rows,", ncol(x$r), "columns\n")
}


#' Check tables
#'
#' @description
#' A short description...
#'
#' @param ... tables
#' @return format tables
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
#' @return a list with 7 elements:
#' \item{r}{default: spearman correlation}
#' \item{p.value}{default: p-value of spearman correlation}
#' \item{p.adjust}{default p.adjust.method = NULL}
#' @export
#'
input_corr <- function(filename) {
    # r <- read.csv(paste0(filename, "_r.csv"), row.names = 1, check.names = FALSE)
    # p.value <- read.csv(paste0(filename, "_p.csv"), row.names = 1, check.names = FALSE)
    # p.adjust <- read.csv(paste0(filename, "_p_adj.csv"), row.names = 1, check.names = FALSE)
    if (!grepl(".corr", filename)) filename <- paste0(filename, ".corr")
    return(readRDS(filename))
}


#' Fast correlation
#'
#' @param totu t(otutab)
#' @param totu2 t(otutab) or NULL
#' @param method spearman or pearson
#'
#' @export
#'
#' @return a list with 2 elements:
#' \item{r}{default: spearman correlation}
#' \item{p.value}{default: p-value of spearman correlation}
#'
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab[1:100, ]) -> totu
#' f_cor(totu, method = "spearman") -> corr
f_cor <- function(totu, totu2 = NULL, method = c("pearson", "spearman")) {
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


par_cor <- function(totu, totu2 = NULL, threads = 1, method = c("spearman", "pearson")) {
    i <- NULL
    method <- match.arg(method, c("pearson", "spearman"))

    if (method == "spearman") {
        totu <- apply(totu, 2, rank)
    } else if (method == "pearson") {
        totu <- totu
    } else {
        stop('method must be one of "spearman","pearson"')
    }
    r_p <- \(rx, ry){
        lxy <- sum((rx - mean(rx)) * (ry - mean(ry)))
        lxx <- sum((rx - mean(rx))^2)
        lyy <- sum((ry - mean(ry))^2)
        r <- lxy / sqrt(lxx * lyy)
        p <- 0
        n <- length(rx)
        t <- (r * sqrt(n - 2)) / sqrt(1 - r^2)
        p <- -2 * expm1(pt(abs(t), (n - 2), log.p = TRUE))
        return(c(r, p))
    }

    nc <- ncol(totu)

    pcutils::lib_ps("foreach", "doSNOW", "snow")
    pb <- utils::txtProgressBar(max = nc, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    cl <- snow::makeCluster(threads)
    doSNOW::registerDoSNOW(cl)
    if (is.null(totu2)) {
        corr <- foreach::`%dopar%`(foreach::foreach(i = 1:nc, .options.snow = opts), {
            corr1 <- matrix(rep(0, 2 * nc), nrow = 2, ncol = nc)
            for (j in 1:nc) {
                if (j > i) corr1[, j] <- r_p(totu[, i], totu[, j])
            }
            corr1
        })
        simplify2array(corr) -> corr
        rr <- corr[1, , ]
        rr <- rr + t(rr)
        diag(rr) <- 1

        pp <- corr[2, , ]
        pp <- pp + t(pp)
        rownames(rr) <- rownames(pp) <- colnames(totu)
        colnames(rr) <- colnames(pp) <- colnames(totu)
    } else {
        if (method == "spearman") {
            totu2 <- apply(totu2, 2, rank)
        } else if (method == "pearson") {
            totu2 <- totu2
        } else {
            stop('method must be one of "spearman","pearson"')
        }
        corr <- foreach::`%dopar%`(foreach::foreach(i = 1:nc, .options.snow = opts), {
            corr1 <- matrix(rep(0, 2 * ncol(totu2)), nrow = 2, ncol = ncol(totu2))
            for (j in 1:ncol(totu2)) {
                corr1[, j] <- r_p(totu[, i], totu2[, j])
            }
            corr1
        })
        simplify2array(corr) -> corr
        rr <- t(corr[1, , ])
        pp <- t(corr[2, , ])

        rownames(rr) <- rownames(pp) <- colnames(totu)
        colnames(rr) <- colnames(pp) <- colnames(totu2)
    }
    snow::stopCluster(cl)
    gc()
    pcutils::del_ps("doSNOW", "snow", "foreach")

    return(list(r = rr, p.value = pp))
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

#' Calculate distance for otutab
#'
#' @param totu t(otutab)
#' @param norm hellinger normalization in tax (default:TRUE)
#' @param method Dissimilarity index, partial match to "bray", "euclidean"...
#'
#' @return 1-dist, similarity
#' @export
#' @seealso \code{\link[vegan]{vegdist}}
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab) -> totu
#' par_sim(totu) -> sim_corr
par_sim <- function(totu, method = "bray", norm = FALSE) {
    otu <- t(totu)
    lib_ps("vegan", library = FALSE)

    if (norm) otu <- vegan::decostand(otu, "hellinger", 1) %>% as.data.frame()

    vegan::vegdist(otu, method = method) %>% as.matrix() -> dist

    sim <- 1 - dist
    p <- dist
    p[p != 0] <- 0

    return(list(r = sim, p.value = p))
}

# flag determined by the correlation table from one table or two
t_flag <- \(corr){
    if (!nrow(corr) == ncol(corr)) {
        return(FALSE)
    }
    if (!all(rownames(corr) == colnames(corr))) {
        return(FALSE)
    }
    return(TRUE)
}

#' p.adjust apply on a correlation table (matrix or data.frame)
#'
#' @param pp table
#' @param method see \code{\link[stats]{p.adjust}}
#' @param mode "all" for all values; "rows" adjust each row one by one; "columns" adjust each column one by one.
#'
#' @return matrix
#' @export
#'
#' @examples
#' matrix(abs(rnorm(100, 0.01, 0.1)), 10, 10) -> pp
#' p.adjust.table(pp)
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
