# =======1.calculate========
#' Calculate spearman correlation for one or two t(otutab)
#'
#' @param totu t(otutab)
#' @param totu2 t(otutab2) or NULL
#' @param method spearman, pearson
#' @param filename the prefix of saved .corr file or FALSE
#' @param p.adjust.method see \code{\link[stats]{p.adjust}}
#' @param threads threads, default: 1.
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
                      p.adjust.method = NULL,threads=1) {
  if (!is.null(totu2)) {
    tls <- check_tabs(totu, totu2)
    totu <- tls[[1]]
    totu2 <- tls[[2]]
  } else {
    totu <- check_tabs(totu)
  }

  if(method%in%c("spearman", "pearson")) corr <- f_cor(totu, totu2, method = method)
  else if (method=="sparcc"){
    if(!is.null(totu2))warning("sparcc only take the totu, ignore the totu2.")
    corr=par_sparcc(totu,threads = threads)
  }

  # if(threads>1){
  # corr<-par_cor(totu,totu2,threads = threads,method = "spearman")
  # }

  # if (is.null(totu2)) diag(corr$r) <- 0

  if (!is.null(p.adjust.method)) {
    p.adjust <- p.adjust.table(corr$p.value, p.adjust.method)
    res=list(r = corr$r, p.value = corr$p.value, p.adjust = p.adjust)
  } else {
    res=list(r = corr$r, p.value = corr$p.value)
  }

  class(res)="corr"

  # save the correlation result
  if(is.logical(filename)){if(filename)filename=paste0("c_net_",date())}
  if (is.character(filename)) {
    saveRDS(res,file = paste0(filename,".corr"))
  }

  # if (is.character(filename)) {
  #   # if(!dir.exists("net_res"))dir.create("net_res/")
  #   write.csv(round(corr$r, 5), paste0(filename, "_r.csv"))
  #   write.csv(round(corr$p.value, 5), paste0(filename, "_p.csv"))
  #   write.csv(round(p.adjust, 5), paste0(filename, "_p_adj.csv"))
  # }

  return(res)
}

#' Print method for 'corr' objects
#'
#' @param x 'corr' object
#' @param ... additional arguments
#'
#' @exportS3Method
#' @method print corr
print.corr <- function(x, ...) {
  cat("Correlation table:")
  cat("Table dimensions:", nrow(x$r), "rows,", ncol(x$r), "columns\n")
}


#' Check tables
#' @param ... tables
#' @return format tables
#' @export
check_tabs=function(...){
  tables=list(...)
  if(all(class(tables[[1]])=="list"))tables=tables[[1]]

  names(tables)=NULL
  comm=Reduce(intersect,lapply(tables, rownames))
  if(length(comm)<2)stop("There are ",length(comm)," common sample! Can not calculate correlation.")
  if(all(lapply(tables, \(i)identical(rownames(i),comm))%>%unlist()))message("All samples matched.")
  else message("Extract, ",length(comm) ," commmon samples.")

  dup=lapply(tables, colnames)%>%do.call(c,.)
  dup=dup[duplicated(dup)]
  if(length(dup)>0)stop("Duplicated colnames found: ",paste0(dup,collapse = ", "), "\nPlease check colnames of input tables.")
  else message("All objects are OK.")

  tables=lapply(tables, \(i)i[comm,])
  if(length(tables)==1)tables=tables[[1]]
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
#' @examples
#' \dontrun{
#' \donttest{
#'  input_corr("test")->corr
#' }}
input_corr <- function(filename) {
  # r <- read.csv(paste0(filename, "_r.csv"), row.names = 1, check.names = FALSE)
  # p.value <- read.csv(paste0(filename, "_p.csv"), row.names = 1, check.names = FALSE)
  # p.adjust <- read.csv(paste0(filename, "_p_adj.csv"), row.names = 1, check.names = FALSE)
  if(!grepl(".corr",filename))filename=paste0(filename,".corr")
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

  pcutils::lib_ps("foreach","doSNOW","snow")
  pb <- utils::txtProgressBar(max = nc, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- snow::makeCluster(threads)
  doSNOW::registerDoSNOW(cl)
  if (is.null(totu2)) {
    corr <- foreach::foreach(i = 1:nc, .options.snow = opts) %dopar% {
      corr1 <- matrix(rep(0, 2 * nc), nrow = 2, ncol = nc)
      for (j in 1:nc) {
        if (j > i) corr1[, j] <- r_p(totu[, i], totu[, j])
      }
      corr1
    }
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
    corr <- foreach::foreach(i = 1:nc, .options.snow = opts) %dopar% {
      corr1 <- matrix(rep(0, 2 * ncol(totu2)), nrow = 2, ncol = ncol(totu2))
      for (j in 1:ncol(totu2)) {
        corr1[, j] <- r_p(totu[, i], totu2[, j])
      }
      corr1
    }
    simplify2array(corr) -> corr
    rr <- t(corr[1, , ])
    pp <- t(corr[2, , ])

    rownames(rr) <- rownames(pp) <- colnames(totu)
    colnames(rr) <- colnames(pp) <- colnames(totu2)
  }
  snow::stopCluster(cl)
  gc()
  pcutils::del_ps("doSNOW","snow","foreach")

  return(list(r = rr, p.value = pp))
}

Hmisc_cor <- function(totu, totu2 = NULL, method = c("spearman", "pearson")[1]) {
  lib_ps("Hmisc",library = FALSE)
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

# comparison
if (FALSE) {
  # a <- matrix(rnorm(73 * 5160), nrow = 73)
  # b <- matrix(rnorm(73 * 6450), nrow = 73)
  # system.time(f_cor(a, b) -> c)
  #
  # method <- "spearman"
  # t(otutab[1:100, ]) -> totu
  # bench::mark(
  #   B = Hmisc::rcorr(totu, type = method)[["r"]],
  #   C = MetaNet::f_cor(totu, method = method)[["r"]],
  #   D = picante::cor.table(totu, cor.method = method)[["r"]],
  #   A = psych::corr.test(totu, method = method)[["r"]]
  # )
  # bench::mark(
  #   B = {
  #     tmp <- Hmisc::rcorr(totu, type = method)[["P"]]
  #     tmp[is.na(tmp)] <- 0
  #     tmp
  #   },
  #   C = MetaNet::f_cor(totu, method = method)[["p.value"]],
  #   D = picante::cor.table(totu, cor.method = method)[["P"]],
  #   A = {
  #     tmp <- psych::corr.test(totu, method = method)[["p"]]
  #     tmp[upper.tri(tmp)] <- 0
  #     tmp <- tmp + t(tmp)
  #     tmp
  #   }
  # )
  #
  # # two tables
  # totu2 <- t(otutab[101:150, ])
  # bench::mark(
  #   B = Hmisc_cor(totu, totu2, method = method)[["r"]],
  #   C = MetaNet::f_cor(totu, totu2, method = method)[["r"]],
  #   # D=picante::cor.table(totu,totu2,cor.method = method)[["r"]],
  #   A = psych::corr.test(totu, totu2, method = method)[["r"]]
  # )
}

#' SparCC correlation for a otutab
#'
#' @param totu an t(otutab)
#' @param threads default: 1
#' @param bootstrap default:100, recommend not less than 100
#'
#' @return list contain a correlation matrix and a bootstrap p_value matrix
#' @export
#'
#' @examples
#' \dontrun{
#' \donttest{
#' data("otutab", package = "pcutils")
#' t(otutab) -> totu
#' par_sparcc(totu[, 1:50],4) -> sparcc_corr
#' }}
par_sparcc <- function(totu, threads = 1, bootstrap = 100) {
  # sparcc
  lib_ps("SpiecEasi",library = FALSE)
  totu <- as.data.frame(totu)
  set.seed(123)
  totu.sparcc <- SpiecEasi::sparcc(totu, iter = 10, inner_iter = 5)
  sparcc0 <- totu.sparcc$Cor # sparcc correlation
  rownames(sparcc0) <- colnames(sparcc0) <- colnames(totu)

  #parallel
  reps <- bootstrap
  #main function
  loop=function (i)
  {
    totu.boot <- sample(totu, replace = TRUE)
    totu.sparcc_boot <- SpiecEasi::sparcc(totu.boot, iter = 10,
                                          inner_iter = 5)
    sparcc_boot <- totu.sparcc_boot$Cor
    sparcc_boot
  }
  {
    if(threads>1){
      pcutils::lib_ps("foreach","doSNOW","snow")
      pb <- utils::txtProgressBar(max =reps, style = 3)
      opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
      cl <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cl)
      res <- foreach::foreach(i = 1:reps,.options.snow = opts,
                              .packages = c()) %dopar% {
                                loop(i)
                              }
      snow::stopCluster(cl)
      gc()
      pcutils::del_ps("doSNOW","snow","foreach")
    }
    else {
      res <-lapply(1:reps, loop)
    }}
  #simplify method
  pcutils::del_ps("foreach","doSNOW")

  # get the pseudo-pvalue by bootstrap
  p <- sparcc0
  p[p != 0] <- 0
  lapply(res, \(x)(x > sparcc0) * 1) -> s
  p <- apply(simplify2array(s), 1:2, sum)
  p <- p / reps
  colnames(p) <- rownames(p) <- colnames(sparcc0)
  return(list(r = sparcc0, p.value = p))
}


#' Calculate distance for otutab
#'
#' @param totu t(otutab)
#' @param norm hellinger normalization in tax (default:TRUE)
#' @param method Dissimilarity index, partial match to "bray", "euclidean"...
#'
#' @return dist
#' @export
#' @seealso \code{\link[vegan]{vegdist}}
#' @examples
#' data("otutab", package = "pcutils")
#' t(otutab) -> totu
#' par_sim(totu) -> sim_corr
par_sim <- function(totu, method = "bray", norm = TRUE) {
  otu <- t(totu)
  lib_ps("vegan",library = FALSE)

  if (norm) otu <- vegan::decostand(otu, "hellinger", 1) %>% as.data.frame()

  vegan::vegdist(otu, method = method) %>% as.matrix() -> a

  sim <- 1 - a
  p <- a
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
#'
#' @return matrix
#' @export
#'
#' @examples
#' matrix(abs(rnorm(100, 0.01)), 10, 10) -> pp
#' p.adjust.table(pp)
p.adjust.table <- \(pp, method = "BH"){
  pp <- as.matrix(pp)
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
  return(pp)
}
