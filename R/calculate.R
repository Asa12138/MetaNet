#=======1.calculate========
#' Calculate spearman correlation for one or two t(otutab)
#'
#' @param totu t(otutab)
#' @param totu2 t(otutab) or NULL
#' @param parallel open parallel mode?(default:F)
#' @param threads parallel mode threads
#' @param filename the prefix of saved files
#'
#' @return a list with 2 elements:
#' \item{r}{spearman correlation}
#' \item{p.value}{p.adjust.method = 'NULL}
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab) -> totu
#' c_net_cal(totu) -> corr
#' metadata[,3:10] -> env
#' c_net_cal(totu, env) -> corr2
c_net_cal <- function(totu, totu2 = NULL, filename = "occor",threads=4) {
  corr<-par_cor(totu,totu2,threads = threads,method = "spearman",p.adjust.methods=NULL)
  r <- round(corr$r, 5)
  p.value <- round(corr$p.value, 5)
  if(is.null(totu2))diag(r)=0
  # save the correlation result
  if(is.character(filename)){
    write.csv(r, paste0(filename, "_r.csv"))
    write.csv(p.value, paste0(filename, "_p.csv"))
  }
  return(list(r = r, p.value = p.value))
}

#' Parallel calculate spearman correlation for one t(otutab)
#'
#' @param totu t(otutab)
#' @param totu2 t(otutab) or NULL
#' @param threads parallel mode threads, default 2
#' @param method spearman or not
#' @param p.adjust.methods how to adjust p-value (default:NULL, e.g.fdr)
#'
#' @return a list with 2 elements:
#' \item{r}{spearman correlation}
#' \item{p.value}{p.adjust.method = 'NULL}
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab[1:100, ]) -> totu
#' par_cor(totu,4) -> corr
#'
par_cor <- function(totu, totu2 = NULL,threads=2,method = "spearman",p.adjust.methods=NULL){
  lib_ps("foreach")
  lib_ps("doSNOW")
  if(method=="spearman")totu <- apply(totu,2,rank)
  r_p <- \(rx,ry){
    n <- length(rx)
    lxy <- sum((rx-mean(rx))*(ry-mean(ry)))
    lxx <- sum((rx-mean(rx))^2)
    lyy <- sum((ry-mean(ry))^2)
    r <- lxy/sqrt(lxx*lyy)
    t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
    p <- -2 * expm1(pt(abs(t), (n - 2), log.p = TRUE))
    return(c(r,p))
  }

  nc <- ncol(totu)
  pb <- txtProgressBar(max =nc, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  if(is.null(totu2)){
    corr <- foreach (i = 1:nc,.options.snow=opts) %dopar%{
      corr1 <- matrix(rep(0,2*nc),nrow = 2,ncol=nc)
      for(j in 1:nc) {
        if(j > i) corr1[,j] <- r_p(totu[,i],totu[,j])
      }
      corr1
    }
    simplify2array(corr)->corr
    rr <- corr[1,,]
    rr <- rr+t(rr)
    diag(rr) <- 1

    pp <- corr[2,,]
    pp <- pp+t(pp)
    rownames(rr) <-rownames(pp) <- colnames(totu)
    colnames(rr) <-colnames(pp) <- colnames(totu)
  }
  else{
    corr <- foreach (i = 1:nc,.options.snow=opts) %dopar%{
      corr1 <- matrix(rep(0,2*nc),nrow = 2,ncol=ncol(totu2))
      for(j in 1:ncol(totu2)) {
        corr1[,j] <- r_p(totu[,i],totu2[,j])
      }
      corr1
    }
    simplify2array(corr)->corr
    rr <- t(corr[1,,])
    pp <- t(corr[2,,])

    rownames(rr) <-rownames(pp) <- colnames(totu)
    colnames(rr) <-colnames(pp) <- colnames(totu2)
  }
  stopCluster(cl)

  if(!is.null(p.adjust.methods))pp<-p.adjust.table(pp,p.adjust.methods)
  return(list(r = rr,p.value = pp))
}


#' SparCC correlation for a otutab
#'
#' @param totu an t(otutab)
#' @param threads default: 4
#' @param filename saved file name
#' @param bootstrap default:100, recommend not less than 100
#'
#' @return list contain a correlation matrix and a bootstrap p_value matrix
#' @export
#'
#' @examples
#' par_sparcc(totu)->sparcc_corr
par_sparcc<-function(totu,filename="sparcc",threads=4,bootstrap=100){
  #sparcc
  lib_ps("SpiecEasi")

  set.seed(123)
  totu.sparcc <- sparcc(totu,iter = 10,inner_iter = 5)
  sparcc0 <- totu.sparcc$Cor  #sparcc correlation
  rownames(sparcc0)=colnames(sparcc0)=colnames(totu)

  #100 times bootstrap to get random matrix
  set.seed(123)
  reps = bootstrap
  tmpd=tempdir()

  #parallel
  lib_ps("foreach","doSNOW")
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  foreach (rep = 1:reps,
           .options.snow = opts,
           .packages = c("SpiecEasi")) %dopar%{
             totu.boot <- sample(totu, replace = TRUE)  #bootstrap
             totu.sparcc_boot <- sparcc(totu.boot,iter = 10,inner_iter = 5)  #sparcc same to above calculation
             sparcc_boot <- totu.sparcc_boot$Cor
             colnames(sparcc_boot)=rownames(sparcc_boot) =colnames(totu.boot)
             write.table(sparcc_boot, paste(tmpd,'/sparcc_boot', rep, '.txt', sep = ''), sep = '\t', col.names = NA, quote = FALSE)  #save
           }
  stopCluster(cl)
  #get the pseudo-pvalue by bootstrap
  p <- sparcc0
  p[p!=0] <- 0
  for (i in 1:reps) {
    p_boot <- read.delim(paste(tmpd,'/sparcc_boot', i, '.txt', sep = ''), sep = '\t', row.names = 1)
    p[abs(p_boot)>=abs(sparcc0)] <- p[abs(p_boot)>=abs(sparcc0)] + 1
  }
  p <- p / reps

  # save the correlation result
  if(is.character(filename)){
    write.csv(sparcc0, paste0(filename, "_r.csv"))
    write.csv(p, paste0(filename, "_p.csv"))
  }
  #unlink(tmpd,recursive = T,force = T)
  return(list(r = sparcc0, p.value = p))
}


#' p.adjust apply on a table (matrix or data.frame)
#'
#' @param pp table
#' @param method see \code{\link[stats]{p.adjust}}
#'
#' @return matrix
#' @export
#'
#' @examples
#' matrix(abs(rnorm(100,0.01)),10,10)->pp
#' p.adjust.table(pp)
p.adjust.table<-\(pp,method="BH"){
  pp<-as.matrix(pp)
  #self-cor or two-tables cor?
  flag <- \(corr){
    if (!nrow(corr) == ncol(corr)) {
      return(FALSE)
    }
    if (!all(rownames(corr) == colnames(corr))) {
      return(FALSE)
    }
    return(TRUE)
  }
  if(flag(pp)){
    lp <- lower.tri(pp)
    pa <- pp[lp]
    pa <- p.adjust(pa,method)
    pp[lower.tri(pp, diag = FALSE)] <- pa
    pp[upper.tri(pp,diag = F)]<-0
    pp<-pp+t(pp)
  }
  else{
    pp1<-p.adjust(pp,method)
    pp<-matrix(pp1,nrow(pp),ncol(pp),dimnames = list(rownames(pp),colnames(pp)))
  }
  return(pp)
}

#' Import corr from .csv file
#'
#' @param filename
#'
#' @return a list with 7 elements:
#' \item{r}{spearman correlation}
#' \item{p.value}{p.adjust.method = 'BH'}
#' @export
#'
#' @examples
#' \dontrun{
#' # input_corr("occor")->corr
#' }
input_corr <- function(filename) {
  r <- read.csv(paste0(filename, "_r.csv"), row.names = 1,check.names = F)
  p.value <- read.csv(paste0(filename, "_p.csv"), row.names = 1,check.names = F)
  return(list(r = r, p.value = p.value))
}

mmscale<-function(x,min_s=0,max_s=1){
  if((max(x)-min(x))==0)return(rep((min_s+max_s)/2,length(x)))
  min_s+(x-min(x))/(max(x)-min(x))*(max_s-min_s)
}
