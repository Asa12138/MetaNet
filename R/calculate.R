#=======1.calculate========
#' Calculate spearman correlation for one or two t(otutab)
#'
#' @param totu t(otutab)
#' @param totu2 t(otutab) or NULL
#' @param threads parallel mode threads
#' @param filename the prefix of saved files
#' @param p.adjust.method
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
c_net_cal <- function(totu, totu2 = NULL, filename = "occor",threads=4,p.adjust.method = NULL) {
  #corr<-par_cor(totu,totu2,threads = threads,method = "spearman")
  corr<-f_cor(totu,totu2,method = "spearman")
  if(is.null(totu2))diag(corr$r)=0
  r <- round(corr$r, 5)
  p.value<-corr$p.value
  if(!is.null(p.adjust.method))p.value<-p.adjust.table(corr$p.value,p.adjust.method)
  p.value <- round(p.value, 5)
  # save the correlation result
  if(is.character(filename)){
    write.csv(r, paste0(filename, "_r.csv"))
    write.csv(p.value, paste0(filename, "_p.csv"))
  }
  return(list(r = r, p.value = p.value))
}

#' Fast correlation
#' @param totu t(otutab)
#' @param totu2 t(otutab) or NULL
#' @param method spearman or pearson
#' @export
#' @return a list with 2 elements:
#' \item{r}{spearman correlation}
#' \item{p.value}
#' @examples
#' data("otutab")
#' t(otutab[1:100, ]) -> totu
#' f_cor(totu,method="spearman") -> corr
f_cor<-function (totu,totu2=NULL, method = c("pearson", "spearman")){
  totu=as.matrix(totu)
  if(!is.null(totu2))totu2=as.matrix(totu2)
  method <- match.arg(method)
  r <- stats::cor(totu,totu2, method = method)
  df <- dim(totu)[1] - 2
  t <- r * sqrt(df/(1 - r^2))
  p <- -2 * expm1(pt(abs(t), df, log.p = TRUE))
  return(list(r=r,p.value=p))
}

#' Parallel calculate spearman correlation for one t(otutab)
#'
#' @rdname f_cor
#' @export
par_cor <- function(totu, totu2 = NULL,threads=2,method = c("spearman","pearson")[1],p.value=T){
  lib_ps("foreach")
  lib_ps("doSNOW")
  if(method=="spearman")totu <- apply(totu,2,rank)
  else if(method=="pearson")totu=totu
  else stop('method must be one of "spearman","pearson"')
  r_p_func <- \(pv)\(rx,ry){
    lxy <- sum((rx-mean(rx))*(ry-mean(ry)))
    lxx <- sum((rx-mean(rx))^2)
    lyy <- sum((ry-mean(ry))^2)
    r <- lxy/sqrt(lxx*lyy)
    p=0
    if(pv){
      n <- length(rx)
      t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
      p <- -2 * expm1(pt(abs(t), (n - 2), log.p = TRUE))
    }
    return(c(r,p))
  }
  r_p=r_p_func(p.value)

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

  return(list(r = rr,p.value = pp))
}

Hmisc_cor<-function(totu, totu2 = NULL,method = c("spearman","pearson")[1]){
  lib_ps("Hmisc")
  totu=as.matrix(totu)
  if(!is.null(totu2)){
    totu2=as.matrix(totu2)
    tmp=Hmisc::rcorr(totu,totu2,type = method)
    r=tmp$r[1:ncol(totu),(ncol(totu)+1):ncol(tmp$r)]
    p=tmp$P[1:ncol(totu),(ncol(totu)+1):ncol(tmp$r)]
    return(list(r = r,p.value = p))
  }
  tmp=Hmisc::rcorr(totu,type = method)
  p=tmp$P;p[is.na(p)]=0
  return(list(r = tmp$r,p.value = p))
}
#比较cor
if(F){
  a=matrix(rnorm(73*5160),nrow = 73)
  b=matrix(rnorm(73*6450),nrow = 73)
  system.time(f_cor(a,b)->c)

  method="spearman"
  t(otutab[1:100, ]) -> totu
  bench::mark(
    B=Hmisc::rcorr(totu,type =method )[["r"]],
    C=MetaNet::f_cor(totu,method = method)[["r"]],
    D=picante::cor.table(totu,cor.method = method)[["r"]],

    A=psych::corr.test(totu,method = method)[["r"]]
  )
  bench::mark(
    B={tmp=Hmisc::rcorr(totu,type =method )[["P"]];tmp[is.na(tmp)]=0;tmp},
    C=MetaNet::f_cor(totu,method = method)[["p.value"]],
    D=picante::cor.table(totu,cor.method = method)[["P"]],
    A={tmp=psych::corr.test(totu,method = method)[["p"]];tmp[upper.tri(tmp)]=0;tmp=tmp+t(tmp);tmp}
  )

#两个表
  totu2=t(otutab[101:150, ])
  bench::mark(
    B=Hmisc_cor(totu,totu2,method =method )[["r"]],
    C=MetaNet::f_cor(totu,totu2,method = method)[["r"]],
    #D=picante::cor.table(totu,totu2,cor.method = method)[["r"]],
    A=psych::corr.test(totu,totu2,method = method)[["r"]]
  )
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
  totu=as.data.frame(totu)
  set.seed(123)
  totu.sparcc <- sparcc(totu,iter = 10,inner_iter = 5)
  sparcc0 <- totu.sparcc$Cor  #sparcc correlation
  rownames(sparcc0)=colnames(sparcc0)=colnames(totu)

  #100 times bootstrap to get random matrix
  set.seed(123)
  reps = bootstrap
  #tmpd=tempdir()

  #parallel
  lib_ps("foreach","doSNOW")
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  p_boot=foreach (rep = 1:reps,
           .options.snow = opts,
           .packages = c("SpiecEasi")) %dopar%{
             totu.boot <- sample(totu, replace = TRUE)  #bootstrap
             totu.sparcc_boot <- sparcc(totu.boot,iter = 10,inner_iter = 5)  #sparcc same to above calculation
             sparcc_boot <- totu.sparcc_boot$Cor
             #colnames(sparcc_boot)=rownames(sparcc_boot) =colnames(totu.boot)
             sparcc_boot
             #write.table(sparcc_boot, paste(tmpd,'/sparcc_boot', rep, '.txt', sep = ''), sep = '\t', col.names = NA, quote = FALSE)  #save
           }
  stopCluster(cl)
  #get the pseudo-pvalue by bootstrap
  #p_boot1=simplify2array(p_boot)
  p = sparcc0
  p[p!=0]=0
  lapply(p_boot,\(x)(x>sparcc0)*1)->s
  p=apply(simplify2array(s), 1:2, sum)
  p <- p / reps
  colnames(p)=rownames(p) =colnames(sparcc0)
  # save the correlation result
  if(is.character(filename)){
    write.csv(sparcc0, paste0(filename, "_r.csv"))
    write.csv(p, paste0(filename, "_p.csv"))
  }
  #unlink(tmpd,recursive = T,force = T)
  return(list(r = sparcc0, p.value = p))
}

#' Calulate distance for otutab
#'
#' @param otutab an otutab data.frame, samples are columns, taxs are rows.
#' @param method Dissimilarity index, partial match to "bray", "euclidean"...
#' @param spe_nwk a phylo tree
#'
#' @return dist
#' @export
#' @seealso \code{\link[vegan]{vegdist};\link[picante]{unifrac}}
#' @examples
#' par_sim(totu)->sim_corr
par_sim<-function(totu,method="bray"){
  otu=t(totu)

  otu=vegan::decostand(otu,"hellinger",1)%>%as.data.frame()

  lib_ps("vegan");vegan::vegdist(otu,method = method)%>%as.matrix()->a
  sim=1-a
  p=a
  p[p!=0]=0
  return(list(r = sim, p.value = p))
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

