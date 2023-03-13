#======0.preprocess========

#' Transfer your data
#' @param df dataframe
#' @param method "cpm","minmax","acpm","total","log", "max", "frequency", "normalize", "range", "rank", "rrank",
#' "standardize", "pa", "chi.square", "hellinger", "log", "clr", "rclr", "alr"
#' @param margin 1 for row and 2 for column(default: 2)
#' @param ... additional
#' @export
#' @examples
#'data(otutab)
#'trans(otutab,method="cpm")
#' @seealso \code{\link[vegan]{decostand}}
trans<-function(df,method = "normalize",margin=2,...){
  all=c("cpm","minmax","acpm","total", "log1","max", "frequency", "normalize", "range", "rank", "rrank",
        "standardize", "pa", "chi.square", "hellinger", "log", "clr", "rclr", "alr")
  if (!method%in%all)stop("methods should be one of ",all)
  if(method=="cpm"){
    df=apply(df, margin, \(x){x*10**6/sum(x)})
  }
  else if(method=="minmax"){
    df=apply(df,margin,mmscale,...)
  }
  else if(method=="acpm"){
    df=asinh(apply(df, margin, \(x){x*10**6/sum(x)}))
  }
  else if(method=="log1"){
    df=log(df+1,...)
  }
  else df=vegan::decostand(df,method = method,margin,...)
  return(data.frame(df,check.names = F))
}

#' Filter your data
#'
#' @param tab dataframe
#' @param sum the rowsum should bigger than sum(default:10)
#' @param exist the exist number bigger than exist(default:1)
#'
#' @return input object
#' @export
#'
#' @examples
#'data(otutab)
#'guolv(otutab)
guolv<-function(tab,sum=0,exist=1){
  tab[rowSums(tab)>sum,]->tab
  tab[rowSums(tab>0)>exist,]->tab
  return(tab)
}

#' Group your data
#'
#' @param otutab dataframe
#' @param group group vector
#' @param margin 1 for row and 2 for column(default: 2)
#' @param act do (default: mean)
#' @rdname hebing
#' @return input object
#' @export
#'
#' @examples
#' data(otutab)
#' hebing(otutab,metadata$Group)
hebing<-function(otutab,group,margin=2,act='mean'){
  if (margin==2) {
    aggregate(t(otutab),FUN=act,by=list(factor(group)))->a
    a[,-1]->a
    data.frame(t(a))->a
    levels(factor(group))->colnames(a)
  }
  else{
    aggregate(otutab,FUN=act,by=list(factor(group)))->a
    a[,-1]->a
    levels(factor(group))->rownames(a)
  }
  return(a)
}


#=======1.calculate========
#' Calculate spearman correlation for one or two t(otutab)
#'
#' @param totu t(otutab)
#' @param totu2 t(otutab) or NULL
#' @param method spearman, pearson
#' @param filename the prefix of saved files
#' @param p.adjust.method see \code{\link[stats]{p.adjust}}
#'
#' @return a list with 2 elements:
#' \item{r}{default: spearman correlation}
#' \item{p.value}{}
#' \item{p.adjust}{default p.adjust.method = NULL}
#' @export
#'
#' @examples
#' data("otutab")
#' t(otutab) -> totu
#' c_net_cal(totu) -> corr
#' metadata[,3:10] -> env
#' c_net_cal(totu, env) -> corr2
c_net_cal <- function(totu, totu2 = NULL,method = "spearman", filename = "occor",p.adjust.method = NULL) {
  corr<-f_cor(totu,totu2,method = method)
  # if(threads>1){
  # corr<-par_cor(totu,totu2,threads = threads,method = "spearman")
  # }
  if(is.null(totu2))diag(corr$r)=0

  if(!is.null(p.adjust.method))p.adjust<-p.adjust.table(corr$p.value,p.adjust.method)
  else p.adjust=corr$p.value

  # save the correlation result
  if(is.character(filename)){
    #if(!dir.exists("net_res"))dir.create("net_res/")
    write.csv(round(corr$r,5), paste0(filename, "_r.csv"))
    write.csv(round(corr$p.value, 5), paste0(filename, "_p.csv"))
    write.csv(round(p.adjust, 5), paste0(filename, "_p_adj.csv"))
  }
  return(list(r = corr$r, p.value = corr$p.value,p.adjust=p.adjust))
}

#' Import corr from .csv file
#'
#' @param filename
#'
#' @return a list with 7 elements:
#' \item{r}{default: spearman correlation}
#' \item{p.value}{}
#' \item{p.adjust}{default p.adjust.method = NULL}
#' @export
#'
#' @examples
#' \dontrun{
#' # input_corr("occor")->corr
#' }
input_corr <- function(filename) {
  r <- read.csv(paste0(filename, "_r.csv"), row.names = 1,check.names = F)
  p.value <- read.csv(paste0(filename, "_p.csv"), row.names = 1,check.names = F)
  p.adjust <- read.csv(paste0(filename, "_p_adj.csv"), row.names = 1,check.names = F)
  return(list(r = r, p.value = p.value,p.adjust=p.adjust))
}

#' Fast correlation
#' @param totu t(otutab)
#' @param totu2 t(otutab) or NULL
#' @param method spearman or pearson
#' @export
#'
#' @return a list with 2 elements:
#' \item{r}{}
#' \item{p.value}{}
#'
#' @examples
#' data("otutab")
#' t(otutab[1:100, ]) -> totu
#' f_cor(totu,method="spearman") -> corr
f_cor<-function (totu,totu2=NULL, method = c("pearson", "spearman")){
  method = match.arg(method,c("pearson", "spearman"))
  totu=as.matrix(totu)
  if(!is.null(totu2)){
    totu2=as.matrix(totu2)
    r <- stats::cor(totu,totu2, method = method)
    df <- dim(totu)[1] - 2
    t <- r * sqrt(df/(1 - r^2))
    p <- -2 * expm1(pt(abs(t), df, log.p = TRUE))}
  else{
    r <- stats::cor(totu, method = method)
    p=r;p[p!=0]=0
    r_tri=r[upper.tri(r)]
    df <- dim(totu)[1] - 2
    t <- r_tri * sqrt(df/(1 - r_tri^2))
    p[upper.tri(p)] <- -2 * expm1(pt(abs(t), df, log.p = TRUE))
    p=p+t(p)
  }
  return(list(r=r,p.value=p))
}

#' Parallel calculate correlation
#'
#' @rdname f_cor
#' @export
par_cor <- function(totu, totu2 = NULL,threads=2,method = c("spearman","pearson")){
  method = match.arg(method,c("pearson", "spearman"))
  lib_ps("foreach")
  lib_ps("doSNOW")
  if(method=="spearman")totu <- apply(totu,2,rank)
  else if(method=="pearson")totu=totu
  else stop('method must be one of "spearman","pearson"')
  r_p <- \(rx,ry){
    lxy <- sum((rx-mean(rx))*(ry-mean(ry)))
    lxx <- sum((rx-mean(rx))^2)
    lyy <- sum((ry-mean(ry))^2)
    r <- lxy/sqrt(lxx*lyy)
    p=0
    n <- length(rx)
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
    if(method=="spearman")totu2 <- apply(totu2,2,rank)
    else if(method=="pearson")totu2=totu2
    else stop('method must be one of "spearman","pearson"')
    corr <- foreach (i = 1:nc,.options.snow=opts) %dopar%{
      corr1 <- matrix(rep(0,2*ncol(totu2)),nrow = 2,ncol=ncol(totu2))
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

#comparison
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

#two tables
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
#' @param bootstrap default:100, recommend not less than 100
#'
#' @return list contain a correlation matrix and a bootstrap p_value matrix
#' @export
#'
#' @examples
#' par_sparcc(totu[,1:50])->sparcc_corr
par_sparcc<-function(totu,threads=4,bootstrap=100){
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

  loop<-\(i){
    totu.boot <- sample(totu, replace = TRUE)  #bootstrap
    totu.sparcc_boot <- sparcc(totu.boot,iter = 10,inner_iter = 5)  #sparcc same to above calculation
    sparcc_boot <- totu.sparcc_boot$Cor
    sparcc_boot
    #write.table(sparcc_boot, paste(tmpd,'/sparcc_boot', rep, '.txt', sep = ''), sep = '\t', col.names = NA, quote = FALSE)  #save
  }

  if(threads==1)lapply(1:reps,FUN = loop)->p_boot
  else{
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
               loop(rep)
             }
    stopCluster(cl)
    }

  #get the pseudo-pvalue by bootstrap
  #
  #p_boot1=simplify2array(p_boot)
  p = sparcc0
  p[p!=0]=0
  lapply(p_boot,\(x)(x>sparcc0)*1)->s
  p=apply(simplify2array(s), 1:2, sum)
  p <- p / reps
  colnames(p)=rownames(p) =colnames(sparcc0)
  # save the correlation result
  # if(is.character(filename)){
  #   #if(!dir.exists("net_res"))dir.create("net_res/")
  #   write.csv(sparcc0, paste0(filename, "_r.csv"))
  #   write.csv(p, paste0(filename, "_p.csv"))
  # }
  #unlink(tmpd,recursive = T,force = T)
  return(list(r = sparcc0, p.value = p))
}

#' Calculate distance for otutab
#'
#' @param totu an t(otutab)
#' @param norm hellinger normalization in tax (default:T)
#' @param method Dissimilarity index, partial match to "bray", "euclidean"...
#'
#' @return dist
#' @export
#' @seealso \code{\link[vegan]{vegdist}}
#' @examples
#' par_sim(totu)->sim_corr
par_sim<-function(totu,method="bray",norm=T){
  otu=t(totu)
  lib_ps("vegan")

  if(norm)otu=vegan::decostand(otu,"hellinger",1)%>%as.data.frame()
  vegan::vegdist(otu,method = method)%>%as.matrix()->a
  sim=1-a
  p=a
  p[p!=0]=0

  return(list(r = sim, p.value = p))
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



