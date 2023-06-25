#========7.stability===========
#' Robust_test for a network
#'
#' @param go_ls a igraph object or igraph list.
#' @param partial how much percent vertexes be removed in total
#' @param step how many nodes be removed each time?
#' @param reps simulation number
#' @param threads threads
#'
#' @return dataframe(robustness class)
#' @export
#'
#' @examples
#' \donttest{
#' \dontrun{
#' data("c_net")
#' robust_test(co_net,step=4)->robust_res
#' plot(robust_res,index="ave_degree",mode=2)
#' }}
robust_test <- function(go_ls, partial = 0.5,step=10,reps=9,threads=1){
  if("igraph"%in%class(go_ls)){
    robustness_res=robust_test_in(go_ls,partial = partial,step=step,reps=reps,threads=threads)
    robustness_res=data.frame(robustness_res,group="Net")
  }
  else {
    if(!is.igraph(go_ls[[1]]))stop("No igraph or igraph-list.")
    robustness_res=lapply(names(go_ls), \(i){
      tmp=robust_test_in(go_ls[[i]],partial = partial,step=step,reps=reps,threads=threads)
      data.frame(tmp,group=i)
    })
    robustness_res=do.call(rbind,robustness_res)
    robustness_res$group=factor(robustness_res$group,levels = names(go_ls))
  }
  class(robustness_res)=c("robust","data.frame")
  return(robustness_res)
}

robust_test_in <- function(go, partial = 0.5,step=10,reps=9,threads=1) {
  lib_ps("igraph",library = FALSE)

  cal_del<-\(go,partial,step,rep){
    nodes <- length(igraph::V(go))
    floor(nodes * partial) -> del_i
    del_i_indexs <- data.frame()
    sequ=seq(0,del_i,step)
    if(sequ[length(sequ)]<del_i)sequ=c(sequ,del_i)
    res=lapply(sequ, \(i){
      # remove i nodes in the network
      remove_node <- sample(1:nodes, i)
      dp <- igraph::delete.vertices(go, remove_node)
      dp=igraph::delete.vertices(dp,igraph::V(dp)[igraph::degree(dp)==0])

      # calculate network parameters
      tmp_ind=(MetaNet::net_par(dp,mode = "n")$n_index)
      data.frame(tmp_ind, i = i)
    })
    del_i_indexs=do.call(rbind,res)

    # for (i in sequ) {
    #   # remove i nodes in the network
    #   remove_node <- sample(1:nodes, i)
    #   dp <- igraph::delete.vertices(go, remove_node)
    #   dp=igraph::delete.vertices(dp,igraph::V(dp)[igraph::degree(dp)==0])
    #
    #   # calculate network parameters
    #   tmp_ind=(MetaNet::net_par(dp,mode = "n")$n_index)
    #   del_i_indexs <- rbind(
    #     del_i_indexs,
    #     data.frame(tmp_ind, i = i)
    #   )
    # }

    return(data.frame(del_i_indexs,"rep"=rep))
  }

  #parallel
  #main function
  loop=function (i)
  {
    cal_del(go, partial, step, i)
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
  del_i_indexs=do.call(rbind,res)

  del_i_indexs$i_ratio=del_i_indexs$i/length(V(go))%>%round(.,3)
  class(del_i_indexs)<-c("robust",class(del_i_indexs))
  return(del_i_indexs)
}

#' Plot robust
#'
#' @param x `robust_test()` result (robust object)
#' @param mode plot mode, 1~3
#' @param indexes indexes selected to show
#' @param use_ratio use the delete nodes ratio rather than nodes number
#' @param ... additional arguments for \code{\link[pcutils]{group_box}}
#'
#' @return a ggplot
#' @method plot robust
#' @exportS3Method
#' @rdname robust_test
plot.robust<-function(x, indexes = c("nat_connectivity", "ave_path_len", "ave_degree"),use_ratio=FALSE,mode=1,...){
  robust_res=x
  lib_ps("reshape2",library = FALSE)
  xlab="Removed_nodes"
  if(use_ratio){
    robust_res$i=robust_res$i_ratio
    xlab="Removed_nodes_ratio"
  }
  robust_res %>%
    dplyr::select(i,group, c("ave_degree","nat_connectivity")) %>%
    reshape2::melt(id.var = c("i","group")) -> pdat

  pdat%>%dplyr::group_by(i,variable,group)%>%dplyr::summarise(mean=mean(value),sd=sd(value),se=sd/sqrt(length(value)))->sdd

  if(mode==1){
    p <- ggplot(sdd, aes(x=i, y=mean,col=group)) +
      geom_line()+
      geom_errorbar(data = sdd,aes(ymax=mean+se,ymin=mean-se))+
      #geom_smooth(se = FALSE,method = "loess",formula = 'y ~ x') +
      facet_wrap(. ~ variable, scales = "free")+labs(x=xlab,y=NULL)+
      theme_bw()+
      theme(legend.position = "none")
  }
  if(mode==2){
    lib_ps("ggpmisc",library = FALSE)
    p=ggplot(sdd,aes(x=i, y=mean,col=group))+
      geom_point(size=0.2,alpha=0.4)+
      geom_smooth(se = FALSE,method = "loess",formula = 'y ~ x')+
      ggpmisc::stat_poly_eq(
        aes(label = paste(after_stat(eq.label), after_stat(adj.rr.label), sep = '~~~~~')),
        formula = y ~ x,  parse = TRUE,label.x = "right",
        size = 3,
      )+
      facet_wrap(. ~ variable, scales = "free")+labs(x=xlab,y=NULL)+theme_bw()+
      theme(legend.position = "none")
  }
  if(mode==3){
    robust_res%>%dplyr::select(i,rep,group,nat_connectivity)->pdat
    #robust_res%>%filter(i==0)%>%select(num_nodes,group)%>%distinct()%>%left_join(pdat,.,by="group")->pdat
    #pdat%>%mutate(i=i/num_nodes)->pdat
    pdat%>%
      #filter(i<0.4)%>%
      dplyr::group_by(rep,group)%>%
      dplyr::summarise(slope=coefficients(lm(nat_connectivity~i))[2])->slope

    p=pcutils::group_box(slope["slope"],group="group",metadata = slope,alpha = TRUE,...)+theme_bw()
  }
  return(p)
}

#' Cohesion calculate
#'
#' @param otutab otutab
#' @param reps iteration time
#' @param threads threads
#' @param mycor a correlation matrix you want to use, skip the null model build when mycor is not NULL, default: NULL
#'
#' @return Cohesion object: a list with two dataframe
#' @export
#' @references 1. Herren, C. M. & McMahon, K. Cohesion: a method for quantifying the connectivity of microbial communities. (2017) doi:10.1038/ismej.2017.91.
#' @examples
#' \donttest{
#' \dontrun{
#' Cohesion(otutab[1:50,])->cohesion_res
#' plot.cohesion(cohesion_res,group = "Group",metadata = metadata,mode = 1)
#' plot.cohesion(cohesion_res,group = "Group",metadata = metadata,mode = 2)
#' }}
Cohesion<-function(otutab,reps= 200,threads=1,mycor=NULL){
  d <- t(otutab)
  rel.d <- d / rowSums(d)

  if(is.null(mycor)){
    # Create observed correlation matrix
    cor.mat.true <- cor(rel.d)
    nc <- ncol(rel.d)

    loop=\(which.taxon){
      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()

      for(i in 1:reps){
        #Create empty matrix of same dimension as rel.d
        perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)

        #For each otu
        for(j in 1:dim(rel.d)[2]){
          # Replace the original taxon vector with a permuted taxon vector
          perm.rel.d[, j ] <- sample(rel.d[ ,j ])
        }

        # Do not randomize focal column
        perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]

        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d,perm.rel.d[, which.taxon])

        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null)

      }

      # Save the median correlations between the focal taxon and all other taxa
      apply(perm.cor.vec.mat, 1, median)
      }
    {
      if(threads>1){
        pcutils::lib_ps("foreach","doSNOW","snow")
        pb <- utils::txtProgressBar(max =nc, style = 3)
        opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
        cl <- snow::makeCluster(threads)
        doSNOW::registerDoSNOW(cl)
        res <- foreach::foreach(i = 1:nc,.options.snow = opts) %dopar% {
                                  loop(i)
                                }
        snow::stopCluster(cl)
        gc()
        pcutils::del_ps("doSNOW","snow","foreach")
      }
      else {
        res <-lapply(1:nc, loop)
      }}

    simplify2array(res)->med.tax.cors

    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }
  else obs.exp.cors.mat=mycor

  diag(obs.exp.cors.mat) <- 0

  # Calculate connectedness by averaging positive and negative observed - expected correlations
  pn.mean=\(vector,mode="p"){
    if(mode=="p")pos.vals=vector[which(vector > 0)]
    else if(mode=="n") pos.vals=vector[which(vector < 0)]
    p.mean <- mean(pos.vals)
    if(length(pos.vals) == 0) p.mean <- 0
    return(p.mean)
  }
  connectedness.pos <- apply(obs.exp.cors.mat, 2, pn.mean,mode="p")
  connectedness.neg <- apply(obs.exp.cors.mat, 2, pn.mean,mode="n")

  # Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
  cohesion.pos <- rel.d %*% connectedness.pos
  cohesion.neg <- rel.d %*% connectedness.neg

  #### Combine vectors into one list and print
  output <- list(data.frame(neg=connectedness.neg, pos=connectedness.pos),
                 data.frame(neg=cohesion.neg, pos=cohesion.pos))

  names(output) <- c("Connectedness", "Cohesion")

  class(output)=c("cohesion",class(output))
  return(output)
}


#' Plot cohesion
#'
#' @param x `Cohesion()` result (cohesion object)
#' @param mode plot mode, 1~2
#' @param group group name in colnames(metadata)
#' @param metadata metadata
#' @param ... additional arguments for \code{\link[pcutils]{group_box}} (mode=1) or \code{\link[pcutils]{group_box}} (mode=2)
#'
#' @return a ggplot
#' @method plot cohesion
#' @exportS3Method
#' @rdname Cohesion
plot.cohesion=function(x,group,metadata,mode=1,...){
  cohesion_res=x
  if(mode==1)p=pcutils::stackplot(abs(t(cohesion_res$Cohesion)),group = group,metadata = metadata,...)
  if(mode==2){
    co=cohesion_res$Cohesion%>%dplyr::transmute(`neg:pos`=neg/pos)
    p=pcutils::group_box(co,group = group,metadata = metadata,p_value2 = TRUE,...)+
      ylab("neg:pos cohesion")+theme_bw()
  }
  p
}

#' Vulnerability
#'
#' @param go_ls a igraph object or igraph list.
#' @param threads threads
#'
#' @return a vector
#' @export
#' @description
#' \deqn{Vi=\frac{E-Ei}{E}}
#' E is the global efficiency and Ei is the global efficiency after the removal of the node i and its entire links.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' vulnerability(co_net)
#' }}
vulnerability=function(go_ls,threads=1){
  if("igraph"%in%class(go_ls))vulnerability_res=vul_max(go_ls,threads = threads)
  else {
    if(!"igraph"%in%class(go_ls[[1]]))stop("No igraph or igraph-list.")
    vulnerability_res=lapply(names(go_ls), \(i){
      vul_max(go_ls[[i]],threads = threads)
    })
    names(vulnerability_res)=names(go_ls)
  }
  class(vulnerability_res)=c("vulnerability",class(vulnerability_res))
  return(vulnerability_res)
}

vul_max<-function(go,threads=1){
  stopifnot(is.igraph(go))
  if(is.null(V(go)$name))V(go)$name=V(go)

  #parallel
  #main function
  loop=function (i) igraph::global_efficiency(igraph::delete_vertices(go, i))
  {
    if(threads>1){
      pcutils::lib_ps("foreach","doSNOW","snow")
      pb <- utils::txtProgressBar(max =length(V(go)$name), style = 3)
      opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
      cl <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cl)
      res <- foreach::foreach(i = V(go)$name,.options.snow = opts) %dopar% {
                                loop(i)
                              }
      snow::stopCluster(cl)
      gc()
      pcutils::del_ps("doSNOW","snow","foreach")
    }
    else {
      res <-lapply(V(go)$name, loop)
    }}
  #simplify method
  tmpv=do.call(c,res)

  (igraph::global_efficiency(go)-tmpv)/igraph::global_efficiency(go)
}

#' Plot vulnerability
#'
#' @param x `vulnerability()` result (vulnerability object)
#' @param ... add
#'
#' @return a ggplot
#' @method plot vulnerability
#' @exportS3Method
#' @rdname vulnerability
plot.vulnerability=function(x,...){
  vulnerability_res=x
  if(is.vector(vulnerability_res)){
    vulnerability_res=list(vulnerability=vulnerability_res)
  }
  pdat=data.frame(group=factor(names(vulnerability_res),levels = names(vulnerability_res)),
                  vulnerability=sapply(vulnerability_res,max))
  ggplot(pdat,aes(group,vulnerability))+
    geom_col(aes(fill=group),...)+
    geom_text(aes(group,vulnerability*1.05,label=round(vulnerability,3)))+
    theme_bw()
}


#' Robustness after remove 50%% nodes or some hubs, need the metanet contains "cor" attribute.
#'
#' @param go_ls a igraph object or igraph list.
#' @param keystone remove 70%% keystones instead of remove 50%% nodes (default: False)
#' @param reps simulation time
#' @param threads threads
#'
#' @export
#'
#' @examples
#' \donttest{
#' \dontrun{
#' robustness(co_net)->robustness_res
#' plot.robustness(robustness_res)
#' modu_dect(co_net) -> co_net_modu
#' zp_analyse(co_net_modu,mode = 2)->co_net_modu
#' robustness(co_net_modu,keystone=TRUE)->robustness_res
#' plot.robustness(robustness_res)
#' }}
robustness=function(go_ls,keystone=FALSE,reps=9,threads=1){
  if("igraph"%in%class(go_ls))robustness_res=robustness_in(go_ls,keystone=keystone,reps=reps,threads = threads)
  else {
    if(!"igraph"%in%class(go_ls[[1]]))stop("No igraph or igraph-list.")
    robustness_res=lapply(names(go_ls), \(i){
      tmp=robustness_in(go_ls[[i]],keystone=keystone,reps=reps,threads = threads)
      data.frame(tmp,group=i)
    })
    robustness_res=do.call(rbind,robustness_res)
    robustness_res$group=factor(robustness_res$group,levels = names(go_ls))
  }
  class(robustness_res)=c("robustness",class(robustness_res))
  return(robustness_res)
}

robustness_in<-function(go,keystone=FALSE,reps=9,threads=1){
  nodes <- length(V(go))
  floor(nodes * 0.5) -> del_i
  if(keystone){
    get_v(go)->tmp_v
    if(!"module"%in%colnames(tmp_v))stop("no modules, please `modu_dect()` first")
    if(!"roles"%in%colnames(tmp_v))stop("no roles, please `zp_analyse()` first")
    tmp_v%>%dplyr::filter(roles=="Module hubs")%>%dplyr::pull(name)->hubs
    floor(length(hubs)* 0.7) -> del_i
  }

  #parallel
  #main function
  loop=\(i){
    if(!keystone)remove_node <- sample(1:nodes, del_i)
    if(keystone)remove_node <- sample(hubs, del_i)

    dp <- igraph::delete.vertices(go, remove_node)
    dead_s="init"
    #calculated the abundance-weighted mean interaction strength (wMIS) of nodes,<=0 will dead
    #repeat until all >0
    while(length(dead_s)>0){
      get_e(dp)->edge_list
      edge_list%>%dplyr::select(from,cor)%>%
        rbind(.,dplyr::select(edge_list,to,cor)%>%dplyr::rename(from=to))%>%
        dplyr::group_by(from)%>%dplyr::summarise(w_degree=sum(cor))->w_degree
      w_degree%>%dplyr::filter(w_degree<=0)%>%dplyr::pull(from)->dead_s
      dp=igraph::delete.vertices(dp,dead_s)
    }
    tmp_ind=(MetaNet::net_par(dp,mode = "n")$n_index)
    data.frame(tmp_ind, reps = i)
  }
  {
    if(threads>1){
      pcutils::lib_ps("foreach","doSNOW","snow")
      pb <- utils::txtProgressBar(max =reps, style = 3)
      opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
      cl <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cl)
      res <- foreach::foreach(i = 1:reps,.options.snow = opts) %dopar% {
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
  del_i_indexs=do.call(rbind,res)
  return(del_i_indexs)
}

#' Plot robustness
#'
#' @param x `robustness()` result (robustness object)
#' @param indexes indexes selected to show
#' @param ... additional arguments for \code{\link[pcutils]{group_box}}
#'
#' @return a ggplot
#' @method plot robustness
#' @exportS3Method
#' @rdname robustness
plot.robustness=function(x,indexes="num_nodes",...){
  robustness_res=x
  if(!"group"%in%colnames(robustness_res))robustness_res$group="Network"
  pcutils::group_box(robustness_res[,indexes,drop=FALSE],group = "group",
                     metadata = robustness_res,...,facet = FALSE)+theme_bw()
}
