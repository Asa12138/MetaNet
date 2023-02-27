#========7.stability===========
#' Robustness_test for a network
#'
#' @param go a igraph object
#' @param partial how much percent vertexes be removed in total
#' @param step how many nodes be removed each time?
#' @param reps simulation number
#' @param threads threads
#'
#' @return dataframe(robustness class)
#' @export
#'
#' @examples
#' data("c_net")
#' robustness_test(co_net,step=4)->robust_res
#' class(robust_res)
#' plot(robust_res,index="ave_degree",mode=2)
robustness_test <- function(go, partial = 0.5,step=10,reps=9,threads=3) {
  lib_ps("igraph")
  cal_del<-\(go,partial,step,rep){
    nodes <- length(igraph::V(go))
    floor(nodes * partial) -> del_i
    del_i_indexs <- data.frame()
    sequ=seq(0,del_i,step)
    if(sequ[length(sequ)]<del_i)sequ=c(sequ,del_i)

    for (i in sequ) {
      # remove i nodes in the network
      remove_node <- sample(1:nodes, i)
      dp <- igraph::delete.vertices(go, remove_node)
      dp=igraph::delete.vertices(dp,igraph::V(dp)[igraph::degree(dp)==0])

      # calculate network parameters
      tmp_ind=(MetaNet::net_par(dp,mode = "n")$n_index)
      del_i_indexs <- rbind(
        del_i_indexs,
        data.frame(tmp_ind, i = i)
      )
    }
    return(data.frame(del_i_indexs,"rep"=rep))
  }

  #parallel
  #main function
  loop=function(i){
    cal_del(go,partial,step,i)
  }
  {
    if(threads>1){
      lib_ps("foreach","doSNOW")
      pb <- utils::txtProgressBar(max =reps, style = 3)
      opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
      cl <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cl)
      res <- foreach::foreach(i = 1:reps,.options.snow = opts) %dopar% {
                                loop(i)
                              }
      snow::stopCluster(cl)
      gc()
    }
    else {
      res <-lapply(1:reps, loop)
    }}
  #simplify method
  del_i_indexs=do.call(rbind,res)

  class(del_i_indexs)<-c("robustness",class(del_i_indexs))
  return(del_i_indexs)
}


#' Plot robustness
#'
#' @param robust_res robustness_test() result (robustness class)
#' @param mode plot mode
#' @param indexs indexs selected to show
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @return a ggplot
#' @method plot robustness
#' @export
#' @rdname robustness_test
plot.robustness<-function(robust_res, indexs = c("nat_connectivity", "ave_path_len", "ave_degree"),mode=1){
  lib_ps("reshape2")
  robust_res %>%
    dplyr::select(i, all_of(indexs)) %>%
    reshape2::melt(id.var = "i") -> pdat

  pdat%>%group_by(i,variable)%>%summarise(mean=mean(value),sd=sd(value),se=sd/sqrt(length(value)))->sdd

  if(mode==1){
    p <- ggplot(sdd, aes(x=i, y=mean,col=variable)) +
      geom_line()+
      geom_errorbar(data = sdd,aes(ymax=mean+se,ymin=mean-se))+
      #geom_smooth(se = FALSE,method = "loess",formula = 'y ~ x') +
      facet_wrap(. ~ variable, scales = "free")
  }
  if(mode==2){
    p<-ggplot(sdd,aes(x=i, y=mean,col=variable))+
      geom_point(size=0.2,alpha=0.4)+
      geom_smooth(se = F,method = "loess",formula = 'y ~ x')+
      ggpmisc::stat_poly_eq(
        aes(label = paste(after_stat(eq.label), after_stat(adj.rr.label), sep = '~~~~~')),
        formula = y ~ x,  parse = TRUE,label.x = "right",
        size = 3,
      )+
      facet_wrap(. ~ variable, scales = "free")
  }
  p=p+labs(x="Removed_nodes",y=NULL)+theme(legend.position = "none")
  return(p)
}

#' Cohesion calculate
#'
#' @param otutab otutab
#' @param reps iteration time
#' @param threads threads
#' @param mycor a correlation matrix you want to use, skip the null model build
#'
#' @return a list with two dataframe
#' @export
#' @references 1. Herren, C. M. & McMahon, K. Cohesion: a method for quantifying the connectivity of microbial communities. (2017) doi:10.1038/ismej.2017.91.
#' @examples
#' Cohesion(otutab[1:50,])->a
#' pctax::stackplot(abs(t(a$Cohesion)),metadata,groupID = "Group")
#' a$Cohesion%>%transmute(`neg:pos`=neg/pos)%>%pcutils::group_box(.,"Group",metadata)
Cohesion<-function(otutab,reps= 200,threads=4,mycor=NULL){
  d <- t(otutab)
  rel.d <- d / rowSums(d)

  if(is.null(mycor)){# Create observed correlation matrix
    cor.mat.true <- cor(rel.d)
    lib_ps("foreach")
    lib_ps("doSNOW")
    nc <- ncol(rel.d)
    pb <- txtProgressBar(max =nc, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    cl <- makeCluster(threads)
    registerDoSNOW(cl)
    a=foreach (which.taxon = 1:nc,.options.snow=opts) %dopar%{
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
    simplify2array(a)->med.tax.cors
    stopCluster(cl)

    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }
  else obs.exp.cors.mat=mycor

  diag(obs.exp.cors.mat) <- 0
  # Calculate connectedness by averaging positive and negative observed - expected correlations
  pn.mean=function(vector,mode="p"){
    if(mode=="p")pos.vals=vector[which(vector > 0)]
    else pos.vals=vector[which(vector < 0)]
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
  return(output)
}

#' Vulnerability
#'
#' @param go igraph
#' @param threads threads
#'
#' @return a vector
#' @export
#' @description
#' \deqn{Vi=\frac{E-Ei}{E}}
#' E is the global efficiency and Ei is the global efficiency after the removal of the node i and its entire links.
vulnerability<-function(go,threads=4){
  if(is.null(V(go)$name))V(go)$name=V(go)

  if(threads==1){
    lapply(V(go)$name,\(i)igraph::global_efficiency(igraph::delete_vertices(go,i)))->tmpls
    do.call(c,tmpls)->tmpv
  }
  else{
    #parallel
    lib_ps("foreach","doSNOW")
    nc=length(V(go)$name)
    pb <- txtProgressBar(max =nc, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    cl <- makeCluster(threads)
    registerDoSNOW(cl)
    tmpv <- foreach(i = V(go)$name,.combine = "c",
                    .options.snow = opts,
                    .packages = c("igraph"))%dopar%{
                      igraph::global_efficiency(igraph::delete_vertices(go,i))
                    }
    stopCluster(cl)
    gc()
  }
  (igraph::global_efficiency(go)-tmpv)/igraph::global_efficiency(go)
}

#' Robustness after remove 50\% nodes or some hubs
#'
#' @param go igraph
#' @param keystone remove 70\% keystone (default:False)
#' @param reps simulation time
#' @param threads threads
#'
#' @export
#'
#' @examples
#' robustness(co_net)
#' modu_dect(co_net) -> co_net_modu
#' zp_analyse(co_net_modu,mode = 2)->co_net_modu
#' robustness(co_net_modu,keystone=T)
robustness<-function(go,keystone=F,reps=9,threads=4){
  nodes <- length(V(go))
  floor(nodes * 0.5) -> del_i
  if(keystone){
    get_v(go)->tmp_v
    if(!"module"%in%colnames(tmp_v))stop("no modules, please modu_dect() first")
    tmp_v%>%filter(roles=="Module hubs")%>%pull(name)->hubs
    floor(length(hubs)* 0.7) -> del_i
  }

  #parallel
  lib_ps("foreach","doSNOW")
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  del_i_indexs<- foreach (i = 1:reps,.combine = "rbind",
                          .options.snow = opts,
                          .packages = c("igraph","dplyr")) %dopar%{
                            if(!keystone)remove_node <- sample(1:nodes, del_i)
                            if(keystone)remove_node <- sample(hubs, del_i)

                            dp <- igraph::delete.vertices(go, remove_node)
                            dead_s="init"
                            #calculated the abundance-weighted mean interaction strength (wMIS) of nodes,<=0 will dead
                            #repeat until all >0
                            while(length(dead_s)>0){
                              igraph::as_data_frame(dp)->edge_list
                              edge_list%>%select(from,cor)%>%rbind(.,select(edge_list,to,cor)%>%rename(from=to))%>%
                                group_by(from)%>%summarise(w_degree=sum(cor))->w_degree
                              w_degree%>%filter(w_degree<=0)%>%pull(from)->dead_s
                              dp=delete.vertices(dp,dead_s)
                            }

                            tmp_ind=(MetaNet::net_par(dp,mode = "n")$n_index)
                            data.frame(tmp_ind, reps = i)
                          }
  stopCluster(cl)
  return(del_i_indexs)
}
