#试试画循环通路网络
#======data======
#Nitrogen states

#' Load N-cycle data
#'
#' @return list
#' @export
#'
#' @references
#' 1. Tu, Q., Lin, L., Cheng, L., Deng, Y. & He, Z. NCycDB: a curated integrative database for fast and accurate metagenomic profiling of nitrogen cycling genes. Bioinformatics 35, 1040–1048 (2019).
#' 2. Kuypers, M. M. M., Marchant, H. K. & Kartal, B. The microbial nitrogen-cycling network. Nat Rev Microbiol 16, 263–276 (2018).
load_N_data=function(){
  lib_ps("dplyr",library = FALSE)
N_df=data.frame(X=c((+5):(-3),-3),Y=c(rep(1,9),3),
                `oxidation state`=c("+5","+4","+3","+2","+1","0","-1","-2","-3","-3"),
                N=c("NO3-",NA,"NO2-","NO","N2O","N2","NH2OH","N2H4","NH4+","organic_N"),
                N_name=c("paste(NO[3]^'-')",NA,"NO[2]^'-'","NO","paste(N[2],O)","N[2]","paste(NH[2],OH)","paste(N[2],H[4])","NH[4]^'+'","organic_N"),
                check.names = FALSE)
#Pathways
N_path1=tibble::tribble(
        ~from,         ~to,              ~Pathway,              ~Pathway2, ~curvature, ~up_down,
       "NO3-",      "NO2-",     "Denitrification",     "Denitrification1",         1L,   "down",
       "NO3-",      "NO2-",                "ANRA",                "ANRA1",         7L,     "up",
       "NO3-",      "NO2-",                "DNRA",                "DNRA1",         6L,     "up",
       "NO2-",      "NO3-",       "Nitrification",       "Nitrification3",         1L,     "up",
       "NO2-",        "NO",     "Denitrification",     "Denitrification2",         1L,   "down",
       "NO2-",      "NH4+",                "ANRA",                "ANRA2",         3L,     "up",
       "NO2-",      "NH4+",                "DNRA",                "DNRA2",         2L,     "up",
         "NO",       "N2O",     "Denitrification",     "Denitrification3",         1L,   "down",
         "NO",      "co-N",             "Anammox",             "Anammox1",         8L,   "down",
        "N2O",        "N2",     "Denitrification",     "Denitrification4",         1L,   "down",
         "N2",      "NH4+",          "N-fixation",           "N-fixation",         1L,   "down",
       "N2H4",        "N2",             "Anammox",             "Anammox2",         4L,   "down",
      "NH2OH",      "NO2-",       "Nitrification",       "Nitrification2",         1L,     "up",
       "NH4+",     "NH2OH",       "Nitrification",       "Nitrification1",         1L,     "up",
       "NH4+",      "co-N",             "Anammox",             "Anammox1",         5L,   "down",
       "NH4+", "organic_N",   "Organic-synthesis",    "Organic-synthesis",         0L,   "down",
  "organic_N",      "NO2-", "Organic-degradation", "Organic-degradation1",         4L,   "down",
  "organic_N",      "NH4+", "Organic-degradation",  "Organic-degradation",         0L,   "down",
       "co-N",      "N2H4",             "Anammox",             "Anammox1",         0L,   "down"
  )%>%as.data.frame()

#Nitrogen metabolism related genes
{
  N_genes=tibble::tribble(
                  ~Pathway2, ~`Gene_families`,    ~X,    ~Y,
           "Nitrification1",               "amoA_A",  -1.5,   0.8,
           "Nitrification1",               "amoB_A",    -2,   0.7,
           "Nitrification1",               "amoC_A",  -2.5,   0.8,
           "Nitrification1",               "amoA_B",  -1.5,   0.6,
           "Nitrification1",               "amoB_B",    -2,   0.5,
           "Nitrification1",               "amoC_B",  -2.5,   0.6,
           "Nitrification2",                  "hao",     1,   0.4,
           "Nitrification3",                 "nxrA",   4.2,   0.8,
           "Nitrification3",                 "nxrB",   3.8,   0.8,
         "Denitrification1",                 "napA",   4.4,   1.9,
         "Denitrification1",                 "napB",     4,   1.9,
         "Denitrification1",                 "napC",   3.6,   1.9,
         "Denitrification1",                 "narG",   4.4,   1.7,
         "Denitrification1",                 "narH",     4,   1.7,
         "Denitrification1",                 "narJ",   3.6,   1.7,
         "Denitrification1",                 "narI",   4.4,   1.5,
         "Denitrification1",                 "narZ",     4,   1.5,
         "Denitrification1",                 "narY",   3.6,   1.5,
         "Denitrification1",                 "narV",   4.2,   1.3,
         "Denitrification1",                 "narW",   3.8,   1.3,
         "Denitrification2",                 "nirK",   2.7,   1.3,
         "Denitrification2",                 "nirS",   2.3,   1.3,
         "Denitrification3",                 "norB",  1.65,   1.3,
         "Denitrification3",                 "norC",   1.3,   1.3,
         "Denitrification4",                 "nosZ",   0.5,   1.3,
                    "ANRA1",                 "nasA",   4.4, -0.55,
                    "ANRA1",                 "nasB",     4, -0.55,
                    "ANRA1",                   "NR",   3.7, -0.55,
                    "ANRA1",                 "narB",   4.2, -0.75,
                    "ANRA1",                 "narC",   3.8, -0.75,
                    "ANRA2",                 "nirA",     0,  -0.8,
                    "DNRA1",                 "napA",   4.4,   0.4,
                    "DNRA1",                 "napB",     4,   0.4,
                    "DNRA1",                 "napC",   3.6,   0.4,
                    "DNRA1",                 "narG",   4.3,   0.2,
                    "DNRA1",                 "narH",     4,   0.2,
                    "DNRA1",                 "narJ",  3.72,   0.2,
                    "DNRA1",                 "narI",   4.3,     0,
                    "DNRA1",                 "narZ",     4,     0,
                    "DNRA1",                 "narY",  3.72,     0,
                    "DNRA1",                 "narV",   4.2,  -0.2,
                    "DNRA1",                 "narW",   3.8,  -0.2,
                    "DNRA2",                 "nirB",   0.3,  -0.1,
                    "DNRA2",                 "nirD",     0,  -0.1,
                    "DNRA2",                 "nrfA",  -0.3,  -0.1,
                    "DNRA2",                 "nrfB",   0.3,  -0.3,
                    "DNRA2",                 "nrfC",     0,  -0.3,
                    "DNRA2",                 "nrfD",  -0.3,  -0.3,
               "N-fixation",                 "anfG",  -1.2,   1.6,
               "N-fixation",                 "nifD", -1.55,   1.6,
               "N-fixation",                 "nifH", -1.85,   1.6,
               "N-fixation",                 "nifK",  -1.4,   1.4,
               "N-fixation",                 "nifW",  -1.7,   1.4,
                 "Anammox1",                 "hzsA",  -0.8,  2.25,
                 "Anammox1",                 "hzsB",  -1.2,  2.25,
                 "Anammox1",                 "hzsC",  -0.8,     2,
                 "Anammox2",                  "hdh",  -1.2,   1.2,
                 "Anammox2",                  "hzo",  -0.9,   1.2,
        "Organic-synthesis",                 "glsA",  -3.3,   2.8,
        "Organic-synthesis",                 "glnA",  -3.7,   2.8,
        "Organic-synthesis",                 "asnB",  -3.3,  2.65,
        "Organic-synthesis",                 "ansB",  -3.7,  2.65,
        "Organic-synthesis",            "gs_K00264", -3.45,   2.5,
        "Organic-synthesis",            "gs_K00265", -3.45,  2.35,
        "Organic-synthesis",            "gs_K00266", -3.45,   2.2,
        "Organic-synthesis",            "gs_K00284", -3.45,  2.05,
     "Organic-degradation1",                  "nao",  1.35,     2,
     "Organic-degradation1",                  "nmo",     1,   2.1,
      "Organic-degradation",                 "ureA", -3.25,  1.75,
      "Organic-degradation",                 "ureB", -3.55,  1.75,
      "Organic-degradation",                 "ureC", -3.88,  1.75,
      "Organic-degradation",           "gdh_K00260",  -3.5,   1.6,
      "Organic-degradation",           "gdh_K00261",  -3.5,  1.45,
      "Organic-degradation",           "gdh_K00262",  -3.5,   1.3,
      "Organic-degradation",           "gdh_K15371",  -3.5,  1.15
     )
  N_genes=dplyr::left_join(N_genes,N_path1%>%dplyr::distinct(Pathway,Pathway2))
}

N_path1%>%dplyr::left_join(.,rbind(N_df,list(-2,2,NA,"co-N",NA))[,c("N","X","Y")],by=c("from"="N"))%>%dplyr::rename(X1=X,Y1=Y)%>%
  dplyr::left_join(.,rbind(N_df,list(-2,2,NA,"co-N",NA))[,c("N","X","Y")],by=c("to"="N"))%>%dplyr::rename(X2=X,Y2=Y)->N_path

N_path%>%dplyr::mutate(Y1=ifelse(up_down=="down",Y1+0.1,ifelse(up_down=="up",Y1-0.1,Y1)))%>%
  dplyr::mutate(Y2=ifelse(up_down=="down",Y2+0.1,ifelse(up_down=="up",Y2-0.1,Y2)))->N_path

N_path[which(N_path$to=="organic_N"),c("X2","Y2")]=c(-3.1,2.9)
N_path[which(N_path$from=="organic_N"),c("X1","Y1")]=rep(c(-2.9,2.9),each=2)

return(list(N_path=N_path,N_genes=N_genes,N_df=N_df))
}
#===ggplot====

#' Plot the N-cycling pathway and genes
#'
#' @param my_N_genes dataframe, "Gene_families","type" should in colnames of my_N_genes
#' @param just_diff logical, just plot the different genes?
#' @param path_col colors of pathways
#' @param type_col colors of types
#' @param fill_alpha alpha, default 0.5
#' @param arrow_size arrow_size, default 0.1
#' @param line_width line_width, default 1
#' @param title title, default "Nitrogen cycling"
#' @param legend.position default c(0.85,0.15)
#'
#' @return ggplot
#' @export
#'
#' @examples
#' N_data=load_N_data()
#' my_N_genes=data.frame(
#'     `Gene_families`=sample(N_data$N_genes$Gene_families,10,replace = FALSE),
#'     change=rnorm(10),check.names = FALSE)
#' my_N_genes=dplyr::mutate(my_N_genes,type=ifelse(change>0,"up",ifelse(change<0,"down","none")))
#' plot_N_cycle(my_N_genes,just_diff = FALSE,fill_alpha=0.2)
#' \donttest{
#' \dontrun{
#' ggsave(filename = "test.pdf",width = 14,height = 10)
#' }}
plot_N_cycle=function(my_N_genes=NULL,just_diff=FALSE,path_col=NULL,type_col=c(up="red",down="blue",none=NA),fill_alpha=0.5,
                      arrow_size=0.1,line_width=1,title="Nitrogen cycling",legend.position = c(0.85,0.15)){
if(!is.null(my_N_genes))if(!all(c("Gene_families","type")%in%colnames(my_N_genes)))stop('"Gene_families","type" should in colnames of my_N_genes')

N_data=load_N_data()
N_path=N_data$N_path
N_genes=N_data$N_genes
N_df=N_data$N_df

if(is.null(path_col)){
  path_col=c("Denitrification"="#EB443E","DNRA"="#A6830F","Anammox"="#62CB68","ANRA"="#8AC0DD","N-fixation"="#F2DF4E",
           "Organic-degradation"="#F5AB80","Organic-synthesis"="#8B71DD","Nitrification"="#3C57AA")
  }

p1=ggplot()+geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = Pathway),linewidth=line_width,
                        data = filter(N_path,curvature==0),curvature=0,arrow=arrow(length = unit(arrow_size, "inches"),type="closed"))+
    geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = Pathway),linewidth=line_width,
               data = filter(N_path,curvature==1),curvature=-0.3,arrow=arrow(length = unit(arrow_size, "inches"),type="closed"))+
    geom_curve(aes(x = X1, y = Y1-0.1, xend = X2, yend = Y2-0.1,color = Pathway),linewidth=line_width,
               data = filter(N_path,curvature==2),curvature=0.45,arrow=arrow(length = unit(arrow_size, "inches"),type="closed"))+
    geom_curve(aes(x = X1, y = Y1-0.15, xend = X2, yend = Y2-0.15,color = Pathway),linewidth=line_width,
               data = filter(N_path,curvature==3),curvature=0.6,arrow=arrow(length = unit(arrow_size, "inches"),type="closed"))+
    geom_curve(aes(x = X1, y = Y1-0.1, xend = X2, yend = Y2-0.1,color = Pathway),linewidth=line_width,
               data = filter(N_path,curvature==6),curvature=0.8,arrow=arrow(length = unit(arrow_size, "inches"),type="closed"))+
    geom_curve(aes(x = X1, y = Y1-0.15, xend = X2, yend = Y2-0.15,color = Pathway),linewidth=line_width,
               data = filter(N_path,curvature==7),curvature=1.3,arrow=arrow(length = unit(arrow_size, "inches"),type="closed"))+
    geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = Pathway),linewidth=line_width,
               data = filter(N_path,curvature==4),curvature=0.2,arrow=arrow(length = unit(arrow_size, "inches"),type="closed"))+
    geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = Pathway),linewidth=line_width,
               data = filter(N_path,curvature==5),curvature=0.2)+
    geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = Pathway),linewidth=line_width,
               data = filter(N_path,curvature==8),curvature=-0.2)+
    scale_color_manual(name="Pathway",values = path_col,guide=guide_legend(nrow=4))+ylim(-2,3.5)+
    mytheme+scale_x_reverse(breaks =(+5):(-3),labels = c(N_df$`oxidation state`[1:9]))

p2=p1+geom_text(aes(x=X,y=Y,label=paste('bold(',N_name,')')),data = N_df,parse=TRUE,size=5)+
    labs(x="Oxidation state",y=NULL)+
    theme(legend.position = "right", panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), axis.line.y = element_blank(),
          axis.ticks.y = element_blank(), axis.text.y = element_blank())

message("recommend ggsave(width = 14,height = 10)")

if(is.null(my_N_genes)){
  p3=p2+geom_text(aes(x=X,y=Y,label=`Gene_families`,color=Pathway),data =N_genes,show.legend = FALSE)
}
else{
  if(!just_diff) {my_N_genes_change=dplyr::left_join(N_genes,my_N_genes)}
  else {my_N_genes_change=dplyr::left_join(my_N_genes,N_genes)}

  p3=p2+ggnewscale::new_scale_color()+
    geom_label(aes(x=X,y=Y,label=`Gene_families`,color=Pathway,fill=type),alpha=fill_alpha,data =my_N_genes_change)+
    scale_color_manual(values = path_col,guide=guide_none())
  if(!is.numeric(my_N_genes_change%>%pull(type)))p3=p3+scale_fill_manual(values = type_col,na.value = NA,guide=guide_legend(direction = "horizontal"))
  else p3=p3+scale_fill_gradient(low = "white",high = "red",na.value = NA,guide=guide_colorbar(direction = "horizontal"))
}

p3+coord_fixed(ratio = 1.02)+
  geom_text(data = data.frame(x = 4, y = 3.5, label = title),
            mapping = aes(x = x, y = y, label = label),
            size = 10, fontface = 2, inherit.aes = FALSE)+
  theme(legend.position = legend.position,legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
}

