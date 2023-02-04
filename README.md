# MetaNet

network analysis for metagenomic data (big data)

## Install

`remotes::install_github('Asa12138/MetaNet',dependencies=T)`\
`library(MetaNet)`

## usage

    library(MetaNet)
    data(otutab)

    #1.Calculate spearman correlation for one t(otutab)
    t(otutab) -> totu  #492 OTUs
    system.time(c_net_cal(totu,threads = 4) -> corr)
    #2.347s, 
    system.time(psych::corr.test(totu,method = "spearman")->tmp)
    #>3min

    #2.Construct a network 
    c_net_build(corr,r_thres = 0.6,p_thres = 0.01) -> co_net
    #check the nodes information
    as.data.frame(vertex_attr(co_net))
    #annotation
    co_net <- c_net_set(co_net, t(otutab), taxonomy %>% select(Phylum))
    co_net <- anno_vertex(co_net, taxonomy)

    #3.plot
    c_net_plot(co_net)
    #change coordinate
    c_net_lay(co_net)->coors
    c_net_plot(co_net,coors)
    #ggplot style
    to.ggig(co_net,coors = coors)->ggig
    plot(ggig)
    #export to Gephi
    write_graph(co_net,file = "test.graphml",format = "graphml")
    #import from Gephi
    input_gephi("~/Documents/R/desert/desert_code/net_f1/Untitled.graphml")->gephi
    c_net_plot(gephi$go,coors = gephi$coors)

    #4.topology
    net_par(co_net)

    rand_net(co_net)
    fit_power(co_net)
    smallworldness(co_net)

    extract_sub_net(co_net,otutab,save_net = "../testnet")

    #5.modules detection
    modu_dect(co_net) -> co_net_modu
    graph.attributes(co_net_modu)
    modu_plot(co_net_modu, n_modu = 50)
    g_lay(co_net_modu,group ="module" ,zoom2 = 5,layout2 =nicely())->oridata
    modu_plot(co_net_modu,coors = oridata)
    g_lay_nice(co_net_modu,group ="module")->oridata
    modu_plot(co_net_modu,coors = oridata)

    zp_analyse(co_net_modu,mode = 2)->co_net_modu
    zp_plot(co_net_modu)

    #6.ecological
    robustness_test(co_net,step=8)->robust_res
    plot(robust_res,index="nat_connectivity",mode=2)+mytheme

    robustness(co_net)
    vulnerability(co_net)

    Cohesion(otutab[1:150,])->a
    stackplot(abs(t(a$Cohesion)),metadata,groupID = "Group")
