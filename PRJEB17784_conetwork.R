setwd("F:\\Zhaolab2020\\BioNetwork\\MicrobiomeNetwork\\MNDnetwork\\PRJEB17784\\co-network")

####Step 1 - Installation and biom format conversion####
#install.packages("devtools")
library(devtools)
#install_github('zdk123/SpiecEasi', host = "api.github.com", dependencies = TRUE)
library(SpiecEasi)

#install.packages("BiocManager")
#BiocManager::install("phyloseq")
library(phyloseq)

#install_github("hallucigenia-sparsa/seqtime") 
library(seqtime)

####Step 2 - Data import and preprocessing####
#明确arctic_soils_json.biom的路径
input.path="F:\\Zhaolab2020\\BioNetwork\\MicrobiomeNetwork\\MNDnetwork\\PRJEB17784\\co-network"
output.path="F:\\Zhaolab2020\\BioNetwork\\MicrobiomeNetwork\\MNDnetwork\\PRJEB17784\\co-network\\output"
species.path=file.path(input.path,"Control_filtered_species.txt")

library(SpiecEasi)
data0=read.table(species.path, head=TRUE, row.names=1) #species abundance table, sample ids as rownames
# print(data0)
# print(rownames(data0))
# print(colnames(data0))
# print(dim(data0))
# data<-t(data0)
# print(rownames(data))
# print(colnames(data))
# print(dim(data))

pargs=list(seed=666,ncores=1)
# print(pargs)
# print(class(pargs))
# ?list

spiec = spiec.easi(as.matrix(data), method='glasso',sel.criterion = "stars",lambda.min.ratio=1e-3, nlambda=30, pulsar.params=pargs)
# print(spiec)
# print(getRefit(spiec))  ###getOptNet / getRefit: the optimal (StARS-refit) network
# ?getRefit

library(Matrix)
ig.spiec <- adj2igraph(getRefit(spiec))

library(igraph)
## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(as.matrix(data), 1))+6
am.coord <- layout.fruchterman.reingold(ig.spiec)

par(mfrow=c(1,1))
plot(ig.spiec, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Control")

V(ig.spiec)
V(ig.spiec)$label <- seq(colnames(data))
E(ig.spiec)
get.edgelist(ig.spiec)
print(ig.spiec)


diameter(ig.spiec)
transitivity(ig.spiec)
mean_distance(ig.spiec)
mean(centr_betw(ig.spiec)$res)
#?mean_distance
fast_greedy <- cluster_fast_greedy(ig.spiec)
#modularity(fast_greedy)
modularity(ig.spiec, membership(fast_greedy))

pr <- page.rank(ig.spiec)
key_species <- data.frame(Object=1:ncol(data), PageRank=pr$vector)
print(key_species)
key_species$Species <- colnames(data)
print(key_species)
library(dplyr)
key_species_ranked <- arrange(key_species, desc(PageRank))      #按PageRank降序排列
#print(key_species_ranked[1:10, c(1,3,2)])
write.csv(key_species_ranked[, c(1,3,2)], file = "02spiec.easi_key_species_ranked_Control.csv", row.names = TRUE)


###<---------------------------------------------------------------------------->###
cor=cov2cor(apply(getOptCov(spiec),2,as.numeric))
###将相关矩阵转化为协方差矩阵。
row.names(cor)=colnames(data) ##协方差矩阵的行名=数据的列名
colnames(cor)=colnames(data)  ##协方差矩阵的列名=数据的列名
cor[lower.tri(cor, diag = T)]=NA
#?lower.tri  #lower.tri(x, diag = FALSE), X是matrix或者其它R对象。diag是逻辑变量，是否应该包括对角线。
ind = which(is.na(cor)==F,arr.ind = T)  ##返回协方差矩阵中不为NA的值的index。
#?which  #Give the TRUE indices of a logical object, allowing for array indices.给出逻辑对象的真索引，允许使用数组索引。
print(ind)
print(colnames(cor))
print(colnames(cor)[ind[,1]])
result = data.frame(start = colnames(cor)[ind[,1]],end = colnames(cor)[ind[,2]],Cor = cor[ind],stringsAsFactors = F)

print(result)
write.csv(result, file = "02spiec.easi_result_Control.csv", row.names = TRUE)


# ####新的测试
# se.mb.mNetwork <- spiec.easi(as.matrix(data), method='mb', lambda.min.ratio=1e-2,
#                           nlambda=20, pulsar.params=list(rep.num=50))
# se.gl.mNetwork <- spiec.easi(as.matrix(data), method='glasso', lambda.min.ratio=1e-2,
#                           nlambda=20, pulsar.params=list(rep.num=50))
# sparcc.mNetwork <- sparcc(as.matrix(data))
# ## Define arbitrary threshold for SparCC correlation matrix for the graph
# sparcc.graph <- abs(sparcc.mNetwork$Cor) >= 0.1
# diag(sparcc.graph) <- 0
# library(Matrix)
# sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
# ## Create igraph objects
# ig.mb     <- adj2igraph(getRefit(se.mb.mNetwork))
# ig.gl     <- adj2igraph(getRefit(se.gl.mNetwork))
# ig.sparcc <- adj2igraph(sparcc.graph)
# 
# library(igraph)
# ##统计图的特性########################################################################################
# #参考：https://mattzheng.blog.csdn.net/article/details/51444536
# 
# V(ig.mb)  ##+ 89/89 vertices, named, from 7484ad2
# E(ig.mb)  ##+ 135/135 edges from 7484ad2 (vertex names)
# ig.mb.clusters <- clusters(ig.mb,mode="weak")
# print(ig.mb.clusters)
# ig.mb.member<-walktrap.community(ig.mb,weights=E(ig.mb)$weight,step=4)
# print(ig.mb.member)  ##IGRAPH clustering walktrap, groups: 16, mod: 0.61
# 
# ##基于中间中心度的社团发现
# ig.mb.edge.betweenness.member<-edge.betweenness.community(ig.mb,weight=E(ig.mb)$weight,directed=F)
# print(ig.mb.edge.betweenness.member) #IGRAPH clustering edge betweenness, groups: 13, mod: 0.63
# 
# 
# V(ig.gl)   ##+ 89/89 vertices, named, from 74d97e7
# E(ig.gl)   ##+ 211/211 edges from 74d97e7 (vertex names):
# 
# V(ig.sparcc)   ##+ 89/89 vertices, named, from 75215de
# E(ig.sparcc)   ##+ 184/184 edges from 75215de (vertex names)
# #####################################################################################################
# 
# ## set size of vertex proportional to clr-mean
# vsize    <- rowMeans(clr(as.matrix(data), 1))+10
# am.coord <- layout.fruchterman.reingold(ig.mb)
# 
# par(mfrow=c(1,3))
# plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
# plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
# plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")