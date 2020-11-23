library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)

?getOptInd

data(amgut1.filt)
#print(amgut_data)

depths <- rowSums(amgut1.filt)
print(depths)
print(class(depths))
print(dim(depths))

amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
print(d)
n <- nrow(amgut1.filt.cs)
print(n)
e <- d

### 2.合成数据
set.seed(10010)
graph <- SpiecEasi::make_graph('cluster', d, e)  ##从高斯图模型生成图拓扑结构。
##make_graph(method, D, e, enforce = TRUE, ...), D是生成图的方法，N是节点数，e是边数。 
Prec  <- graph2prec(graph) 
#?graph2prec #将图拓扑结构转化为相应的邻接矩阵，再转化为精度矩阵。
Cor   <- cov2cor(prec2cov(Prec))
#?prec2cov ##将精度矩阵转化为协方差矩阵。
#?cov2cor #cov2cor有效地将协方差矩阵缩放成相应的相关矩阵。

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)
#?synth_comm_from_counts(comm, mar = 2, distr, Sigma = cov(comm)) 
##从count数据拟合OTU边界的参数，模拟这些属性的新群落。
#comm 群落：计数矩阵；mar 群落数据矩阵的样本边界（1表示行；2表示列）。
#distr 拟合的分布（Zero-Inflated Negative Binomial Distribution）；sigma 协方差结构；n是样本数目。
print(X)
print(dim(X))

se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)
#?spiec.easi #data,没有归一化的OUT计数表。
# Applying data transformations...
# Selecting model with pulsar using stars...
# Fitting final estimate with mb...
# done
print(se)
print(se$est)
print(class(se)) #Model selection for penalized graphical models using the Stability Approach to Regularization Selection ('StARS')
print(se$est$path)
print(class(se$est$path))  ##列表
print(length(se$est$path))
print(se$est$path[1])  ##127 x 127 sparse Matrix of class "dgCMatrix"

huge::huge.roc(se$est$path, graph, verbose=FALSE)
stars.pr(getOptMerge(se), graph, verbose=FALSE)
# stars selected final network under: getRefit(se)
?huge::huge.roc  #huge.roc(path, theta, verbose = TRUE),根据真实图结构，绘制图路径的ROC曲线。
print(getOptMerge(se))
?getOptMerge  ##Get the optimal network, and related structures, when StARS is run.
?stars.pr     ##沿着stars的置信路径，绘制ROC或Precision-Recall曲线。




