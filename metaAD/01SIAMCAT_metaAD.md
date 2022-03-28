加载包：**

```R
#####models汇总#####
#update: 2021.7.7 by ljx

#加载包
library(tidyverse)
library(SIAMCAT)
library(randomForest)
library(dplyr)
library(tidyr)
library(reshape2)  
library(ggplot2)
library(cowplot)
library(RColorBrewer)
```

##### 2.6.1 数据填充，缺失补0

```R
#####一、数据的填充：缺失的补充为0即可#####
#input:
#' @param feat_list去除Unclassified数据后的，行为样本，列为特征(已经进行归一化)

my_pair.table=function(feat_list){
  library(tidyverse)
  
  #1.先进行其余特征补充
  feat_name=lapply(feat_list,function(data){colnames(data)})
  feat_ID=feat_name%>%unlist()%>%as.vector()%>%unique()#所有特征
  b=lapply(feat_list, function(data){
    add=setdiff(feat_ID,colnames(data))
  })
  
  my_add=function(b,feat){
    if(length(b)!=0){
      b=as.character(b)#去除一些属性标签
      x_add=data.frame(matrix(0,dim(feat)[1],length(b)))
      colnames(x_add)=b#所有没识别出来的全部为0
      xtest=cbind(feat,x_add)#合并在一起啦，按列合并(行名仍保留)：样本名
    }else{
      xtest=feat
    }
    return(xtest)
  }
  
  #补充所有特征
  data_add=list()
  for (i in 1:length(feat_list)) {
    data_add[[i]]=my_add(b[[i]],feat_list[[i]])
  }
  #列表名需要单独设置
  names(data_add)=names(feat_list)
  
  ##2.预处理：将列向量按照feature_ID重新排序：F1,F2,...
  feat_seq=sort(feat_ID)
  data_add=lapply(data_add, function(data){new_data=data[,c(feat_seq)]})
  
  return(data_add)
}
```

-   取特征集合部分可参考：https://bioconductor.org/packages/release/bioc/vignettes/SIAMCAT/inst/doc/SIAMCAT_ml_pitfalls.html



##### 2.6.2 单个项目的训练

```R
#####二、单个项目模型的训练#####
#siamcat流程需要feat、meta独立输入，行名对应顺序排列
#input：
#' @param method是调用siamcat包的统计方法：默认lasso(字符串)
#' @param num.folds=5,num.resample=3(默认：3次5折)

#output：
#model_back返回的是siamcat对象，aucinfo只有一个auc值
my_siamcat=function(train_feat,train_meta,method="lasso",num.folds=5,num.resample=3){
  library(SIAMCAT)
  library('randomForest')
  #1.读为siamcat对象（phyloseq对象）
  siamcat.train <- siamcat(feat=t(train_feat), meta=as.data.frame(train_meta),
                           label='Group', case='Case')
  #注意此函数必须是metadata含两列信息以上，不然就会报错
  
  #2.丰度过滤
  siamcat.train <- filter.features(
    siamcat.train,
    filter.method = 'abundance',
    cutoff = 0.001,
    rm.unmapped = TRUE,
    verbose=2
  )
  #默认去除所有unmapping的特征(但在这之前就采取了计算归一化)
  #每个特征在所有样本中都没超过0.001的话就去除该特征
  #(比较每个特征在该样本中最大丰度即可)
  
  #3.标准化(使得原始很多0变为非0--sparity=0%)
  siamcat.train <- normalize.features(
    siamcat.train,
    norm.method = "log.std",
    norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1),
    verbose = 2
  )
  
  #4.交叉验证分割数据
  siamcat.train <-  create.data.split(
    siamcat.train,
    num.folds = num.folds,
    num.resample = num.resample
  )
  
  #5.训练模型
  siamcat.train<- train.model(
    siamcat.train,
    method = method
  )
  ## Trained lasso models successfully.
  
  #6.模型评估
  #选择模型后：用于全部数据集
  siamcat.train <- make.predictions(siamcat.train)
  siamcat.train <-  evaluate.predictions(siamcat.train)
  
  auc=as.numeric(eval_data(siamcat.train)$auroc)
  cat(auc,'\n')
  #可以得到auc值(整个训练项目的值)
  
  list(model_back=siamcat.train,aucinfo=auc)
}
```

-   **model_back**中存储的是siamcat.train对象。

    

##### 2.6.3 外部模型的交叉验证

```R
####三、外部模型交叉验证#####
#' @param plot=F不画出auc交叉图
#output：
#model_back返回的是siamcat对象，aucinfo只有一个auc值
#与my_siamcat输出结构一致
my_external.vali_siamcat=function(siamcat.obj,test_feat,test_meta,plot=F){
  library(SIAMCAT)
  siamcat.train=siamcat.obj$model_back
  
  #1.读为siamcat对象（phyloseq对象）
  siamcat.test <- siamcat(feat=t(test_feat), meta=as.data.frame(test_meta),
                          label='Group', case='Case')
  siamcat.test <- normalize.features(siamcat.test,
                                     norm.param=norm_params(siamcat.train),
                                     feature.type='original',
                                     verbose = 2)
  #2.代入模型
  siamcat.test <- make.predictions(
    siamcat = siamcat.train,
    siamcat.holdout = siamcat.test,
    normalize.holdout = FALSE)
  #normalize.holdout = FALSE是因为前方已经normlize
  siamcat.test <- evaluate.predictions(siamcat.test)
  
  auc=as.numeric(eval_data(siamcat.test)$auroc)
  cat(auc,'\n')
  #可以得到auc值(整个训练项目的值)
  
  #3.画该两个项目的auc.prc图像
  if(plot==T){
    model.evaluation.plot('train'=siamcat.train,
                          'test'=siamcat.test,
                          colours=c('dimgrey', 'orange'))
  }
  
  list(model_back=siamcat.test,aucinfo=auc)
}
```

##### 2.6.4 筛选特征训练模型

```R
#####四、筛选特征训练模型#####
#与my_siamcat输出结构一致
#' @param siamcat是train后的对象，默认都是normalized后的
#' @param consens.thres是筛选变量percentage的阈值(lasso:0.5,指重复建模后系数非0比例)
#' @param max.show控制最大选择特征数(与上面阈值percentage选出个数比较，取min)

my_siamcat.sel=function(siamcat,
                        consens.thres=0.5,
                        max.show=20,
                        verbose = 1,method='lasso'){
  #1.定义函数
  library(SIAMCAT)
  library('randomForest')
  siamcat.select.features <-function(
    siamcat,
    consens.thres,
    max.show,
    verbose = 1) {
    
    library(SIAMCAT)
    library('randomForest')
    
    label <- label(siamcat)
    model.type <- model_type(siamcat)
    feature.type <- feature_type(siamcat)
    models <- models(siamcat)
    feature.weights <- feature_weights(siamcat)
    weight.matrix <- weight_matrix(siamcat)
    if (verbose > 2) message("+ model.interpretation.select.features")
    
    if (model.type != "randomForest") {
      sel.idx = which(feature.weights$percentage > consens.thres)
      median.sorted.features <-
        sort(feature.weights$median.rel.weight[sel.idx],
             decreasing = TRUE,
             index.return = TRUE)
      
      if (length(sel.idx) > max.show) {
        warning(paste0("WARNING: restricting amount of features",
                       " to be plotted to ", max.show))
        #当根据percenage筛选Top feature超过限度时，根据相对权重中位数排序
        #选出其中的top（绝对值：从大到小）
        median.sorted.features.abs <- sort(
          abs(feature.weights$median.rel.weight),
          decreasing = TRUE,
          index.return = TRUE)
        idx <- head(median.sorted.features.abs$ix, n = max.show)
        median.sorted.features <- sort(
          feature.weights$mean.rel.weight[idx],
          decreasing = TRUE,
          index.return = TRUE)
        sel.idx <- idx[median.sorted.features$ix]
      } else {
        warning(paste0("OK:An amount of features",
                       " to be plotted to ", length(sel.idx)))
        sel.idx = sel.idx[median.sorted.features$ix]
      }
    } else {
      #随机森林中利用相对weight进行排序，与重要性有啥区别与联系？
      median.sorted.features <-
        sort(feature.weights$median.rel.weight,
             decreasing = FALSE,
             index.return = TRUE)
      # take the feature with median higher than consens.threshold
      consens.thres=0#修改:consens.thres=0.5太高了！
      #修改后相当于全选
      sel.idx <-
        median.sorted.features$ix[which(median.sorted.features$x >=
                                          consens.thres)]
      
      sel.idx <- tail(sel.idx, n = max.show)
      warning(paste0("WARNING:Rf restricting amount of features",
                     " to be plotted to ", max.show))
    }
    
    if (verbose > 2)
      message(paste(
        "+++ generating plot for a model with",
        length(sel.idx),
        "selected features"
      ))
    if (verbose > 2)
      message("+ finished model.interpretation.select.features")
    
    feat.norm<-get.norm_feat.matrix(siamcat)[sel.idx,]
    norm_feat(siamcat)[[1]] <- feat.norm
    
    #得到新的norm_feat即可
    return(siamcat)
  }
  siamcat.train <- siamcat.select.features(siamcat,consens.thres,max.show,verbose)
  
  #2.将特征选择+训练好的siamcat重新建模
  #训练模型
  siamcat.train<- train.model(
    siamcat.train,
    method = method
  )
  
  #3.模型评估
  #3.1选择模型后：用于test数据集
  siamcat.train <- make.predictions(siamcat.train)
  
  #3.2得到一个eval_data数据
  siamcat.train <-  evaluate.predictions(siamcat.train)
  auc=as.numeric(eval_data(siamcat.train)$auroc)
  cat(auc,'\n')
  #可以得到auc值(整个训练项目的值)
  
  list(model_back=siamcat.train,aucinfo=auc)
}
```

##### 2.6.5 可视化

```
#####五、可视化#####
#heatmap矩阵图（一种疾病多项目间的交叉验证）+ barplot
#(一)所有数据训练，储存矩阵图的数据格式#####
#输出：models是所有特征模型，models_top是top特征(如果top=null,那就是models)
#' @param top是否选择top特征，前提是之前已经训练过models
#' @param models是用所有特征进行训练的列表
#' @param label是指将method=label作为因子进行可视化
#' @param new：新增数据的project_id(此时models非空，把新增数据的位置NA即可)
#' 注意lodo数据集如若有更新，则是全部重新计算，无需new
#' @param auc.matrix=NULL是调出原来auc的交叉验证结果（连接mysql）

my_result=function(feat,meta,method,label,top=NULL,models=NULL,num.folds=5, num.resample=3,new=NULL,auc.matrix=NULL){
  #update5.23加入一个选择交互new(新增数据)：新增数据的project_id
  #update5.29加入输出交叉验证的模型，以便于画出ROC曲线，cross_models
  #注意此时cross_models只保存新增数据的交叉验证模型
  # feat = datas_AD
  # meta = meta_list
  # method = 'lasso'
  # label ='lasso_top20'
  # models = results_AD$models
  # top = 20
  # models=NULL
  # num.folds=5
  # num.resample=3
  # new=NULL
  # auc.matrix=NULL
  if(!is.null(new)){
    message('you add new datas to train new models!')
    #1.除新增数据未训练模型外，其余models可以直接调用
    #可以是多个新增数据
    if(is.null(models)){stop('models.obj can not be NULL!')}
    if(is.null(auc.matrix)){stop('auc.matrix.obj can not be NULL!')}
    for (i in new) {
      if(is.na(models[[i]])){
        #若新增数据集的models对象为NA时，需要重新训练
        models[[i]]=my_siamcat(train_feat = feat[[i]],train_meta =meta[[i]],method,num.folds,num.resample)
        
      }
    }
    
    models_top=models#与下面代码统一处理
    #注意：若为外部调用，此时的models要么是all.feat/top.feat(均作为models输入)
    if(!is.null(top)){
      for(i in new){
        models_top[[i]]=my_siamcat.sel(siamcat=models_top[[i]]$model_back,
                                       consens.thres=0.5,
                                       max.show=top,
                                       verbose = 1,method=method)
      }
      
    }else{
      message("you don't choose top features!")
    }
    
    #2.组成一个auc矩阵
    task_id=names(feat)
    add.id=c()#单独增加的行，列
    for (k in new) {
      add.id=c(add.id,which(task_id==k))
    }
    
    cross.siamcat <- matrix(NA,length(task_id),length(task_id))#储存AUC矩阵
    cross_models=list()#train-test命名列表
    for (i in 1:length(task_id)) {
      for (j in 1:length(task_id)) {
        if(i==j){
          cross.siamcat[i,i]=models_top[[i]]$aucinfo#直接调用
        }else{
          #i!=j
          if(i%in%add.id | j%in%add.id){
            #i行是训练集
            #j列是验证集
            
            cross_models[[paste(task_id[i],task_id[j],sep='-')]]=my_external.vali_siamcat(models_top[[i]],
                                                                                          test_feat = feat[[j]],
                                                                                          test_meta = meta[[j]])
            cross.siamcat[i,j]=cross_models[[paste(task_id[i],task_id[j],sep='-')]]$aucinfo
          }else{
            #不涉及新增数据集的数据（auc直接提取mysql->auc.matrix）
            old.i=which(rownames(auc.matrix)==task_id[i])
            old.j=which(colnames(auc.matrix)==task_id[j])
            #cross_models[[paste(task_id[i],task_id[j],sep='-')]]=
            cross.siamcat[i,j]=auc.matrix[old.i,old.j]
          }
        }
      }
    }
    
    
  }else{
    #new为空，则是训练全部feat的auc
    #1.模型对象
    if(is.null(models)){
      #训练模型
      models=list()
      for (i in 1:length(feat)) {
        models[[i]]=my_siamcat(train_feat = feat[[i]],train_meta =meta[[i]],method,num.folds,num.resample)
      }
    }
    models_top=models#与下面代码统一处理
    
    if(!is.null(top)){
      if(is.null(models)){
        stop('you must input models this time!')
        #当没有输入models时，会自动训练
      }
      models_top=lapply(models, function(data){
        data_top=my_siamcat.sel(siamcat=data$model_back,
                                consens.thres=0.5,
                                max.show=top,
                                verbose = 1,method=method)
      })
    }else{
      message("you don't choose top features!")
    }
    
    #2.组成一个auc矩阵
    task_id=names(feat)
    cross.siamcat <- matrix(NA,length(task_id),length(task_id))#储存AUC矩阵
    cross_models=list()
    for (i in 1:length(task_id)) {
      for (j in 1:length(task_id)) {
        if(i==j){
          cross.siamcat[i,i]=models_top[[i]]$aucinfo
        }else{
          cross_models[[paste(task_id[i],task_id[j],sep='-')]]=my_external.vali_siamcat(models_top[[i]],
                                                                                        test_feat = feat[[j]],
                                                                                        test_meta = meta[[j]])
          cross.siamcat[i,j]=cross_models[[paste(task_id[i],task_id[j],sep='-')]]$aucinfo
        }
      }
    }
  }
  
  #my_visual.R+final_models.R
  #3.数据准备
  tb = round(cross.siamcat, 2)#保留两位小数，为了图像美观
  tb=as.data.frame(tb)
  colnames(tb) = task_id
  rownames(tb) = task_id
  tb=cbind(TR=rownames(tb),method=label,tb)#行是训练集，列是验证集
  #行为训练集，列为验证集
  names(models)=names(feat)
  names(models_top)=names(feat)
  list(result=tb,models=models,models_top=models_top,cross_models=cross_models)
}
```

##### 2.6.6 模型运行

```R
setwd("F:\\Zhaolab2020\\gut-brain-axis\\metaAD\\04humann")
inputdir="F:\\Zhaolab2020\\gut-brain-axis\\metaAD\\04humann"
workdir="F:\\Zhaolab2020\\gut-brain-axis\\metaAD\\scripts"

pathways = t(read.table(file.path(inputdir, "metaAD_pathabundance-relab-filter.tsv"), sep="\t", header=TRUE, row.names = 1, stringsAsFactors = FALSE))
gene_families = t(read.table(file.path(inputdir, "metaAD_genefamilies-relab-filter.unstratified.tsv"), sep="\t", header=TRUE, row.names = 1, stringsAsFactors = FALSE))
species = t(read.table(file.path(workdir,"00profile", "merged_abundance_table_species.txt"), sep="\t", header=TRUE, row.names = 1, stringsAsFactors = FALSE))
species = species/100.0
print(species[1:5, 1:5])
auc_stat = tibble()
metadata = read.csv(file.path(workdir, "00metadata", "metadata_merge_filter.csv"), header=TRUE, row.names = 1, stringsAsFactors = FALSE)
for(feat_type in c('species', 'gene_families', 'pathways', 'species+gene_families', 'species+pathways', 'gene_families+pathways', 'species+gene_famlies+pathways')){
  #feat_type = 'species+gene_famlies+pathways'
  if(feat_type=='species'){
    feats = species
  }else if(feat_type=='gene_families'){
    feats = gene_families
  }else if(feat_type=='pathways'){
    feats = pathways
  }else if(feat_type=='species+gene_families'){
    feats = merge(species, gene_families, by.x='row.names', by.y="row.names")
    row.names(feats) = feats[,1]
    feats = feats[,-1]
  }else if(feat_type=='species+pathways'){
    feats = merge(species, pathways, by.x='row.names', by.y="row.names")
    row.names(feats) = feats[,1]
    feats = feats[,-1]
  }else if(feat_type=='gene_families+pathways'){
    feats = merge(gene_families, pathways, by.x='row.names', by.y="row.names")
    row.names(feats) = feats[,1]
    feats = feats[,-1]
  }else if(feat_type=='species+gene_famlies+pathways'){
    feats = merge(species, gene_families, by.x='row.names', by.y="row.names")
    row.names(feats) = feats[,1]
    feats = feats[,-1]
    feats = merge(feats, pathways, by.x='row.names', by.y="row.names")
    row.names(feats) = feats[,1]
    feats = feats[,-1]
  }
  print(feats[1:5, 1:5])
  print(metadata[1:5, ])
  #feats = t(feats)
  meta_feats = merge(metadata, feats, by.x='row.names', by.y="row.names")
  row.names(meta_feats) = meta_feats[,1]
  meta_feats = meta_feats[,-1]
  print(meta_feats[1:5, 1:10])
  write.table(t(meta_feats), file=paste0('meta_', feat_type,'.txt'), row.names = TRUE, col.names = TRUE)
  allfeat = list()
  topfeat = list()
  for(stage in c('SCS', 'SCD', 'MCI', 'AD')){
    #stage='AD'
    print(paste0("begin to deal with ",feat_type, '-',  stage))
    meta_feats_filter = meta_feats[(meta_feats$group=='NC') | (meta_feats$group==stage), ]
    meta_feats_filter$group[meta_feats_filter$group==stage] <- 'Case'
    colnames(meta_feats_filter)[4] = 'Group'
    print(meta_feats_filter[1:5, 1:5])
    meta_df = meta_feats_filter[, 1:4]
    feat_df = meta_feats_filter[, 5:length(meta_feats)]
    datas_AD = list()
    meta_list = list()
    datas_AD[[stage]] = feat_df
    meta_list[[stage]] = meta_df
    results_AD=my_result(datas_AD,meta_list,method='lasso',label = 'lasso_all')
    results_AD_top=my_result(datas_AD,meta_list,method = 'lasso',label =
                               'lasso_top20',models = results_AD$models,top = 20)
    # save(results_AD, file = file.path(workdir, 'train', paste0(model, "-all_models.RData")))
    # save(results_AD_top, file = file.path(workdir, ))
    allfeat[[stage]] <- results_AD
    topfeat[[stage]] <- results_AD_top
    
  }
  save(allfeat, file = file.path(workdir, 'train', paste0(feat_type, "all_models.RData")))
  save(topfeat, file = file.path(workdir, 'train', paste0(feat_type, "top_features_models.RData")))
  
  load(file.path(workdir, 'train', paste0(feat_type, "all_models.RData")))
  load(file.path(workdir, 'train', paste0(feat_type, "top_features_models.RData")))
  #auc_stat = tibble()
  allfeat$AD$result$AD[1]
  for(stage in c('SCS', 'SCD', 'MCI', 'AD')){
      for(model in c('lasso')){
          for(feat_num in c('all', 'top20')){
              
              if(feat_num == 'all'){
              auroc = allfeat[[stage]]$result[[stage]][1]
              tibble(stage=stage,
                     model=model, 
                     feature=feat,
                     AUC=c(auroc))
              auc_stat <- bind_rows(auc_stat, 
                                    tibble(stage=c(stage),
                                           model=c(model),
                                           feat_type=c(feat_type),
                                           feat_num=c(feat_num),
                                           AUC=c(auroc)))
              }else if(feat_num == 'top20'){
              auroc = topfeat[[stage]]$result[[stage]][1]
              auc_stat <- bind_rows(auc_stat, 
                    tibble(stage=c(stage),
                           model=c(model), 
                           feat_type=c(feat_type),
                           feat_num=c(feat_num),
                           AUC=c(auroc)))
              }
          }
       }
    }
  print(auc_stat) 
  
  ####################模型的可解释性。#################################
  ####模型的解释性。
  allfeat$AD$models[[1]]$model_back
  topfeat$AD$models_top[[1]]$model_back
  outdir = file.path(workdir, feat_type)
  if(! file.exists(outdir)){
    dir.create(outdir)
  }
  for(stage in c('SCS', 'SCD', 'MCI', 'AD')){
    for(model in c('lasso')){
      for(feat in c('all', 'top20')){
        if(feat == 'all'){
          all_siamcat = allfeat[[stage]]$models[[1]]$model_back
          model.interpretation.plot(
            all_siamcat,
            fn.plot = file.path(workdir, feat_type, paste(stage, model, feat_num, 'interpretation.pdf', sep='-')),
            consens.thres = 0.5,
            limits = c(-3, 3),
            heatmap.type = 'zscore',
          )
        }else if(feat == 'top20'){
          top_siamcat = topfeat[[stage]]$models_top[[1]]$model_back
          model.interpretation.plot(
            top_siamcat,
            fn.plot = file.path(workdir, feat_type, paste(stage, model, feat_num, 'interpretation.pdf', sep='-')),
            consens.thres = 0.5,
            limits = c(-3, 3),
            heatmap.type = 'zscore',
          )
        }
      }
    }
  }
}

write.csv(auc_stat, file='auc_stat.csv',  row.names = TRUE)
#rm(list = ls())
workdir="F:\\Zhaolab2020\\gut-brain-axis\\metaAD\\scripts"
auc_stat = read.csv(file.path(workdir, 'auc_stat.csv'), header=1, row.names =1)
print(auc_stat)
setwd("F:\\Zhaolab2020\\gut-brain-axis\\metaAD\\04humann")
feat_type <- auc_stat$feat_type
print(unique(feat_type))
auc_stat$feat_type <- factor(feat_type, levels = c('species', 'gene_families', 'pathways', 'species+gene_families', 'species+pathways', 'gene_families+pathways', 'species+gene_famlies+pathways'))
print(auc_stat)
#pdf("all_feat_and_models_auc.pdf", width=18, height=12)
png("all_feat_and_models_auc.png", width=3000, height=1800, res=300)
p1 <- ggplot(auc_stat, aes(x=feat_type, y=AUC)) + 
geom_bar(aes(fill=feat_num), stat="identity",position = "dodge") + 
geom_text(aes(label=AUC, group = feat_num), size = 3, position = position_dodge(width = 1)) +
labs(x = 'Group', y = 'AUC', title=paste0('predict by lasso')) +
facet_wrap( ~ stage, ncol=2) + 
scale_fill_manual(values = c("#8DD3C7", "#FB8072")) +
scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1)) +
theme(
    axis.text.x = element_text(size=10, hjust = 1, angle=30,  family="serif"),
    axis.text.y = element_text(size=10, hjust = 0.5, family="serif"),
    axis.title.x = element_text(size=12, hjust = 0.5, family="serif"),
    axis.title.y = element_text(size=12, hjust = 0.5, family="serif"),
    plot.title = element_text(size = 15, hjust=0.5, family="serif",),
    axis.line = element_line(colour = "black")
)
#print(p2)
print(p1)
dev.off() 
```

在服务器上的运行：

**Bug1:** 可以使用R4的环境，安装siamcat时可能出现错误`ERROR: configuration failed for package ‘stringi’`

解决办法参考：https://stackoverflow.com/questions/31942322/how-to-install-stringi-from-local-file-absolutely-no-internet-access [**行不通**]

```shell
conda install -c conda-forge r-stringi
```

然后继续安装其他包：

```R
options(repos=structure(c(CRAN="https://mirrors.sjtug.sjtu.edu.cn/cran/")))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SIAMCAT")
install.packages("randomForest")
install.packages(c("dplyr", "tidyr", "reshape2", "ggplot2", "cowplot", "RColorBrewer", "optparse"))
```



##### 2.6.7 解析models.Rdata文件

```R
library(SIAMCAT)

workdir = "F:\\Zhaolab2020\\gut-brain-axis\\metaAD\\scripts\\06predict"
feat_types = c('species', 'gene_families', 'pathways', 'species+gene_families', 'species+pathways', 'gene_families+pathways', 'species+gene_famlies+pathways')
feat_type = 'species'

model_data = load(file.path(workdir, "train", "speciesall_models.RData"))
print(allfeat)
print(allfeat$SCS)
print(allfeat$SCS$result)  ##输出为auroc。
print(allfeat$SCS$models)
print(allfeat$SCS$models$SCS)
print(allfeat$SCS$models$SCS$model_back)
print(label(allfeat$SCS$models$SCS$model_back))  ##返回的是label。

##filt_feat()相关操作
filt_feat = filt_feat(allfeat$SCS$models$SCS$model_back)
print(filt_feat)
print(filt_feat$filt.feat)
print(dim(filt_feat$filt.feat))  ##维数为(322,101),包括101个样本，322个特征。

##norm_feat()相关操作
norm_feat <- norm_feat(allfeat$SCS$models$SCS$model_back)
print(norm_feat)
print(norm_feat$norm.feat)  ##标准化以后的特征[322*101].
print(norm_feat$norm.param$retained.feat) ##保留的特征。
print(norm_feat$norm.param$feat.mean) ##特征均值。
print(norm_feat$norm.param$feat.adj.sd)  ##调整后的标准差（standard deviation）。

##data_split()相关操作。
data_split <- data_split(allfeat$SCS$models$SCS$model_back)
print(data_split)
print(data_split$training.folds[[1]][[1]])
print(data_split$test.folds[[1]][[1]])
print(data_split$num.resample)
print(data_split$num.folds)

##model_list()相关操作
model_list = model_list(allfeat$SCS$models$SCS$model_back)
print(model_list)
print(model_list$models[[15]])
print(model_list$model.type)
print(model_list$feature.type)

##feature_weights()相关操作
feature_weights = feature_weights(allfeat$SCS$models$SCS$model_back)
print(feature_weights) ##返回矩阵 【322*7】 平均权重、中位数权重、权重标准差、(相对)、权重百分数。

##pred_matrix()相关操作
pred_matrix = pred_matrix(allfeat$SCS$models$SCS$model_back) 
print(pred_matrix)  ##预测矩阵，包括3次重复中，每个样本的预测概率。

##eval_data()相关操作
eval_data = eval_data(allfeat$SCS$models$SCS$model_back) 
print(eval_data$auroc)
print(eval_data$auroc.all)
print(eval_data$prc$recall)
print(eval_data$prc$precision)
print(eval_data$auprc)
print(eval_data$auprc.all)
```

##### 2.6.8 绘制auroc曲线和PR_curve曲线

```R
model.evaluation.plot(SCS_NC=allfeat$SCS$models$SCS$model_back, 
                      SCD_NC=allfeat$SCD$models$SCD$model_back,
                      MCI_NC=allfeat$MCI$models$MCI$model_back,
                      AD_NC=allfeat$AD$models$AD$model_back,
                      fn.plot = "all_species.model.evaluatation.plot.pdf",
                      colours = c('blue', 'green', 'orange', 'red'))
```

结果如下图所示：

![image-20220309155234564](https://jialh.oss-cn-shanghai.aliyuncs.com/img/image-20220309155234564.png)

##### 2.6.9 模型解释

model.interpretation.plot: 模型可解释性绘图结果如下。

-   **第53行**

    ```R
    color.scheme <- rev(colorRampPalette(brewer.pal(brewer.pal.info[color.scheme, "maxcolors"], color.scheme))(100))
    ```

-   **第57行：model.interpretation.select.features**

    ```R
    sel.idx <- model.interpretation.select.features(feature.weights = feature.weights, 
            model.type = model.type, consens.thres = consens.thres, 
            label = label, max.show = max.show, verbose = verbose)
        num.sel.f <- length(sel.idx)
        mean.agg.pred <- rowMeans(pred_matrix(siamcat))
        srt.idx <- sort(label$label + mean.agg.pred, index.return = TRUE)$ix
        
    ##model.interpretation.select.features()中的内容。
    ```

-   【**<font color='red'>model.interpretation.plot函数相关解释，有待进一步补充。</font>**】



运行`model.interpretation.plot`函数的改进版本（`F:\Zhaolab2020\gut-brain-axis\metaAD\scripts\06check_model_data.R`）：

```R
source("F:\\Zhaolab2020\\gut-brain-axis\\metaAD\\scripts\\06predict\\01model.interpretation.plot.R")
cases <- c('SCS', 'SCD', 'MCI', 'AD')
sel.feats <- list()
for(case in cases){
    sc.obj <- allfeat[[case]]$models[[case]]$model_back
    model_interpretation <- model.interpretation.plot.modify(
        sc.obj,
        fn.plot = paste0('all_species_interpretation_', case, '_zscore.modify.pdf'),
        consens.thres = 0.5,
        limits = c(-3, 3),
        heatmap.type = 'zscore'
    )
    sel.feats[[case]] <- model_interpretation[["sel.feat"]]
    
}

print(sel.feats)

###绘制韦恩图
install.packages("VennDiagram")
library(VennDiagram)
#?venn.diagram
venn.diagram(list(SCS_NC=sel.feats$SCS,SCD_NC=sel.feats$SCD,MCI_NC=sel.feats$MCI,AD_NC=sel.feats$AD),
             filename="all_species.model.evaluatation.sel.feats.pdf",
             lwd=1,lty=2,
             fill=c('blue', 'green', 'orange', 'red'),
             cat.col=c('blue', 'green', 'orange', 'red'),reverse=TRUE)
```

结果为：

<img src="https://jialh.oss-cn-shanghai.aliyuncs.com/img/all_species.model.evaluatation.sel.feats_1.png" alt="all_species.model.evaluatation.sel.feats_1" style="zoom: 15%;" />

**相关特征在进化分支上，是否存在一致性？
