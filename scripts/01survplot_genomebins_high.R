library(survminer) # 加载包
library(survival) # 加载包
library(broom)

rm(list = ls())

infile = file.path("F:\\Zhaolab2020\\virome\\TGScompare\\02genome_bins", "high_genomes_modify.txt")
genomes = read.table(infile, sep='\t', header = TRUE, stringsAsFactors = FALSE)
print(unique(genomes$quality))
print(genomes$Sample)
high_genomes <- genomes[(genomes$quality=='high-quality') & (grep('A_', genomes['Sample'])), ]
levels(high_genomes$quality) <- c(levels(high_genomes$quality), 1) 
high_genomes$quality[high_genomes$quality=='high-quality']  <- 1 
high_genomes$quality <- as.numeric(as.character(high_genomes$quality))
dim(high_genomes)
print(high_genomes$Sample)
#View(high_genomes)

group <- high_genomes$group
print(unique(group))
high_genomes$group <- factor(group, levels = c("mNGS","PacBio", "Nanopore","mNGS+PacBio","mNGS+Nanopore"))

attach(high_genomes) # 绑定数据集
print(N50)
print(quality)
print(unique(Sample))
Surv(N50, quality) # 创建生存对象

uniq_group = unique(group)
print(uniq_group)

fit <- survfit(Surv(N50, quality) ~ group,  # 创建生存对象 
               data = high_genomes) # 数据集来源
fit # 查看拟合曲线信息

summary(fit, censored=TRUE)

ggsurvplot(fit, data = high_genomes)

ggsurvplot(fit, data = high_genomes, surv.scale='default', risk.table = TRUE,
           title = 'high-quality genomes',
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题
           legend.labs = c("mNGS", "PacBio","Nanopore","mNGS+PacBio", "mNGS+Nanopore"), # 指定图例分组标签
           surv.median.line = "hv", xlab='N50', ylab='bin ratio') # 增加中位生存时间
