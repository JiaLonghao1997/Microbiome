#如果没有安装就先安装。
#install.packages("survminer") # 安装survminer包
#install.packages("survival") # 安装survival包
library(survminer) # 加载包
library(survival) # 加载包

data(lung) # 加载lung数据集
#View(lung) # 查看数据集

attach(lung) # 绑定数据集
print(unique(time))
print(unique(status))
print(unique(sex))
Surv(time,status) # 创建生存对象

fit <- survfit(Surv(time,status) ~ sex,  # 创建生存对象 
               data = lung) # 数据集来源
fit # 查看拟合曲线信息

summary(fit)

ggsurvplot(fit, data = lung)

ggsurvplot(fit, data = lung,
           surv.median.line = "hv") # 增加中位生存时间
