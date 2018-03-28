library(openxlsx)
library(actuar)
library(psych)
library(MASS)
library(mnormt)
library(mFilter)
library(tseries)
library(urca)
library(forecast)



options(digits=16)
#环境与指标预处理
##设定工作空间

setwd("C:/Users/Administrator/Desktop/方老师宁波电力景气指数项目/hierarchical/markov switching hierarchical dynamic factors model")



##读取数据（xlsx文件中，xls文件不行）

#y_index=read.xlsx(xlsxFile = "宁波宏观电力景气指标.xlsx",sheet = 1,detectDates = TRUE, colNames = TRUE)
y_index=read.xlsx(xlsxFile = "宁波宏观电力景气指标.xlsx",sheet = 1,detectDates = TRUE, colNames = TRUE)

#gdp=read.xlsx(xlsxFile = "宁波宏观电力景气指标.xlsx",sheet = 5,detectDates = TRUE, colNames = TRUE)

#y_GDP=read.xlsx(xlsxFile = "宁波市GDP.xlsx",sheet = 2,detectDates = TRUE, colNames = TRUE)

#y_index=read.xlsx(xlsxFile = "全宁波数据V1.0.xlsx",sheet = 1,detectDates = TRUE, colNames = TRUE)



##正负指标调整（逆指标调整为负数）
# j=2
# for(j in 2:ncol(y_index))
# {
#   if(lm(y_index[,j]~y_index[,1])$coefficients[2]<0)
#   {
#     y_index[,j]=0-y_index[,j]
#   }
# }




##去除日期列
y_index=y_index[,-1]
#gdp=gdp[,-1]

#cov(gdp,y_index[,3])

# 
# 
#regression=glm(gdp~y_index[,1])
# 
#summary(regression)
# # 
# 
# ##标准化
# j=1
# for(j in 1:ncol(y_index))
# {
#   y_index[,j]=(y_index[,j]-min(y_index[,j]))/(max(y_index[,j])-min(y_index[,j]))*40+60
# }

##周期因子与随机因素项，后续计算使用的数据序列

y_cycle=matrix(0,nrow=nrow(y_index),ncol=ncol(y_index))

for(j in 1:ncol(y_index))
{
  y_index[,j]=(y_index[,j]-mean(y_index[,j]))/sd(y_index[,j])

}

##季节调整（ts调整为时间序列，stl进行分解，seasadj去除季节因素）

j=1
for(j in 1:ncol(y_index))
{
  #y_index[,j]=as.vector(seasadj(stl(ts(y_index[,j],frequency = 12),s.window = "periodic")))
  y_index[,j]=as.vector(seasadj(stl(ts(y_index[,j],frequency = 12),s.window="periodic")))

}
# 
# 
# 
# 
# # y_GDP[,4]=as.vector(seasadj(decompose(ts(y_GDP[,4],frequency = 4),type="multiplicative")))
# # y_index[,1]=as.vector(seasadj(stl(ts(y_index[,1],frequency = 12),s.window = "periodic")))
# # 
# # write.csv(y_index[,1],file = "全行业用电量去季节化.csv")
# 
# 
# 
##hp滤波（月数据为129600，季度为1600，年为6.25）
#y_cycle=y_index
j=1
for(j in 1:ncol(y_index))
{
  y_cycle[,j]=hpfilter(y_index[,j],type="lambda",freq=129600)$cycle
}




#主成分变量作为变量初始值

##主成分分析（每一列为变量，从变量中提取）

rotated_principal_object=principal(y_cycle,nfactors=1,rotate="oblimin")

#fa.parallel(y_cycle,fa="both",n.iter=100,main="Scree plots with parallel analysis")  


##获得主成分变量


weights_for_principal=t(as.matrix(rotated_principal_object$weights))
principal_variable=matrix(0,nrow=nrow(y_cycle),ncol=1)

i=1
j=1

for(i in 1:nrow(y_cycle))
{
  for(j in 1:ncol(y_cycle))
  {
    principal_variable[i,]=principal_variable[i,]+weights_for_principal[,j]*y_cycle[i,j]
  }
}
#write.table(x = principal_variable,file = "principal.csv",sep = ",")


#通过OLS和自回归得到模型系数初始值

##观测方程系数初始值，为各个指标对主成份的OLS回归系数
observe_weights=matrix(0,nrow=ncol(y_cycle),ncol=1)

# j=1
# for(j in 1:ncol(y_cycle))
# {
#   observe_weights[j,]=rnorm(1,mean=0,sd=1)
# }

reg_model=rep(list(0),ncol(y_cycle))


j=1
for(j in 1:ncol(y_cycle))
{
  reg_model[[j]]=lm(y_cycle[,j]~principal_variable)
  observe_weights[j,]=reg_model[[j]]$coefficients[2]
}



#状态方程系数初始值与特质波动系数初始值，为主成份和特质波动的自回归系数

##特质波动自回归系数初始值
idiosyncracy_autoreg_weights=matrix(0,nrow=ncol(y_cycle),ncol=1)
#idiosyncracy_autoreg_weights=mvrnorm(ncol(y_cycle),mu = c(rep(0,10)),Sigma = diag(c(rep(1,10))))

idiosyncracy_autoreg_model=rep(list(0),ncol(y_cycle))
idiosyncracy_matrix=matrix(0,nrow=nrow(y_cycle),ncol = ncol(y_cycle))

j=1
for(j in 1:ncol(y_cycle))
{
  idiosyncracy_matrix[,j]=reg_model[[j]]$residuals
}

j=1
for(j in 1:ncol(y_cycle))
{
  idiosyncracy_autoreg_model[[j]]=arma(idiosyncracy_matrix[,j],order = c(1,0))
  idiosyncracy_autoreg_weights[j,]=idiosyncracy_autoreg_model[[j]]$coef[1]
}

##主成份自回归系数，也就是状态方程系数初始值

factor_autoreg_model=arma(principal_variable,order = c(1,0))

factor_autoreg_weight=factor_autoreg_model$coef[1]



#factor_autoreg_weight=rnorm(1,mean = 0,sd=1)


#经济周期两种不同状态下均值的先验分布

#factor_miu=as.matrix(mvrnorm(n=1, mu = c(-4,4), Sigma = matrix(c(4^2,0,0,4^2),2,2)))


#特质波动自回归噪声方差的先验分布

noise_sigma=as.matrix((1/rgamma(ncol(y_cycle),shape=1,rate=1))^(1/2))

#行列调整（行为指标，列为时间）
y_cycle=t(y_cycle)
idiosyncracy_matrix=t(idiosyncracy_matrix)
principal_variable=t(principal_variable)

#markov链初始转移概率，分别为0到0,0到1,1到0和1到1的转移概率(先列后行)

# p_miu=matrix(0,2,2)
# 
# p_miu[1,1]=rbeta(1,1,1)
# p_miu[2,2]=rbeta(1,1,1)
# 
# p_miu[2,1]=1-p_miu[1,1]
# p_miu[1,2]=1-p_miu[2,2]
# 
# 
# #经济所属状态识别变量
# 
# state_variable=matrix(0,nrow=1,ncol=ncol(y_cycle))



