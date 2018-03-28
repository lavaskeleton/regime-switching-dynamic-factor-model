library(openxlsx)
library(actuar)
library(psych)
library(MASS)
library(mnormt)
library(mFilter)
library(tseries)
library(urca)
library(forecast)

setwd("C:/Users/Administrator/Desktop/方老师宁波电力景气指数项目/hierarchical/markov switching hierarchical dynamic factors model")


source('initialization without rigime switching.R',encoding = "UTF-8")



#12000次循环gibbs抽样
gibbs_counter=1


#记录矩阵
#p_miu_saver=matrix(0,nrow=2,ncol=2*35000) 

factor_autoreg_weight_saver=matrix(0,nrow=1,ncol=35000)

factor_saver=matrix(0,nrow=ncol(y_cycle),ncol=35000)

#factor_miu_saver=matrix(0,nrow=2,ncol=35000)

observe_weights_saver=matrix(0,nrow=ncol(y_index),ncol=35000)

idiosyncracy_autoreg_weights_saver=matrix(0,nrow=ncol(y_index),ncol=35000)

noise_sigma_saver=matrix(0,nrow=ncol(y_index),ncol=35000)

for(gibbs_counter in 1:35000)
{
  source('kalman filter without state.R',encoding = "UTF-8")
  
  #source('hamilton filter for ningbo.R',encoding = "UTF-8")
  
  source('observe weights without state.R',encoding = "UTF-8")
  
  source('idiosyncracy factor weights and noise sigma without state.R',encoding = "UTF-8")
  
  source('common factor weights without state.R',encoding = "UTF-8")
  
  #source('transition probability matrix for ningbo.R',encoding = "UTF-8")
  
  #p_miu_saver[,((gibbs_counter-1)*2+1):(gibbs_counter*2)]=p_miu
  
  factor_autoreg_weight_saver[,gibbs_counter]=factor_autoreg_weight
  
  #factor_miu_saver[,gibbs_counter]=factor_miu
  
  factor_saver[,gibbs_counter]=t(final_factor)
  
  observe_weights_saver[,gibbs_counter]=observe_weights
  
  idiosyncracy_autoreg_weights_saver[,gibbs_counter]=idiosyncracy_autoreg_weights
  
  noise_sigma_saver[,gibbs_counter]=noise_sigma
  
}

# p_miu_average=matrix(0,nrow=2,ncol=1)
# p_miu_sd=matrix(0,nrow=2,ncol=1)
# p_miu_median=matrix(0,nrow=2,ncol=1)

factor_autoreg_weight_average=matrix(0,nrow=1,ncol=1)
factor_autoreg_weight_sd=matrix(0,nrow=1,ncol=1)
factor_autoreg_weight_median=matrix(0,nrow=1,ncol=1)

# factor_miu_average=matrix(0,nrow=2,ncol=1)
# factor_miu_sd=matrix(0,nrow=2,ncol=1)
# factor_miu_median=matrix(0,nrow=2,ncol=1)

factor_average=matrix(0,nrow=1,ncol=ncol(y_cycle))

observe_weights_average=matrix(0,nrow=ncol(y_index),ncol=1)
observe_weights_sd=matrix(0,nrow=ncol(y_index),ncol=1)
observe_weights_median=matrix(0,nrow=ncol(y_index),ncol=1)

idiosyncracy_autoreg_weights_average=matrix(0,nrow=ncol(y_index),ncol=1)
idiosyncracy_autoreg_weights_sd=matrix(0,nrow=ncol(y_index),ncol=1)
idiosyncracy_autoreg_weights_median=matrix(0,nrow=ncol(y_index),ncol=1)

noise_sigma_average=matrix(0,nrow=ncol(y_index),ncol=1)
noise_sigma_sd=matrix(0,nrow=ncol(y_index),ncol=1)
noise_sigma_median=matrix(0,nrow=ncol(y_index),ncol=1)


# p_miu_average[1,]=mean(p_miu_saver[1,5001:35000])
# p_miu_average[2,]=mean(p_miu_saver[2,5001:35000])
# 
# p_miu_sd[1,]=sd(p_miu_saver[1,5001:35000])
# p_miu_sd[2,]=sd(p_miu_saver[2,5001:35000])
# 
# p_miu_median[1,]=median(p_miu_saver[1,5001:35000])
# p_miu_median[2,]=median(p_miu_saver[2,5001:35000])

# 
# factor_miu_average[1,]=mean(factor_miu_saver[1,5001:35000])
# factor_miu_average[2,]=mean(factor_miu_saver[2,5001:35000])
# 
# factor_miu_sd[1,]=sd(factor_miu_saver[1,5001:35000])
# factor_miu_sd[2,]=sd(factor_miu_saver[2,5001:35000])
# 
# factor_miu_median[1,]=median(factor_miu_saver[1,5001:35000])
# factor_miu_median[2,]=median(factor_miu_saver[2,5001:35000])

for(i in 1:nrow(y_index))
{
  factor_average[,i]=mean(factor_saver[i,5001:35000])
}
# 
# factor_added=matrix(0,nrow=1,ncol=ncol(y_cycle))
# factor_added[,1]=factor_average[,1]
# for(j in 2:ncol(y_cycle))
# {
#   factor_added[,j]=factor_added[,j-1]+factor_average[,j-1]
# }



# 
# dates=seq.Date(as.Date("2013-01-01"), by="1 month", length.out=nrow(y_index))
# dates2=seq.Date(as.Date("2013-01-01"), by="1 quarter", length.out=20)


# plot(x=dates[-(1:4)],y=factor_average[,-((ncol(y_cycle)-3):ncol(y_cycle))]*5,type="l",col="red")
# plot(x=dates,y=y_cycle[1,-(1:3)]/100,type="l",col="blue")
# lines(x=dates[-(1:4)],y=y_index[-(1:4),1]/1500,type="l",col="blue")
# lines(x=dates[-(1:4)],y=y_index[-(1:4),2]/1000,type="l",col="green")
# lines(x=dates[-(1:4)],y=y_index[-(1:4),3]/1000,type="l",col="yellow")





write.table(t(factor_average),file="factor.csv",sep=",")



#write.table(y_GDP,file="seasoned_GDP.csv",sep=",")

cor(y_cycle[,2],y_cycle[,3])
cor(factor_added[,-((ncol(y_cycle)-3):ncol(y_cycle))],y_index[-(1:4),3])
cor(factor_average[,-1],y_cycle[3,-ncol(y_cycle)])


