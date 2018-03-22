

source('read file and initialization for ningbo.R',encoding = "UTF-8")



#12000次循环gibbs抽样
gibbs_counter=1


#记录矩阵
p_miu_saver=matrix(0,nrow=2,ncol=2*50000) 

factor_autoreg_weight_saver=matrix(0,nrow=1,ncol=50000)

factor_miu_saver=matrix(0,nrow=2,ncol=50000)

observe_weights_saver=matrix(0,nrow=10,ncol=50000)

idiosyncracy_autoreg_weights_saver=matrix(0,nrow=10,ncol=50000)

noise_sigma_saver=matrix(0,nrow=10,ncol=50000)

for(gibbs_counter in 1:50000)
{
  source('kalman filter for ningbo.R',encoding = "UTF-8")
  
  source('hamilton filter for ningbo.R',encoding = "UTF-8")
  
  source('observe weights for ningbo.R',encoding = "UTF-8")
  
  source('idiosyncracy factor weights and noise sigma for ningbo.R',encoding = "UTF-8")
  
  source('common factor weights and miu for ningbo.R',encoding = "UTF-8")
 
  source('transition probability matrix for ningbo.R',encoding = "UTF-8")
  
  p_miu_saver[,((gibbs_counter-1)*2+1):(gibbs_counter*2)]=p_miu
  
  factor_autoreg_weight_saver[,gibbs_counter]=factor_autoreg_weight
  
  factor_miu_saver[,gibbs_counter]=factor_miu
  
  observe_weights_saver[,gibbs_counter]=observe_weights
  
  idiosyncracy_autoreg_weights_saver[,gibbs_counter]=idiosyncracy_autoreg_weights
  
  noise_sigma_saver[,gibbs_counter]=noise_sigma
  
}





plot(x=date,y=state_area1_average,type="l",col="red")
lines(x=date,y=state_area2_average,col="blue")
lines(x=date,y=state_area3_average,col="green")

matplot(y=t(y_cycle),type="l", xaxt='n')
axis(1,at=1:39,labels = date)





