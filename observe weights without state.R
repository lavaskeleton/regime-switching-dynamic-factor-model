options(digits=16)

#观测方程系数更新

##共同因子修正
factor_for_observe_weights=matrix(0,nrow=nrow(y_cycle),ncol=(ncol(y_cycle)-1))



for(i in 1:nrow(y_cycle))
{
  for(t in 1:(ncol(y_cycle)-1))
  {
    factor_for_observe_weights[i,t]=final_factor[,t+1]-idiosyncracy_autoreg_weights[i,]*final_factor[,t]
  }  
}



observe_weights_miu=matrix(0,nrow=nrow(y_cycle),ncol=1)
observe_weights_sigma2=matrix(1,nrow=nrow(y_cycle),ncol=1)





i=1
for(i in 1:nrow(y_cycle))
{
  observe_weights_miu[i,]=solve(1+(t(factor_for_observe_weights[i,])%*%factor_for_observe_weights[i,])/noise_sigma[i,]^2)%*%((t(factor_for_observe_weights[i,])%*%y_cycle_for_factor_estimation[i,])/noise_sigma[i,]^2)
  observe_weights_sigma2[i,]=solve(1+(t(factor_for_observe_weights[i,])%*%factor_for_observe_weights[i,])/noise_sigma[i,]^2)
  
  observe_weights[i,]=rnorm(1,mean=observe_weights_miu[i,],sd=sqrt(observe_weights_sigma2[i,]))
}

