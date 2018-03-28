options(digits=16)


lag_idiosyncracy_matrix=y_cycle[,-1]-observe_weights%*%final_factor[,-1]
advanced_idiosyncracy_matrix=y_cycle[,-ncol(y_cycle)]-observe_weights%*%final_factor[,-ncol(y_cycle)]

#特质波动自回归系数
idiosyncracy_autoreg_weights_miu=matrix(0,nrow=nrow(y_cycle),ncol=1)
idiosyncracy_autoreg_weights_sigma2=matrix(1,nrow=nrow(y_cycle),ncol=1)



i=1
for(i in 1:nrow(y_cycle))
{
  idiosyncracy_autoreg_weights_miu[i,]=solve(1+(t(lag_idiosyncracy_matrix[i,])%*%lag_idiosyncracy_matrix[i,])/noise_sigma[i,]^2)%*%((t(lag_idiosyncracy_matrix[i,])%*%advanced_idiosyncracy_matrix[i,])/noise_sigma[i,]^2)
  idiosyncracy_autoreg_weights_sigma2[i,]=solve(1+(t(lag_idiosyncracy_matrix[i,])%*%lag_idiosyncracy_matrix[i,])/noise_sigma[i,]^2)
  
  idiosyncracy_autoreg_weights[i,]=rnorm(1,mean=idiosyncracy_autoreg_weights_miu[i,],sd=sqrt(idiosyncracy_autoreg_weights_sigma2[i,]))
}


#噪声方差

for(i in 1:nrow(y_cycle))
{
  noise_sigma[i,]=rinvgamma(1,shape = (2+ncol(y_cycle))/2, scale = 1+0.5*(advanced_idiosyncracy_matrix[i,]-lag_idiosyncracy_matrix[i,]*idiosyncracy_autoreg_weights[i,])%*%t(advanced_idiosyncracy_matrix[i,]- lag_idiosyncracy_matrix[i,]*idiosyncracy_autoreg_weights[i,]))
}
