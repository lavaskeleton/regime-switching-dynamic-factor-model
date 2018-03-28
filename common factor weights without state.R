options(digits=16)


#共同因子修正
factor_for_autoreg_weight=matrix(0,nrow=1,ncol=ncol(y_cycle))


for(t in 1:ncol(y_cycle))
{
  #factor_for_autoreg_weight[,t]=final_factor[,t]-factor_miu[state_variable[,t]+1,]
  factor_for_autoreg_weight[,t]=final_factor[,t]
}

lag_factor_for_autoreg_weight=factor_for_autoreg_weight[,-ncol(y_cycle)]
advanced_factor_for_autoreg_weight=factor_for_autoreg_weight[,-1]

#共同因子自回归系数更新
#均值
factor_autoreg_weight_miu=solve(1+t(lag_factor_for_autoreg_weight)%*%lag_factor_for_autoreg_weight)%*%(t(lag_factor_for_autoreg_weight)%*%advanced_factor_for_autoreg_weight)
#方差
factor_autoreg_weight_sigma2=solve(1+t(lag_factor_for_autoreg_weight)%*%lag_factor_for_autoreg_weight)


factor_autoreg_weight=rnorm(1,mean=factor_autoreg_weight_miu,sd=sqrt(factor_autoreg_weight_sigma2))



#两种状态下的均值更新
# 
# lag_state_matrix_for_miu_estimation=matrix(0,nrow=1,ncol=(ncol(y_cycle)-1))
# advanced_state_matrix_for_miu_estimation=matrix(0,nrow=1,ncol=(ncol(y_cycle)-1))
# 
# for(t in 1:(ncol(y_cycle)-1))
# {
#   lag_state_matrix_for_miu_estimation[,t]=state_variable[,t]
#   advanced_state_matrix_for_miu_estimation[,t]=state_variable[,t+1]
# }
# 
# 
# state_matrix_for_miu_estimation=rbind((1-advanced_state_matrix_for_miu_estimation)-factor_autoreg_weight*(1-lag_state_matrix_for_miu_estimation),advanced_state_matrix_for_miu_estimation-factor_autoreg_weight*lag_state_matrix_for_miu_estimation)
# 
# state_matrix_for_miu_estimation=t(state_matrix_for_miu_estimation)
# 
# priori_miu_sigma2=diag(2)
# 
# miu_for_miu_estimation=solve(priori_miu_sigma2+t(state_matrix_for_miu_estimation)%*%state_matrix_for_miu_estimation)%*%(t(state_matrix_for_miu_estimation)%*%t(factor_for_state_estimation))
# 
# sigma2_for_miu_estimation=solve(priori_miu_sigma2+t(state_matrix_for_miu_estimation)%*%state_matrix_for_miu_estimation)
# 
# 
# 
# factor_miu=as.matrix(mvrnorm(1,mu=miu_for_miu_estimation,Sigma = sigma2_for_miu_estimation))
