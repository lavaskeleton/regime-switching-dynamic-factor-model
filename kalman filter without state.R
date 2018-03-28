options(digits=16)


#通过kalman滤波更新共同因子序列


#模型形式转换
##观察值的加权差分
y_cycle_for_factor_estimation=matrix(0,nrow=nrow(y_cycle),ncol=(ncol(y_cycle)-1))
i=1
for(i in 1:nrow(y_cycle))
{
  y_cycle_for_factor_estimation[i,]=y_cycle[i,-1]-idiosyncracy_autoreg_weights[i,]*y_cycle[i,-ncol(y_cycle)]
}


#量测方程系数
observe_weights_for_factor_estimation=matrix(0,nrow=nrow(y_cycle),ncol=2)

i=1
for(i in 1:nrow(y_cycle))
{
  observe_weights_for_factor_estimation[i,]=cbind(observe_weights[i,],-observe_weights[i,]*idiosyncracy_autoreg_weights[i,])
}

##量测方程噪声协方差矩阵
observe_disturbance_variance=diag(noise_sigma[1:nrow(y_cycle)]^2)



##转移方程截距
# observe_intercept=matrix(0,nrow=2,ncol=(ncol(y_cycle)-1))
# 
# t=1
# for(t in 1:(ncol(y_cycle)-1))
# {
#   observe_intercept[,t]=matrix(c(factor_miu[state_variable[,t+1]+1,]-factor_autoreg_weight*factor_miu[state_variable[,t]+1,],0),2,1)
# }

##转移方程系数
factor_weights_for_factor_estimation=cbind(rbind(factor_autoreg_weight,1),matrix(0,nrow = 2, ncol = 1))


##转移方程噪声协方差矩阵
factor_disturbance_variance=matrix(c(1,0,0,0),2,2)


##更新方程参数矩阵
update_parameter=matrix(0,nrow=nrow(y_cycle),ncol=nrow(y_cycle)-2)





##上一期共同因子
if(gibbs_counter==1)
{
  factor_for_former=rbind(principal_variable[,2],principal_variable[,1])
}else if(gibbs_counter>=1)
{
  factor_for_former=rbind(final_factor[,2],final_factor[,1])
}

##上一期共同因子估计误差协方差矩阵
if(gibbs_counter==1)
{
  factor_error_variance_for_former=matrix(c(1,0,0,1),nrow=2,ncol=2)

}else if(gibbs_counter>1)
{
  factor_error_variance_for_former=factor_error_variance_estimation[(1:2),(1:2)]
}




#所有时间的共同因子估计,kalman filter

##共同因子储存变量
factor_estimation=matrix(0,nrow=2,ncol=(ncol(y_cycle)-1))

factor_error_variance_estimation=matrix(0,nrow=2,ncol=2*(ncol(y_cycle)-1))




t=1

for(t in 1:(ncol(y_cycle)-1))
{    
  ##记录状态变量与状态估计误差
  factor_estimation[,t]=factor_for_former
  factor_error_variance_estimation[,((t-1)*2+1):(t*2)]=factor_error_variance_for_former
  
  if(t==(ncol(y_cycle)-1))
  {
    break
  }
  
  ##预测阶段
  
  #factor_for_former=observe_intercept[,t]+factor_weights_for_factor_estimation%*%factor_for_former
  factor_for_former=factor_weights_for_factor_estimation%*%factor_for_former
  
  factor_error_variance_for_former=factor_weights_for_factor_estimation%*%factor_error_variance_for_former%*%t(factor_weights_for_factor_estimation)+factor_disturbance_variance
  
  ##更新阶段
  update_parameter=observe_weights_for_factor_estimation%*%factor_error_variance_for_former%*%t(observe_weights_for_factor_estimation)+observe_disturbance_variance
  
  factor_for_former=factor_for_former+factor_error_variance_for_former%*%t(observe_weights_for_factor_estimation)%*%solve(update_parameter)%*%(y_cycle_for_factor_estimation[,t]-observe_weights_for_factor_estimation%*%factor_for_former)
  factor_error_variance_for_former=factor_error_variance_for_former-factor_error_variance_for_former%*%t(observe_weights_for_factor_estimation)%*%solve(update_parameter)%*%observe_weights_for_factor_estimation%*%factor_error_variance_for_former
}


##反向更新得到最终的共同因子估计

update_factor_estimation=matrix(0,nrow=2,ncol=(ncol(y_cycle)-1))
update_factor_error_variance_estimation=matrix(0,nrow=2,ncol=2*(ncol(y_cycle)-1))

####R
update_update_parameter=matrix(0,nrow=1,ncol=1)
####η
update_update_parameter1=matrix(0,nrow=1,ncol=1)

final_factor=matrix(0,nrow=1,ncol=ncol(y_cycle))

t=ncol(y_cycle)-1

update_factor_estimation[,t]=factor_estimation[,t]
final_factor[,t+1]=as.numeric(factor_estimation[1,t])


for(t in (ncol(y_cycle)-1):1)
{    
  ###更新阶段
  ####R
  update_update_parameter=factor_weights_for_factor_estimation[1,]%*%factor_error_variance_estimation[,((t-1)*2+1):(t*2)]%*%as.matrix(factor_weights_for_factor_estimation[1,])+factor_disturbance_variance[1,1]   
  ####η
  #update_update_parameter1=(final_factor[,t+1]-observe_intercept[1,t]-as.numeric(factor_weights_for_factor_estimation[1,]%*%factor_estimation[,t]))
  update_update_parameter1=(final_factor[,t+1]-as.numeric(factor_weights_for_factor_estimation[1,]%*%factor_estimation[,t]))
  
  update_factor_estimation[,t]=factor_estimation[,t]+factor_error_variance_estimation[,((t-1)*2+1):(t*2)]%*%factor_weights_for_factor_estimation[1,]*as.numeric(update_update_parameter1)/as.numeric(update_update_parameter)
  
  update_factor_error_variance_estimation[,((t-1)*2+1):(t*2)]=factor_error_variance_estimation[,((t-1)*2+1):(t*2)]-factor_error_variance_estimation[,((t-1)*2+1):(t*2)]%*%as.matrix(factor_weights_for_factor_estimation[1,])%*%factor_weights_for_factor_estimation[1,]%*%factor_error_variance_estimation[,((t-1)*2+1):(t*2)]/as.numeric(update_update_parameter)
  
  ###生成共同因子
  final_factor[,t]=rnorm(1,mean=update_factor_estimation[1,t],sd=sqrt(update_factor_error_variance_estimation[1,(t-1)*2+1]))
  
}






#plot(c_estimation,type = "l")

