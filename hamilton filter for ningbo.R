
options(digits=16)

#通过Hamilton滤波更新共同因子所处状态


#Hamilton滤波
##变量声明

###共同因子修正
factor_for_state_estimation=matrix(0,nrow=1,ncol=(ncol(y_cycle)-1))

for(t in 1:(ncol(y_cycle)-1))
{
  factor_for_state_estimation[,t]=final_factor[,t+1]-factor_autoreg_weight*final_factor[,t]
}




###在上一期共同因子条件下上一期状态的概率
p_s=matrix(0,nrow=2,ncol=1)

p_s[1,]=dnorm(final_factor[,1],mean=factor_miu[1,],sd=1)/(dnorm(final_factor[,1],mean=factor_miu[1,],sd=1)+dnorm(final_factor[,1],mean=factor_miu[2,],sd=1))
p_s[2,]=dnorm(final_factor[,1],mean=factor_miu[2,],sd=1)/(dnorm(final_factor[,1],mean=factor_miu[1,],sd=1)+dnorm(final_factor[,1],mean=factor_miu[2,],sd=1))




###上一期与当期状态的联合概率（先列后行）,包括预测与更新的
p_joint_forecast=matrix(0,nrow=2,ncol=2)
p_joint_update=matrix(0,nrow=2,ncol=2)

###共同因子状态概率储存变量
p_s_forecast=matrix(0,nrow=2,ncol=ncol(y_cycle))
p_s_update=matrix(0,nrow=2,ncol=ncol(y_cycle))



###上一期到当期指标的概率
f_y_y=0

###在当期预测状态概率下取到当期指标观测值的概率
f_y_s=matrix(0,nrow=2,ncol=2)


##hamilton filter，所有时间状态的估计
t=1
for(t in 1:(ncol(y_cycle)))
{
  ###储存状态概率预测值与更新值
  p_s_forecast[,t]=rowSums(p_joint_forecast)
  p_s_update[,t]=p_s
  
  if(t==ncol(y_cycle))
  {
    break
  }
  
  ###预测阶段
  p_joint_forecast[,1]=p_miu[,1]*p_s[1,]
  p_joint_forecast[,2]=p_miu[,2]*p_s[2,]
    
  ###更新阶段
  f_y_s[1,1]=dnorm(factor_for_state_estimation[,t],mean=(1-factor_autoreg_weight)*factor_miu[1,],sd=1)
  f_y_s[2,1]=dnorm(factor_for_state_estimation[,t],mean=factor_miu[2,]-factor_autoreg_weight*factor_miu[1,],sd=1)
  f_y_s[1,2]=dnorm(factor_for_state_estimation[,t],mean=factor_miu[1,]-factor_autoreg_weight*factor_miu[2,],sd=1)
  f_y_s[2,2]=dnorm(factor_for_state_estimation[,t],mean=(1-factor_autoreg_weight)*factor_miu[2,],sd=1)
  
  f_y_y=p_joint_forecast[1,1]*f_y_s[1,1]+p_joint_forecast[1,2]*f_y_s[1,2]+p_joint_forecast[2,1]*f_y_s[2,1]+p_joint_forecast[2,2]*f_y_s[2,2]
  
  p_joint_update[1,1]=f_y_s[1,1]*p_joint_forecast[1,1]/f_y_y
  p_joint_update[2,1]=f_y_s[2,1]*p_joint_forecast[2,1]/f_y_y
  p_joint_update[1,2]=f_y_s[1,2]*p_joint_forecast[1,2]/f_y_y
  p_joint_update[2,2]=f_y_s[2,2]*p_joint_forecast[2,2]/f_y_y
  
  p_s[1,]=p_joint_update[1,1]+p_joint_update[1,2]
  p_s[2,]=p_joint_update[2,1]+p_joint_update[2,2]
  

  
}



###生成最后一期状态变量
state_variable[,(ncol(y_cycle)-1)]=sample(x=c(0,1),size = 1,prob = p_s_update[,ncol(y_cycle)])

##反向更新生成状态变量

###最终状态变量概率
p_s_final=matrix(0,nrow=2,ncol=(ncol(y_cycle)-1))
t=ncol(y_cycle)-1

for(t in (ncol(y_cycle)-1):1)
{
  p_s_final[1,t]=p_miu[state_variable[,t+1]+1,1]*p_s_update[1,t]/p_s_forecast[state_variable[,t+1]+1,t+1]
  p_s_final[2,t]=p_miu[state_variable[,t+1]+1,2]*p_s_update[2,t]/p_s_forecast[state_variable[,t+1]+1,t+1]

  
  
  state_variable[,t]=sample(x=c(0,1),size = 1,prob = p_s_final[,t])

}

