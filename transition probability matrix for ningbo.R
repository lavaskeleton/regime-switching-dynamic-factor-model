options(digits=16)

#转移概率矩阵更新

counter_0_0=0
counter_0_1=0
counter_1_0=0
counter_1_1=0


for(t in 1:(ncol(y_cycle)-1))
{
  if(state_variable[,t]==0&&state_variable[,t+1]==0)
  {  
    counter_0_0=counter_0_0+1
  }else if(state_variable[,t]==0&&state_variable[,t+1]==1)
  {
    counter_0_1=counter_0_1+1
  }else if(state_variable[,t]==1&&state_variable[,t+1]==0)
  {
    counter_1_0=counter_1_0+1
  }else if (state_variable[,t]==1&&state_variable[,t+1]==1)
  {
    counter_1_1=counter_1_1+1
  }
}


p_miu[1,1]=rbeta(1, shape1 = 8+counter_0_0, shape2 = 2+counter_0_1)

p_miu[2,2]=rbeta(1, shape1 = 8+counter_1_1, shape2 = 2+counter_1_0)

p_miu[2,1]=1-p_miu[1,1]

p_miu[1,2]=1-p_miu[2,2]

