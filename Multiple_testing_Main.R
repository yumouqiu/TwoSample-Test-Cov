############################################################
##
##     The code is for the multipule testing procedure 
##
############################################################

rm(list = ls())

n1s=c(150,150,150)     # sample size for the first population 
n2s=c(150,150,150)     # sample size for the second population
ps=c(100,200,400)      # dimension of covariance matrix

betas = c(0.5,0.6)     # sparsity level
r_1  = 0.8             # signals strength
r_2 = 0.6
r_3 = 0.5
prop_1 = 0.2           # proportion of different signal strength
prop_2 = 0.4
s0=0.5                 # thresholding lower bound
eta=0.05               # thresholding upper bound 1-eta
alpha=0.05             # significant level           

N=1000                  # simulation replication


writedir = getwd()
setwd("./codes_part") 

source("proposed.R")
source("maxthreshold.R")
source("intermediate.R")
source("DGP_Sigma_3.R")



for(beta in betas)
{
  for(iter in 1:length(n1s))
  {
    n1 = n1s[iter]
    n2 = n2s[iter]
    p = ps[iter]
    
    filename="power"
    library(PDSCE)
    library(doParallel)
    library(foreach)
    library(MASS)
    
    
    DGP_Sigma=DGP_Sigma_3(n1,n2,p,beta,r_1,r_2,r_3,prop_1,prop_2)   # Covariance type generation
    Sigma_1 = DGP_Sigma[[1]]
    Sigma_2 = DGP_Sigma[[2]]
    
    sam1=mvrnorm(n1,mu=rep(0,p),Sigma=Sigma_1)        # Data generation 
    sam2=mvrnorm(n2,mu=rep(0,p),Sigma=Sigma_2)        
    
    
    X=rbind(sam1,sam2)                                # Common covariance estimation for bootstrap calibration
    output=pdsoft.cv(x=X)
    Sigma=output$sigma
    
    
    
    parallel_num=50
    cl = makeCluster(parallel_num);
    registerDoParallel(cl);  
    
    result<-foreach(iter=1:N,.packages=c('foreach','parallel'),.combine=rbind) %dopar%
      {
        library(PDSCE)
        library(MASS)
        
        numProposed=0
        numboot=0
        q = p*(p+1)/2
        
        sam1=mvrnorm(n1,mu=rep(0,p),Sigma=Sigma_1)      # Data generation
        sam2=mvrnorm(n2,mu=rep(0,p),Sigma=Sigma_2)
        
        statboot=NULL           # Bootstrap resamples of the testing statistic                      
        B=200
        for(j in 1:B)
        {
          statboot[j]<-maxthreshold(Sigma,p,n1,n2,s_0,eta)
        }
        
        
        Sigma_diff = Sigma_1 - Sigma_2      # Testing statistic by using the original data
        medresults = intermediate(sam1,sam2)
        CovX = medresults[[1]]
        CovY = medresults[[2]]
        M = medresults[[3]]
        
        # The following are the multiple testing procedure based on MTT 
        
        i_j_pair = NULL
        
        for(i in 1:p)
        {
          for(j in i:p)
          {
            i_j_pair = rbind(i_j_pair,c(i,j))
          }
        }
        
        sorted_M = sort(M,decreasing = T,index.return=T)   
        i_j_pair_sorted = i_j_pair[sorted_M$ix,]       
        
        l_s = 1
        for(l in 1:q)
        {
          D_l = i_j_pair_sorted[l:q,]
          M_l = sorted_M$x[l:q]
          threshold_l=proposed(p,M_l,s0,eta)
          p_thr1<-length(which(sort(statboot)>=threshold_l[1]))/B
          if(p_thr1>=alpha)
          {
            l_s = l
            break
          }
        }
        
        c = 0.1
        
        if(l_s == 1)
        {
          Rej = 1
          FP = 0
        }else{
          l_s = l_s - 1
          l_a = min(q,floor((l_s-c)/(1-c)))
          R_a = i_j_pair_sorted[1:l_a,]
          Rej = dim(R_a)[1]
        }
        
        i_j_true = NULL
        for(i in 1:p)
        {
          for(j in i:p)
          {
            if(Sigma_diff[i,j]!=0)
            {
              i_j_true = rbind(i_j_true,c(i,j))
            }
          }
        }
        
        if(l_s>1)
        {
          TP = 0
          for(i in 1:dim(R_a)[1])
          {
            
            R_a_temp = i_j_true[i_j_true[,1] == R_a[i,1],] 
            if(!is.null(R_a_temp))
            {
              R_a_temp = matrix(R_a_temp,ncol=2)
              R_a_temp = R_a_temp[R_a_temp[,2] == R_a[i,2],]
              if(length(R_a_temp)>0)
              {
                TP = TP + 1
              }
            }
            
          }
          FDP = (Rej-TP)/max(Rej,1)
          power = TP/(dim(i_j_true)[1])
          indicator = 0 
          if(FDP > c)
          {
            indicator = 1
          }
        }
        if(l_s == 1)
        {
          FDP = 0 
          power = 0
          indicator = 0
        }
        
      
        ############### maximum procedure ##########################
        # The following are the multiple testing procedure based on the maximum type test
        
        M = medresults[[3]]
        sorted_M = sort(M,decreasing = T,index.return=T)
        i_j_pair_sorted = i_j_pair[sorted_M$ix,]
        
        l_s_max = 1
        for(l in 1:q)
        {
          M_l_max = sorted_M$x[l]
          if(M_l_max<=(-log(8*pi)-2*log(log(1/(1-alpha)))+4*log(p)-log(log(p))))
          {
            l_s_max = l
            break
          }
        }
        
        c = 0.1
        if(l_s_max == 1)
        {
          Rej_max = 1
          FP_max = 0
        }else{
          l_s_max = l_s_max - 1
          l_a_max = min(q,floor((l_s_max-c)/(1-c)))
          R_a_max = i_j_pair_sorted[1:l_a_max,]
          Rej_max = dim(R_a_max)[1]
        }
        
        i_j_true = NULL
        for(i in 1:p)
        {
          for(j in i:p)
          {
            if(Sigma_diff[i,j]!=0)
            {
              i_j_true = rbind(i_j_true,c(i,j))
            }
          }
        }
        
        if(l_s_max>1)
        {
          TP_max = 0
          for(i in 1:dim(R_a_max)[1])
          {
            
            R_a_temp_max = i_j_true[i_j_true[,1] == R_a_max[i,1],] 
            if(!is.null(R_a_temp_max))
            {
              R_a_temp_max = matrix(R_a_temp_max,ncol=2)
              R_a_temp_max = R_a_temp_max[R_a_temp_max[,2] == R_a_max[i,2],]
              if(length(R_a_temp_max)>0)
              {
                TP_max = TP_max + 1
              }
            }
            
          }
          FDP_max = (Rej_max-TP_max)/max(Rej_max,1)
          power_max = TP_max/(dim(i_j_true)[1]+10^(-6))
          indicator_max = 0 
          if(FDP_max > c)
          {
            indicator_max = 1
          }
        }
        if(l_s_max == 1)
        {
          FDP_max = 0 
          power_max = 0
          indicator_max = 0
        }
        
        c(FDP,power,indicator,FDP_max,power_max,indicator_max)
      }
    stopCluster(cl)
    FDR=mean(result[,1],na.rm = T)          # FDR    
    power=mean(result[,2],na.rm = T)        # Power
    prop=mean(result[,3],na.rm = T)         # FDPEx
    
    FDR_max=mean(result[,4],na.rm = T)
    power_max=mean(result[,5],na.rm = T) 
    prop_max=mean(result[,6],na.rm = T) 
   
    write.csv(data.frame(FDR=FDR,power=power,FDPEx=prop,FDR_max=FDR_max,power_max=power_max,FDPEx_max=prop_max),file=paste0(writedir,paste0("/","Summerized_",p,"_",n1,"_",beta,".csv")))
     
    
    
    
  }  
}


