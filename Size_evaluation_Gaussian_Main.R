########################################################################################
##         Main part: used to calculate the size of the global MTT test               ##
##                                                                                    ##
#########################################################################################

rm(list = ls())

n1s=c(80,100,120)        # n1
n2s=c(80,100,120)        # n2
ps=c(277,396,530)         # p
 


writedir = getwd()
setwd("./codes_part") 
source("proposed.R")
source("maxthreshold.R")
source("intermediate.R")
source("DGP_data.R") 
source("DGP_Sigma_1.R")
source("DGP_Sigma_2.R")

for(i in 1:length(n1s))
{
  
n1 = n1s[i]
n2 = n2s[i]
p = ps[i]
  
beta=0.5                     # signal sparsity level, no use in size evaluation
r = 0                        # signal strength
s0=0.5                       # thresholding lower bound
eta=0.05                     # thresholding upper bound 1-eta
alpha=0.05                   # significant level
 
N=1000                       # simulation replication
filename="size"
dist = "Gaussian"            # Data type, for Gamma type, use dist = "Gamma"   
type = "DGP_1"               # Covariance type Design 1, for Design 2, use type = "DGP_2" 

library(PDSCE)
library(doParallel)
library(foreach)
 
# Covariance generation

if(type == "DGP_1")              
{
DGP_Sigma=DGP_Sigma_1(n1,n2,p,beta,r)
Sigma_1 = DGP_Sigma[[1]]
Sigma_2 = DGP_Sigma[[2]]
}
if(type == "DGP_2")
{
  DGP_Sigma=DGP_Sigma_2(n1,n2,p,beta,r)
  Sigma_1 = DGP_Sigma[[1]]
  Sigma_2 = DGP_Sigma[[2]]
}

# Data generation
Data1 = DGP_data(Sigma_1,Sigma_2,n1,n2,p,beta,r,dist)
X_0 = Data1[[1]]
Y_0 = Data1[[2]]

# Estimation the common covariance for the bootstrap calibration
permuation = sample(1:p)
sam1 = X_0[,permuation]
sam2 = Y_0[,permuation] 
X=rbind(sam1,sam2)
output=pdsoft.cv(x=X)
Sigma=output$sigma

# Parallel calculation
parallel_num=50

cl = makeCluster(parallel_num);
registerDoParallel(cl);  


result<-foreach(iter=1:N,.packages=c('foreach','parallel'),.combine=rbind) %dopar%
{
  library(PDSCE)
 
  numProposed=0     # Record rejection times for MTT by using the limit distribution
  numboot=0         # Record rejection times for MTT by using the bootstrap carlibration
  
  Data1=DGP_data(Sigma_1,Sigma_2,n1,n2,p,beta,r,dist)  # Generate the data
  X =Data1[[1]]
  Y =Data1[[2]]
  sam1 = X[,permuation]
  sam2 = Y[,permuation] 
  
  medresults=intermediate(sam1,sam2)     # Calculate the covariance and M_{ij}
  CovX=medresults[[1]]
  CovY=medresults[[2]]
  M=medresults[[3]]                      # vector M_{ij}
  
  threshold=proposed(p,M,s0,eta)         # The thresholding procedure
  
  if(threshold[[2]]<=alpha)              # MTT based on the limiting distribution
  {
    numProposed=numProposed+1
  }
  
  
   statboot=NULL                        # Bootstrap calibration
   B=100
   for(j in 1:B)
   {
     statboot[j]<-maxthreshold(Sigma,p,n1,n2,s_0,eta)
   }
   
   p_thr1<-length(which(sort(statboot)>=threshold[1]))/B
   if(p_thr1<=alpha)
   {
     numboot= numboot+1
   }
   
   c(numProposed,numboot)
}
 
Power_thre=sum(result[,1])/N
Power_boot=sum(result[,2])/N


write.csv(data.frame(MTT=Power_thre,MTT_BT=Power_boot),file=paste0(writedir,paste0("/",filename,"_",dist,"_",type,"_",p,"_",n1,"_",beta,".csv")))
stopCluster(cl)

}