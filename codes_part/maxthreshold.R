###################################################################################################################
##           This function is used to obtain the bootstrap resamples for testing statistic
##                    
##                                    INPUTS
##                   Sigma : estimated common covariance matrix by using data samples from the two populations 
##                      p  : the dimension of covariance matrix
##                     n1  : sample size from the first population
##                     n2  : sample size from the second population 
##                     s_0 : thresholding lower bound
##                     eta : thresholding upper bound 1-eta
##   
##                                   OUTPUTS                                    
##          
##                    T  : the bootstrap resample of testing statistic MTT
##                   
#####################################################################################################################


maxthreshold<-function(Sigma,p,n1,n2,s_0,eta){
  
  size1 = n1
  size2 = n2	
	
  t=eigen(Sigma, symmetric=T)
  ev=sqrt(t$values)  # square root of eigenvalues
  G=as.matrix(t$vectors)
  D=G*0; for(i in 1:p) D[i,i]=ev[i]
  Gamma=G%*%D%*%t(G)    # M1=var(y_t)^{1/2}
  
  
  Z1<-matrix(rnorm(size1*p,0, 1),ncol=p)
  Z2<-matrix(rnorm(size2*p,0, 1),ncol=p)
  sam1=t(Gamma%*%t(Z1))
  sam2=t(Gamma%*%t(Z2))             #  Generate the bootstrap data for the two populations
 
  source("intermediate.R")
  bootres = intermediate(sam1,sam2)
  M = bootres[[3]]                  # Obtain the bootstrap resample M_ij
  
  s=M/(4*log(p))
  s=s[(s>s0)]
  s=s[(s<1-eta)]
  s=c(s,1-eta)
  T=NULL
  for(i in 1:length(s))
  {
    lambdas=sqrt(4*s[i]*log(p))
    phis=dnorm(lambdas)
    Phis=1-pnorm(lambdas)
    mus=p*(p+1)*(lambdas*phis+Phis)
    sigmas=p*(p+1)*(((lambdas)^3+3*lambdas)*phis+3*Phis)
    T0=sum(M[M>4*s[i]*log(p)])
    T=c(T,(T0-mus)/sqrt(sigmas))
  }
  T=max(T)                         # Obtain the bootstrap resample of the maximum of standardized thresholding statistic
  return(T)                 
}
