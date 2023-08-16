###################################################################################################################
##           This function is used to obtain the testing statistic MTT and its theoretical pvalue
##                    
##                                    INPUTS
##          
##                      p  : the dimension of covariance matrix
##                      M  : the vector of M_ij 
##                     s_0 : thresholding lower bound
##                     eta : thresholding upper bound 1-eta
##   
##                                   OUTPUTS                                    
##          
##                     T  :  the maximum of standardized thresholding statistic MTT
##                  pvalue:  the theoretical pvalue of MTT based on the limiting distribution
#####################################################################################################################

proposed<-function(p,M,s0,eta)
{
  library(VGAM)
  
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
  T=max(T)
  a=sqrt(2*log(log(p)))
  b=2*log(log(p))+log(log(log(p)))/2-log(pi)/2+log(1-s0-eta)
  pvalue=1-pgumbel(a*T-b)
  list(T,pvalue)
  
}