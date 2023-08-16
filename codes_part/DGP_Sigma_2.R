###################################################################################################################
##           This function is used to generate the covariance matrix (Design 2) under H_0 (r=0) and H_1 (r!=0)
##                    
##                                    INPUTS 
##                   n1,n2: sample sizes for the two populations
##                      p : dimension of covariance
##                   beta : signal sparsity level 
##                     r  : signal strength level
## 
##                                   OUTPUTS                                    
##          
##                   covariance matrices Sigma_1 and Sigma_2 
#####################################################################################################################

DGP_Sigma_2<-function(n1,n2,p,beta,r)
{
  set.seed(12345)
  D_0=matrix(rep(0,p^2),nrow=p)
  for(i in 1:p)
  {
    D_0[i,i]=sqrt(runif(n=1,min=0.1,max=1))
  }
  # First step: to generate the two Sigma_1 and Sigma_2
  Sigma_1=matrix(rep(0,p^2),nrow=p)
  k_1=floor(p/4)
  
  for(i in 1:p)
  {
    Sigma_1[i,i]=1
    for(j in 1:p)
    {
      if(i!=j)
      { 
        for(k_0 in 1:k_1)
        {
          if((4*(k_0-1)+1<=i)&&(i<=4*k_0))
          {
            
            if((4*(k_0-1)+1<=j)&&(j<=4*k_0))
            {
              Sigma_1[i,j]=0.5
            }
          }
        }  
      }
    }
  }
  Sigma_1=D_0%*%Sigma_1%*%D_0
  Sigma_2=Sigma_1
  n=n1*n2/(n1+n2)
 
  ## under the alternative, add the difference matrix U 
  U=matrix(rep(0,p^2),nrow=p)
  if(r!=0)
  {
    s=sqrt(4*r*log(p)/n)  
    q=p*(p+1)/2
    m_p=floor(q^(1-beta)/2)
    k_0=floor(m_p/p)
    k_1=m_p-p*k_0+k_0*(k_0+1)/2
    
    for(l in 1:k_1)
    {
      U[l+k_0+1,l]=s
      U[l,l+k_0+1]=s
    }
    if(k_0>=1)
    {
      for(k in 1:p)
      {
        for(l in 1:p)
        {
          if(k!=l)
          {
            if(abs(k-l)<=k_0)
            {
              U[k,l]=s
              U[l,k]=s
            }
          }
        }
      }
    }
 
    Sigma_2=Sigma_2+U
    
  }
  list(Sigma_1,Sigma_2)
 
  
}
  
   
 