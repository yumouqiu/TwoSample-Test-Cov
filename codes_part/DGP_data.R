###################################################################################################################
##           This function is used to generate the data under H_0 (r=0) and H_1 (r!=0)
##                    
##                                    INPUTS 
##         Sigma_1,Sigma_2: two population covariance matrices 
##                   n1,n2: sample sizes for the two populations
##                      p : dimension of covariance
##                   beta : signal sparsity level 
##                     r  : signal strength level
##                   dist : "Gaussian"---Gaussian distribution
##                          "Gamma"----Gamma distribution
##   
##                                   OUTPUTS                                    
##          
##                    X1 : n1*p data matrix from the first population
##                    X2 : n2*p data matrix from the second population 
#####################################################################################################################

DGP_data<-function(Sigma_1,Sigma_2,n1,n2,p,beta,r,dist)
{
  library(MASS)
   
  # generate the samples such that X_i=\Gamma_1 Z_i
  if(r==0)
  {
    t=eigen(Sigma_1, symmetric=T)
    ev=sqrt(t$values)  # square root of eigenvalues
    G=as.matrix(t$vectors)
    D=G*0; for(i in 1:p) D[i,i]=ev[i]
    Gamma1=G%*%D%*%t(G)    # M1=var(y_t)^{1/2}
    Gamma2=Gamma1
  }
  if(r!=0)
  {
    t=eigen(Sigma_2, symmetric=T)
    delta=abs(min(min(t$values),0))+0.01
    I=matrix(rep(0,p^2),nrow=p)
    for(i in 1:p)
    {
      I[i,i]=1
    }
    
    Sigma_2=Sigma_2+delta*I
    t=eigen(Sigma_2, symmetric=T)
    ev=sqrt(t$values)  # square root of eigenvalues
    G=as.matrix(t$vectors)
    D=G*0; for(i in 1:p) D[i,i]=ev[i]
    Gamma2=G%*%D%*%t(G)    # M1=var(y_t)^{1/2}
    
    Sigma_1=Sigma_1+delta*I
    t=eigen(Sigma_1, symmetric=T)
    ev=sqrt(t$values)  # square root of eigenvalues
    G=as.matrix(t$vectors)
    D=G*0; for(i in 1:p) D[i,i]=ev[i]
    Gamma1=G%*%D%*%t(G)    # M1=var(y_t)^{1/2}
  }   

  X1=NULL
  X2=NULL
  for(i in 1:n1)
  {
  	if(dist == "Gaussian")
  	{Z1 = rnorm(p)}else{
    Z1 = rgamma(p,4,2)
    Z1 = (Z1-mean(Z1))/sd(Z1)}
    X1 = rbind(X1,t(Gamma1%*%Z1))
 }
  
  for(i in 1:n2)
  {
  	if(dist == "Gaussian"){
    Z2 = rnorm(p)}else{
    Z2 = rgamma(p,4,2)
    Z2 = (Z2-mean(Z2))/sd(Z2)}
    X2=rbind(X2,t(Gamma2%*%Z2))
  }
 
  
  list(X1,X2)
  
}