###################################################################################################################
##           This function is used to calculate the estimated covariance matrix and M_ij
##                    
##                                    INPUTS 
##                    sam1 : n1*p data matrix from the first population
##                    sam2 : n2*p data matrix from the second population 
##   
##                                   OUTPUTS                                    
##          
##                    CovX : estimated p*p covariance matrix for the first population
##                    CovY : estimated p*p covariance matrix for the second population
##                     M   : p*(p+1)/2 vector M_ij
#####################################################################################################################

intermediate<-function(sam1,sam2)
{
  
  n1 = nrow(sam1)
  p = ncol(sam1)
  n2 = nrow(sam2)
   
  CovX = cov(sam1)*(n1-1)/n1
  CovY = cov(sam2)*(n2-1)/n2
  
  Xc = apply(X, 2, function(x) x-mean(x))  # center
  Yc = apply(Y, 2, function(x) x-mean(x))
  theta1 = (t(Xc^2) %*% Xc^2)/n1 - CovX^2  # Theta matrix for X
  theta2 = (t(Yc^2) %*% Yc^2)/n1 - CovY^2  # Theta matrix for Y
  M = (sig1-sig2)^2/(theta1/n1+theta2/n2)  # M matrix
  
  list(CovX, CovY, M[upper.tri(M, diag=T)])
  
}