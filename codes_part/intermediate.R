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
   
  Xmeans=colMeans(sam1)
  Ymeans=colMeans(sam2)
  
  n1=nrow(sam1)
  p=ncol(sam1)
  n2=nrow(sam2)
   
  CovX = cov(sam1)*(n1-1)/n1
  CovY = cov(sam2)*(n2-1)/n2
  
  M = NULL
  
  for(i in 1:p)
  {
  	for(j in i:p)
  	{
  		theta1 = 0
  		theta2 = 0
  		for(t1 in 1:n1)
  		{
  			C_diff = sam1[t1,] - Xmeans
  			theta = C_diff[i]*C_diff[j]-CovX[i,j] 
  			theta1 = theta1 + theta*theta
  		}
  		for(t1 in 1:n2)
  		{
  			C_diff = sam2[t1,] - Ymeans
  			theta = C_diff[i]*C_diff[j]-CovY[i,j] 
  			theta2 = theta2 + theta*theta
  		}
  		theta1 = theta1/n1
  		theta2 = theta2/n2
  		diff = CovX[i,j] - CovY[i,j]
  		M = c(M,diff*diff/(theta1/n1+theta2/n2))
  	}
  }
  
  list(CovX,CovY,M)
  
}