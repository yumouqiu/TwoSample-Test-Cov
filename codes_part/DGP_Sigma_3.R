DGP_Sigma_3<-function(n1,n2,p,beta,r_1,r_2,r_3,prop_1,prop_2)
{
  n=n1*n2/(n1+n2)
  
  #if(r!=0)
  {
    Sigma_1=matrix(rep(0,p^2),nrow=p)
    for(i in 1:p)
    {
      Sigma_1[i,i]=0.6
    }
    Sigma_2=Sigma_1
    q=p*(p+1)/2
    count = floor(q^(1-beta))
    num = 0
    for(i in 1:p)
    {
      if(i%%2 == 1)
      {
      num = num + 1
      if(num <= floor(count*prop_1))
      {
      Sigma_2[(i+1),i] = Sigma_2[(i+1),i] + sqrt(4*r_1*log(p)/n)
      Sigma_2[i,(i+1)] = Sigma_2[(i+1),i]
      }
      if( (num> floor(count*prop_1))& (num<= floor(count*(prop_1+prop_2))))
      {
        Sigma_2[(i+1),i] = Sigma_2[(i+1),i] + sqrt(4*r_2*log(p)/n)
        Sigma_2[i,(i+1)] = Sigma_2[(i+1),i]
      }
      if( (num>floor(count*(prop_1+prop_2)))& (num<=count))
      {
        Sigma_2[(i+1),i] = Sigma_2[(i+1),i] + sqrt(4*r_3*log(p)/n)
        Sigma_2[i,(i+1)] = Sigma_2[(i+1),i]
      }
      }
    }
    
  }
  list(Sigma_1,Sigma_2)
  
  
}


