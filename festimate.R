#' Find the estimation of coefficents for high dimensional compositional data, L1 and L2 penalized functions are used 
#' When using this code, we expect a linear model, and the sum of the coefficietns equals 0.
#' @author Pan Wang (panwang2012@u.northwestern.edu)
#' @references This method is a modification of the work from Lin, Shi, Feng and Li (2014).Variable selection in regression with compositional covariates. Biometrika (2014), pp. 1–13
#' @param yy A vector of response variable y
#' @param Z A matrix of independent variables. It should be log of the compositional data matrix, in which each row should sum to 1, and each component should be between 0 and 1.
#' @param betastart Starting vector for the estimation of coefficient beta.
#' @param alphastart Starting vector for the estimation of lagrange multiplier for sum of the coefficients.
#' @param ac1 Accuracy when getting an estimation of parameters. High accuracy means more accuracy results.
#' @param ac2 Accuracy when getting an estimation of parameters. High accuracy means more accuracy results.
#' @param lambda1 Coefficient for L1 penlized function
#' @param lambda2 Coefficient for L2 penlized function
#' @param mu Multiplier of the square of the sum of the coefficients in Lagrange method.
#' @return Estimation of coefficients
#' @export 
#' 
festimate<-function(yy,Z,betastart,alphastart,ac1,ac2,lambda1,lambda2,mu)
{
  pp=length(Z[1,])
  nn=length(Z[,1])
  beta1=rep(betastart,pp)
  beta2=beta1+100
  alpha1=alphastart
  alpha2=alphastart+100
  vv=rep(0,pp)
  for(kk in 1:pp)
    vv[kk]=sum(Z[,kk]*Z[,kk])/nn ##squared sum of each column
  
  while(sqrt(sum((beta1-beta2)*(beta1-beta2))+(alpha1-alpha2)^2)>ac2)
  {
    beta2=beta1+100   
    alpha2=alpha1
    while(sqrt(sum((beta1-beta2)*(beta1-beta2)))>ac1)
    {
      beta2=beta1
      for(jj in 1:pp)
      {
        tempj=0
        for(i in 1:nn)
        {
          tempj=tempj+1/nn*(yy[i]-sum(Z[i,]*beta1)+Z[i,jj]*beta1[jj])*Z[i,jj]
        }
        tj=tempj-mu*(sum(beta1)-beta1[jj]+alpha2)
        beta1[jj]=1/(vv[jj]+2*lambda2+mu)*fsto(tj,lambda1)
      }    
    }
    alpha2=alpha1
    alpha1=alpha2+sum(beta1)      
  }
  
  vec=c(beta1,betastart,alphastart,ac1,ac2,lambda1,lambda2,mu)
  write.csv(vec,"temptempest.csv")
  file.append("tempest.csv", "temptempest.csv") 
  
  return(beta1)  
}
