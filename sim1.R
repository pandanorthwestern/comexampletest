#' Codes for testing the preformaance of festimate.R using simulation
#' @references These functions below follows the simulation setting from Lin, Shi, Feng and Li (2014).Variable selection in regression with compositional covariates. Biometrika (2014), pp. 1–13
#' @description These functions are basic functions for simulation procedure, return compositional dataset, response vector y, the value of information criterion GIC. 
#' @author Pan Wang (panwang2012@u.northwestern.edu)

#' @description Generate the mean vector for multinormal distribution.
#' @param p Dimension of the data.
#' @param ter Length of nonzero components.
#' @param co co*p is the mean of the nonzero components
#' @return The mean vector for multinormal distribution.
#' @export 
#' 

ftheta<-function(p,ter,co)
{
  ttheta=rep(0,p)
  for(ii in 1:ter)
  {
    ttheta[ii]=log(co*p)
  }
  
  return(ttheta)
}

#' @description Generate correlation matrix for multinormal distribution.
#' @param p Dimension of the data.
#' @param rho cor(xi,xj)=rho^abs(i-j).
#' @return Correlation matrix.
#' @export 

fSIGMA<-function(rho,p)
{
  repzero<-rep(0,p*p)
  SIGMA<-matrix(data = repzero, nrow = p, ncol = p)
  for(ii in 1:p)
  {
    for(jj in 1:p)
      SIGMA[ii,jj]=rho^(abs(ii-jj))
  }
  return(SIGMA)
}

#' @description Generate compositional dataset, and apply log() on each component
#' @param p Dimension of the data. 
#' @param m Number of samples.
#' @return Log of compositional dataset.
#' @export 

fZ<-function(m,p,ter,co,rho)
{
  SIGMAtemp=fSIGMA(rho,p)
  thetatemp=ftheta(p,ter,co)
  library(mvtnorm)
  wtemp=rmvnorm(m, mean = thetatemp, sigma = SIGMAtemp)##m rows and p columns
  rrepzero<-rep(0,m*p)
  xtemp<-matrix(data = rrepzero, nrow = m, ncol = p)
  for(ii in 1:m)
  {
    xexptemp=exp(wtemp[ii,])
    sumexptemp=sum(xexptemp)
    for(jj in 1:p)
      xtemp[ii,jj]=xexptemp[jj]/sumexptemp
  }
  ztemp<-log(xtemp)
  return(ztemp)
}

#' @description Generate response vector.
#' @param zk Log of componsional dataset.
#' @param beta Vector of coefficient. Sum(beta)=0.
#' @param mean mean vector when generating composional dataset.
#' @param sd standard error for error term in linear regression model.
#' @return Response vector.
#' @export 

fY<-function(zk,beta,mean,sd)
{
  q=length(zk[,1])
  epsilontemp=rnorm(q, mean , sd)
  y=zk%*%beta+epsilontemp
  return(y)  
}

#' @description Function of soft thresholding operator.
#' @param t A number of input.
#' @param lambda Threshold.
#' @return Number of fsto(t).
#' @export 
fsto<-function(t,lambda)
{
  if(abs(t)<=lambda)
    return(0)
  else
  {
    if(t>0)
      return(t-lambda)
    else
      return(t+lambda)
  }    
}


#' @description Function for getting the value of information criterion GIC.
#' @param beta Coefficient vector beta.
#' @param Z Log of componsional dataset.
#' @return y response vector.
#' @export 
GIC<-function(beta,y,Z)
{
  ng=length(y)
  numzero=length(which(beta==0))
  lenbeta=length(beta)
  tempsigmasq=0
  for(i in 1:ng)
  {
    tempsigmasq=tempsigmasq+1/ng*(y[i]-sum(Z[i,]*beta))^2   
  }
  GICn=log(tempsigmasq)+(lenbeta-numzero-1)*log(log(ng))/ng*log(max(lenbeta,ng))
  return(GICn)
}
