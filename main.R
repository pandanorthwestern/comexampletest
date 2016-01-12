#' Main function of simulation.
#' @details Should source sim1.R and festimate.R. 
#' @author Pan Wang (panwang2012@u.northwestern.edu)
#' @param veclambda1 a vector of candidate coefficients for L1 penlized funtion.
#' @param veclambda2 a vector of candidate coefficients for L2 penlized funtion.
#' @param m number of samples.
#' @param p dimension of each sample.
#' @return None. Output a csv file with the result summary.
#' @export 
#' 

GICtest_select<-function(veclambda1,veclambda2,m,p,betastart=1,alphastart=1,ac1=0.0001,ac2=0.0001/2,mu=1,ter=5,co=0.5,rho=0.5,setbeta=c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),mean=0,sd=0.5)
{  
  repzero=rep(0,(p-length(setbeta)))
  betaorigin=c(setbeta,repzero)
  Z=fZ(m,p,ter,co,rho)
  yy=fY(Z,betaorigin,mean,sd)
  len1=length(veclambda1)
  len2=length(veclambda2)
  GICrecord=100
  for(i in 1:len1)
    for(j in 1:len2)
    {
      lambda1=veclambda1[i]
      lambda2=veclambda2[j]
      betatemp=festimate(yy,Z,betastart,alphastart,ac1,ac2,lambda1,lambda2,mu)
      GICtemp=GIC(betatemp,yy,Z)
      if(GICtemp<GICrecord)
      {
        GICrecord=GICtemp
        betarecord=betatemp
        lambda1record=lambda1
        lambda2record=lambda2
      }
    }
  
  temp_sq_dif_beta=sum((betarecord-betaorigin)*(betarecord-betaorigin))
  temp_L1_dif_beta=sum(abs(betarecord-betaorigin))
  temp_pre_error=0
  
  ###another independent test sample
  Z=fZ(m,p,ter,co,rho)
  yy=fY(Z,betaorigin,mean,sd)
  
  for(r in 1:m)
    temp_pre_error=temp_pre_error+(yy[r]-sum(Z[r,]*betarecord))^2*1/m
  
  vec=c(m,p,betarecord,lambda1record,lambda2record,mu,GICrecord,m,p,temp_sq_dif_beta,temp_pre_error,temp_L1_dif_beta)
  
  write.csv(vec,"comptemp.csv")
  file.append("comp.csv", "comptemp.csv") 
  
}

