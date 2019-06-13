  #-Instantaneous natural mortality schedule (years-1)

storeM.sim=storePars.sim=storeObs.sim=NULL
  for (i in 1:100) 
  {
    A.sim=round(rtriangle(1, a=RangeA[1], b=A, c=45))    #triangular dist
    growth.pars.sim=rmvnorm(1,mean=c(k.mean,Linf.mean,to.mean),sigma=sigma)    #multivariate normal dist
    if(growth.pars.sim[,2]<6 | growth.pars.sim[,2]>8 | growth.pars.sim[,1]<=0)    #repeat until sensible pars obtained
          { repeat 
              {
                growth.pars.sim=rmvnorm(1,mean=c(k.mean,Linf.mean,to.mean),sigma=sigma)
                if(growth.pars.sim[,2]>6 & growth.pars.sim[,2]<8)break
              }
          }
    k.sim=growth.pars.sim[1]
    Linf.sim=growth.pars.sim[2]
    to.sim=growth.pars.sim[3]
  

    age.mat.sim=round(rtriangle(1,a=age.mat[1],b=age.mat2[2],c=floor(mean(c(age.mat,age.mat2)))))    
    Aver.T.sim=rnorm(1,Aver.T,sd.Aver.T)    #normal dist
  
  
  
  M=function(A,k,Linf,Aver.T,age.mat,to)
  {
    #.number of age classes
    age=1:A      
  
    #.age invariant
      #(Jensen 1996)
    m.Jensen.1=1.5*k
    m.Jensen.1=rep(m.Jensen.1,length(age))
    m.Jensen.2=1.65/age.mat
    m.Jensen.2=rep(m.Jensen.2,length(age))
  
        #(Pauly 1980)  
    m.Pauly=10^(-0.0066-0.279*log10(Linf*100)+0.6543*log10(k)+0.4634*log10(Aver.T))
    m.Pauly=rep(m.Pauly,length(age))
      
        #(Hoenig 1983, combined teleost and cetaceans)    
    m.Hoenig=exp(1.44-0.982*log(A))      
    m.Hoenig=rep(m.Hoenig,length(age))
  
    #.age dependent
         #(Peterson and Wroblewski 1984, weight in grams)
    total.length=Linf*(1-exp(-k*(age-to)))
    fork=(0.9442*total.length*100)-5.7441 # (Kohler et al 1996, length in cm)
    wet.weight=(7.5763e-06)*fork^3.0848 # (Kohler et al 1996)
    m.PetWro=1.92*(wet.weight*1000)^-0.25
  
          #(Chen & Watanabe 1989)
    TM=(-1/k)*log(1-exp(k*to))+to
    ao=1-exp(-k*(TM-to))
    a1=k*exp(-(k*(TM-to)))
    a2=-0.5*(k)^2*exp(-(k*(TM-to)))
    m.Chen=ifelse(age<=TM,k/(1-exp(-k*(age-to))),k/(ao+a1*(age-TM)+a2*((age-TM)^2)))

  
      #.sample M at age from uniform distribution with bounds min and max M at age
    nat.mort=cbind(m.Jensen.1,m.Jensen.2,m.Pauly,m.Hoenig,m.PetWro,m.Chen)
    min.mort=apply(nat.mort, 1, min)
    max.mort=apply(nat.mort, 1, max)
  
    #observed mortality at age
    nat.mort=runif(length(age),min.mort,max.mort)

    
    #estimate mortality at age
    fn_ob=function(theta)
    {
      #predicted mortality
      m_pred=theta[1]*(age)^-theta[2]
      #residual
      epsilon=log(nat.mort)-log(m_pred)
 
    #total negative log likelihood
      nloglike=-1.0 *sum(dnorm(epsilon,0,theta[3]*m_pred,log=T))  

      return(nloglike)
    }
    fit_ob=optim(theta,fn_ob,method="BFGS",hessian=T)
  
    #predicted mortality at age
    nat.mort.pred=fit_ob$par[1]*(age)^-fit_ob$par[2]

    return(list(M=nat.mort.pred,fit=fit_ob$par,obs=nat.mort))
  }
  m.sim=M(A.sim,k.sim,Linf.sim,Aver.T.sim,age.mat.sim,to.sim)$M     #uniform dist    
  Pars.sim=M(A.sim,k.sim,Linf.sim,Aver.T.sim,age.mat.sim,to.sim)$fit
    m.obs=M(A.sim,k.sim,Linf.sim,Aver.T.sim,age.mat.sim,to.sim)$obs  
  storeM.sim=rbind(storeM.sim,m.sim)
  storePars.sim=rbind(storePars.sim,Pars.sim)  
    storeObs.sim=rbind(storeObs.sim,m.obs)
    
}    
colss=1:100
    plot(storeM.sim[1,],col="white")
for (i in 1:nrow(storeM.sim))
{
    lines(storeM.sim[i,],col=colss[i])
}

#MANUAL 
plot(storeObs.sim[1,],ylim=c(0,max(storeObs.sim)))
for(i in 1:nrow(storeObs.sim)) points(storeObs.sim[i,])

lines(theta[1]*(age)^-theta[2],col=2,lwd=3)
lines(0.28*(age)^-0.25,col=3,lwd=3)
lines(0.2*(age)^-0.45,col=4,lwd=3)

par1=rnorm(1000,0.23,0.009)   #parameters for m calculation
par2=rnorm(1000,0.35,0.0275)

for(i in 1:1000)
{
  lines(lines(par1[i]*(age)^-par2[i],col=3,lwd=2))
}  