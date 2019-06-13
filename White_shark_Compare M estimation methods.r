M.dummy.fun=function(A,k,Linf,Aver.T,age.mat,bwt,awt)
{
  #STEP 1. calculate M from different methods (see Kenchington 2013)
  #.age invariant
  #Jensen (1996)
  m.Jensen.1=1.5*k
  m.Jensen.1=rep(m.Jensen.1,length(age))
  
  m.Jensen.2=1.65/age.mat
  m.Jensen.2=rep(m.Jensen.2,length(age))
  
  #Pauly (1980)  
  m.Pauly=10^(-0.0066-0.279*log10(Linf*100)+0.6543*log10(k)+0.4634*log10(Aver.T))
  m.Pauly=rep(m.Pauly,length(age))
  
  #Hoenig (1983), combined teleost and cetaceans    
  m.Hoenig=exp(1.44-0.982*log(A))      
  m.Hoenig=rep(m.Hoenig,length(age))
  
  
  #.age dependent
  #Peterson and Wroblewski 1984 (weight in grams)
  wet.weight=(7.5763e-06)*fork^3.0848 
  m.PetWro=1.92*(wet.weight*1000)^-0.25
  
  #Lorenzen 1996 (weight in grams)
  m.Lorenzen=3*(wet.weight*1000)^-0.288
  
  #Gislason et al (2010) (lengths in grams)
  m.Gislason=1.73*(mid.FL.fem^-1.61)*((Linf*100)^1.44)*k
  if(m.Gislason[1]>1)m.Gislason=rep(NA,length(age))
  
  
  #STEP 2. get mean at age
  nat.mort=data.frame(m.Jensen.2,m.Pauly,m.Hoenig,m.PetWro,m.Lorenzen,m.Gislason)

  
  return(nat.mort)
}

NNN=100

store.m.sim=store.m.sim1=vector('list',length=NNN)

for (i in 1:NNN)
{
  age=first.age:ASim[i]
  #total length
  total.length=Lo.sim+(LinfSim[i]-Lo.sim)*(1-exp(-kSim[i]*(age)))
  
  #fork length
  fork=(0.9442*total.length*100)-5.7441 # (Kohler et al 1996, length in cm)
  mid.FL.fem=fork
  
  
  store.m.sim[[i]]=M.dummy.fun(A=ASim[i],k=kSim[i],Linf=LinfSim[i],Aver.T=TempSim[i],age.mat=AgeMatSim[i],
                               bwt,awt)
  
  
  age=first.age:ASim.1[i]
  #total length
  total.length=Lo.sim+(LinfSim.1[i]-Lo.sim)*(1-exp(-kSim.1[i]*(age)))
  
  #fork length
  fork=(0.9442*total.length*100)-5.7441 # (Kohler et al 1996, length in cm)
  mid.FL.fem=fork
  store.m.sim1[[i]]=M.dummy.fun(A=ASim.1[i],k=kSim.1[i],Linf=LinfSim.1[i],Aver.T=TempSim[i],age.mat=AgeMatSim.1[i],
                                bwt,awt)
}

par(mfcol=c(3,3))
for (i in 1:9)
{
  plot(store.m.sim1[[i]][,1],type='l',ylim=c(0,.25),main="LH1")
  lines(store.m.sim1[[i]][,2],col=2)
  lines(store.m.sim1[[i]][,3],col=3)
  lines(store.m.sim1[[i]][,4],col=4)
  lines(store.m.sim1[[i]][,5],col=5)
  lines(store.m.sim1[[i]][,6],col=6)
}
legend('topright',names(store.m.sim1[[i]]),bty='n',lty=1,col=1:6)

for (i in 1:9)
{
  plot(store.m.sim[[i]][,1],type='l',ylim=c(0,.25),main="LH2")
  lines(store.m.sim[[i]][,2],col=2)
  lines(store.m.sim[[i]][,3],col=3)
  lines(store.m.sim[[i]][,4],col=4)
  lines(store.m.sim[[i]][,5],col=5)
  lines(store.m.sim[[i]][,6],col=6)
}
legend('topright',names(store.m.sim[[i]]),bty='n',lty=1,col=1:6)