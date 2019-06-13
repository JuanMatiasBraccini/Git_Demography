Selec.original=Selectivity.sim
anios=1:length(Selec.original)
#Selectivity.sim=dnorm(1:length(Selec.original),length(Selec.original)/2,10)

Selectivity.sim=1*(1+exp((-log(19)*((anios-30)/((0.95*max(anios))-30)))))^-1

Selectivity.sim=Selectivity.sim/max(Selectivity.sim)

#Eff.prop=c(0.1,0.5,0.9)
#par(mfcol=c(3,3),las=1,mar=c(2,4.2,.15,0.5), oma=c(1.5,2,.5,0.5))

  for (i in 1:3)
  {
    


fit.F=optimize(fn.F,lower=log(0.0001),upper=log(3),A=A.sim,m=m.sim,Selectivity=Selectivity.sim,Linf=Linf.sim,k=k.sim,to=to.sim,Pmax=Pmax.sim,L50=L50.sim,L95=L95.sim,StartRep=StartRep.sim,Rangefec=Meanfec.sim)   
  f.sim=fishMort(A.sim,fit.F$minimum,m.sim,Selectivity.sim,Linf.sim,k.sim,to.sim,Pmax.sim,L50.sim,L95.sim,
                 StartRep.sim,Meanfec.sim)$f
  f.crash=f.sim
iu=1-exp(-f.crash)


  
  #3. Calculate F.current as proportion of F.crash
#   Eff.crash=runif(1,Eff.crash.range[1],Eff.crash.range[2])  #sample of effort crash
#   Eff.p=runif(1,Eff.p.range[1],Eff.p.range[2])              #sample of effort protected
#   Eff.prop=Eff.p/Eff.crash                                  #proportional effort
  
  f.sim=f.crash*Eff.prop[i]

  
  #4. Estimate population parameters   
  fn=function(A,log.r,log.f,m,Selectivity,Linf,k,to,Pmax,L50,L95,StartRep,Rangefec)demographic(A,log.r,log.f,m,Selectivity,
              Linf,k,to,Pmax,L50,L95,StartRep,Rangefec)$epsilon
	fit=optimize(fn,lower=log(0.0001),upper=log(0.99999),A=A.sim,log.f=f.sim,m=m.sim,Selectivity=Selectivity.sim,
               Linf=Linf.sim,k=k.sim,to=to.sim,Pmax=Pmax.sim,L50=L50.sim,L95=L95.sim,StartRep=StartRep.sim,Rangefec=Meanfec.sim) 	

  r.estim=demographic(A=A.sim,fit$minimum,f.sim,m.sim,Selectivity.sim,Linf.sim,k.sim,to.sim,Pmax.sim,L50.sim,
                      L95.sim,StartRep.sim,Meanfec.sim)$r



plot(1:A.sim,Selectivity.sim,type="l",xlab="age",ylab="selectivity")
text(max(A.sim)*0.7,(max(Selectivity.sim)*0.8),paste("eff.prop= ",round(Eff.prop[i],2),sep=""),cex=0.9)
text(max(A.sim)*0.7,(max(Selectivity.sim)*0.6),paste("r= ",round(r.estim,4),sep=""),cex=0.9)
text(max(A.sim)*0.7,(max(Selectivity.sim)*0.9),paste("u.crash= ",round(iu,2),sep=""),cex=0.9)
text(max(A.sim)*0.7,(max(Selectivity.sim)*1),paste("f.crash= ",round(f.crash,2),sep=""),cex=0.9)
text(max(A.sim)*0.7,(max(Selectivity.sim)*0.7),paste("f.sim= ",round(f.sim,2),sep=""),cex=0.9)

}

