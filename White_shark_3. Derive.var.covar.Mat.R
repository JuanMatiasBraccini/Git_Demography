library(mvtnorm)


#growth pars from o'connor 2011
K=0.056
Linf=7.19
Linf=Linf/100
Lo=1.4
Lo=Lo/100

rough.dat=data.frame(age=c(0,0,0,1,rep(2,9),rep(3,11),rep(4,9),rep(5,6),rep(6,8),7,7,7,7,8,8,9,9,11,11,11,11,
                           13,13,15,19,20,20),TL=c(110,130,160,155,140,155,165,166,170,180,190,220,230,170,180,
                          200,210,220,225,240,260,265,270,270,170,180,190,220,230,240,250,290,291,180,220,270,290,
                          340,350,220,230,240,250,310,320,330,340,250,340,380,400,350,370,280,400,385,395,450,460,
                          440,485,450,505,505,485)/(100*100))

Start.pars=list(Linf.est=jitter(Linf,10),K.est=jitter(K,10))  
fit=nls(TL~Lo+(Linf.est-Lo)*(1-exp(-K.est*age)),data=rough.dat,start=Start.pars)
SIGMA=vcov(fit)
aa=rmvnorm(1000,mean=c(K,Linf),sigma=SIGMA)
par(mfcol=c(2,1))
hist(aa[,1],xlab="k",main=paste("CV=",round(sd(aa[,1])/mean(aa[,1]),3)))
hist(aa[,2]*100,xlab="Linf (m)",main=paste("CV=",round(sd(aa[,2])/mean(aa[,2]),3)))


