    #SCRIPT FOR EXPLORING ELASTICITY PATTERNS


library(popbio)     #for solving matrices
library(triangle)
bwt=7.5763e-06;awt=3.0848 #length-weight parameters (Kohler et al 1996)


#Input pars
#LH=1
LH=2
if(LH==1)
  {
    RangeA.1=c(40,60)
    a.1=RangeA.1[2]
    age.mat=c(9,17)
    Rangefec=2:17
  }

if(LH==2)
  {
    Long=70
    RangeA.1=c(Long,round(1.3*Long))
    a.1=RangeA.1[2]
    Rangefec=6:10
    prop.mat=0.38
    age.mat=round(prop.mat*RangeA.1)
  }


first.age=0

Range.reprod_cycle=1:3
prob.rep.cycle=c(0.05, 0.25, 0.70)
sexratio=0.5

k.mean=0.056
Linf.mean=0.0719
Lo=0.014
sigma=matrix(c(0.0002147362, -0.0002039886,-0.0002039886,  0.0001978148),ncol=2,byrow=T)
colnames(sigma)=rownames(sigma)=c("Linf.est","K.est")

Aver.T=14
sd.Aver.T=1.5



    #elasticities
juvenile=(first.age):(age.mat[1]-1)
adult1=age.mat[1]:age.mat[2]
adult2=(age.mat[2]+1):a.1

names.elast.stage=c("juv.surv","adult1.surv","adult2.surv",
                    "juv.fec","adult1.fec","adult2.fec")



n.elas.stages=length(names.elast.stage)



# Natural mortality function
M.fun=function(A,k,Linf,Aver.T,age.mat,bwt,awt,W.G,W.noG)
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
  
  
  
  #STEP 2. get mean at age
  nat.mort=data.frame(m.Jensen.2,m.Pauly,m.Hoenig,m.PetWro,m.Lorenzen)
  N=ncol(nat.mort)
  nn=length(age)
  
  Weights=data.frame(W.J=rep(W.noG,nn),W.P=rep(W.G,nn),W.H=rep(W.noG,nn),
                     W.P=rep(W.noG*1,nn),W.L=rep(W.noG*1,nn))
  
  nat.mort=data.frame(nat.mort,Weights)
  nat.mort=apply(nat.mort, 1, function(x) weighted.mean(x[1:N], x[(N+1):ncol(nat.mort)],na.rm = T))
  
  
  return(nat.mort)
}

#Loop
iterations=10
ELAST=vector('list',length=length(iterations))
A.SIM=LAMBDA=rep(NA,iterations)
for (s in 1:iterations)
{
  A.sim=ceiling(rtriangle(1, a=RangeA.1[1], b=RangeA.1[2], c=RangeA.1[1]))
  #number of age classes
  age=first.age:A.sim
  
  if(LH==1)Age.mat.sim=ceiling(rtriangle(1, a=age.mat[1], b=age.mat[2], c=mean(c(age.mat[1],age.mat[2]))))
  if(LH==2)Age.mat.sim=round(0.38*A.sim)
  
  #growth
  Linf.sim=Linf.mean*100
  Lo.sim=Lo*100
  
  #total length
  total.length=Lo.sim+(Linf.sim-Lo.sim)*(1-exp(-k.mean*age))
  
  #fork length
  fork=(0.9442*total.length*100)-5.7441 # (Kohler et al 1996, length in cm)
  mid.FL.fem=fork
  
  
  #natural mortality
  m.sim=M.fun(A=A.sim,k=k.mean,Linf=Linf.mean,Aver.T=Aver.T,age.mat=Age.mat.sim,
              bwt,awt,W.G=1,W.noG=1)
  
  Fec.sim=round(rtriangle(1,a=Rangefec[1],b=Rangefec[length(Rangefec)],c=round(mean(Rangefec))))
  Rep.cy.sim=sample(Range.reprod_cycle,1,replace=T,prob = prob.rep.cycle)
  
  #reproductive cycle
  Pmax.sim=1/Rep.cy.sim
  
  #survivorship
  S=exp(-m.sim)         
  
  #proportion surviving
  lx=rep(NA,length(age))
  lx[1]=1.0
  for (i in 2:(length(age)))lx[i]=lx[i-1]*S[i]
  
  #reproductive schedules   
  Meanfec.sim=rep(Fec.sim,length(age))
  MF=c(rep(0,(Age.mat.sim-1)),Meanfec.sim[Age.mat.sim:length(Meanfec.sim)])
  mx=MF*sexratio*Pmax.sim
  
  #probability of surviving (for birth-pulse, post-breeding census)
  px=vector(length=length(lx))
  for(i in 2:length(lx)) px[i-1]=(lx[i])/(lx[i-1])
  
  
  #fertility  (for birth-pulse, post-breeding census)
  bx=mx*px
  
  
  #projection matrix
  PX=px
  PX=PX[-length(PX)]
  BX=bx
  n=length(BX)
  Data=matrix(0,nrow=n,ncol=n)
  diag(Data)[-nrow(Data)]=PX
  Data=rbind(matrix(0,nrow=1,ncol=n),Data)
  Data=Data[-(n+1),]
  Data[1,]=BX
  rownames(Data)=colnames(Data)=(first.age+1):n
  Proyec.matrix=Data
  
  #lambdas and elasticities
  LAMBDA[s]=lambda(Proyec.matrix)
  ELAST[[s]]=elasticity(Proyec.matrix)
  A.SIM[s]=A.sim
}



#Extract average elasticities per life stage
Elast.stage=matrix(,ncol=n.elas.stages,nrow=length(ELAST))
nn=1:n.elas.stages
for (e in 1:length(ELAST))
{
  test=ELAST[[e]]
  
    #survival
  Elast.stage[e,nn[1]]=sum(test[2:nrow(test),juvenile])
  Elast.stage[e,nn[2]]=sum(test[2:nrow(test),adult1])
  Elast.stage[e,nn[3]]=sum(test[2:nrow(test),adult2[adult2<=A.SIM[e]]])
  
    #fecundity
  Elast.stage[e,nn[4]]=sum(test[1,juvenile])
  Elast.stage[e,nn[5]]=sum(test[1,adult1])
  Elast.stage[e,nn[6]]=sum(test[1,adult2[adult2<=A.SIM[e]]])
  

}
colnames(Elast.stage)=names.elast.stage 


table.Elast=colMeans(Elast.stage, na.rm = T)
table.Elast=table.Elast/sum(table.Elast)  #standardise
RangeAge.group=rep(c("Juv","Ad1","Ad2"),2)


barplot(table.Elast, ylim=c(0,max(table.Elast)*1.1),mgp = c(2.5, 0.6, 0),
        names.arg= RangeAge.group,xlab="",ylab="",
        axis.lty=1, axes=T,col=2:4,cex.names=1.05,las=1,cex.axis=1.2,cex.lab=1.65,space=c(.0,.99))
box()
abline(v=3.99, col =1,lty=2)
legend('topleft','survival',bty='n',cex=1.25)
legend('topright','fecundity',bty='n',cex=1.25)
mtext("Age group",side=1,line=2,font=1,las=0,cex=1.3)
mtext("Elasticities",side=2,line=2.2,font=1,las=0,cex=1.3)



