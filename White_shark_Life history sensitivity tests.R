#SENSITIVITY OF LIFE HISTORY PARS

#run main script for bringing in parameter values

#---MC simulations----

LH1.pars=list(one=c(40,9,3),two=c(60,9,3),three=c(40,17,3),four=c(60,17,3),
              five=c(40,9,17),six=c(60,9,17),seven=c(40,17,17),eight=c(60,17,17))

LH2.pars=list(one=c(70,27,6),two=c(91,27,6),three=c(70,35,6),four=c(91,35,6),
              five=c(70,27,10),six=c(91,27,10),seven=c(70,35,10),eight=c(91,35,10))


#1. Equilibrium conditions


store.lambda=vector('list',length(LH1.pars))
names(store.lambda)=names(LH1.pars)
  
#select life history scenario
l=1
#l=2

for (p in 1:length(LH1.pars))
{
  
  #1. Set scenarios
  #biological scenario
  scenario=Life.hist.scenarios[[l]]
  
  if(scenario==1)DATA=LH1.pars[[p]]
  if(scenario==2)DATA=LH2.pars[[p]]
  
  
  #1. Draw parameter samples and build functions
  #maximum age
  if(scenario==1)A.sim=DATA[1]
    if(scenario==2)A.sim=DATA[1]
    
    k.sim=0.05495554
  Linf.sim=7.434663
  Lo.sim=Lo
  
  
  #number of age classes
  age=first.age:A.sim
  
  #reproductive cycle
  Reprod_cycle.sim=2
  Pmax.sim=1/Reprod_cycle.sim
  
  
  #fecundity and age at maturity
  if(scenario==1) 
  {
    Meanfec.sim=DATA[3]
      age.mat.sim=DATA[2]
  }
  if(scenario==2) 
  {
    Meanfec.sim=DATA[3]
      age.mat.sim=DATA[2]  
  }
  
  
  #temperature
  Aver.T.sim=14
  
  
  #total length
  total.length=Lo.sim+(Linf.sim-Lo.sim)*(1-exp(-k.sim*age))
  
  #fork length
  fork=(0.9442*total.length*100)-5.7441 # (Kohler et al 1996, length in cm)
  mid.FL.fem=fork
  
  
  #natural mortality
  m.sim=M.fun(A=A.sim,k=k.sim,Linf=Linf.sim,Aver.T=Aver.T.sim,age.mat=age.mat.sim,
              bwt,awt,W.G=1,W.noG=1)
  
  
  #survivorship
  S=exp(-m.sim)         
  
  #proportion surviving
  lx=rep(NA,length(age))
  lx[1]=1.0
  for (i in 2:(length(age)))lx[i]=lx[i-1]*S[i]
  
  #reproductive schedules   
  #mx=Meanfec.sim*sexratio*Pmax.sim *Pred.mat
  MF=c(rep(0,(age.mat.sim-1)),rep(Meanfec.sim,length(age)-(age.mat.sim-1)))
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
  
  
  
  #3. Solve projection matrix
  LAMBDA=lambda(Data)
  store.lambda[[p]]=round(LAMBDA,4)
}


LH1.pars=as.data.frame(matrix(unlist(LH1.pars),nrow=8,byrow=T))
LH2.pars=as.data.frame(matrix(unlist(LH2.pars),nrow=8,byrow=T))
names(LH1.pars)=names(LH2.pars)=c("Max age","Age Mat","Num pups")

LH1.pars$Lambda=unlist(store.lambda)
LH2.pars$Lambda=unlist(store.lambda)


LH1.pars=LH1.pars[order(LH1.pars$Lambda),]
LH2.pars=LH2.pars[order(LH2.pars$Lambda),]

write.csv(LH1.pars,"LH1.pars.csv")
write.csv(LH2.pars,"LH2.pars.csv")