
#------------------White Shark 36 or 70 Demographic Analysis-------------------------

#note: this script runs birth-pulse, post-breeding census Leslie Matrix analyses for white shark
#     to simply compare effect of setting max age to 36 vs 108 and size at maturity to vs 30. 


rm(list=ls(all=TRUE))
library(triangle)
library(mvtnorm)  		#for multivariate normal pdf
library(MASS)       #for fitting distributions and bivariate kernel densities
library(vioplot)    #for violin plots
library(plotrix)      		#needed for graph legends
library(popbio)     #for solving matrices
library(matrixStats) #for colQuantiles function
library(TeachingDemos)      #for grey scale

handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')

setwd(handl_OneDrive("Analyses/Demography/White shark 36 or 70"))



#---DATA SECTION---

 


#---PARAMETERS SECTION---

NAMES=c("Mollet","Smith","Hamady")
Scen.names=1:3

  #first age class
first.age=0  

  #Maximum age
A.Smith=36	
A.Mollet=60
A.new=70     #Hamady et al 2013
A.new=mean(A.new)
A=c(A.Smith,A.Mollet,A.new)

Prop.age.mat=0.38  #from meta analysis of mackerel sharks (see excel spreadsheet)

  #Age maturity
age.mat.Smith=9
age.mat.Mollet=15
age.mat.new=round(A.new*Prop.age.mat)

#elasticities
names.elast.stage=c("juv.surv","adult.surv","juv.fec","adult.fec") 
n.elas.stages=length(names.elast.stage)

  #Natural mortality (Hoenig 1983, combined teleost and cetaceans)
M.Smith=exp(1.44-0.982*log(A.Smith))
M.Mollet=exp(1.44-0.982*log(A.Mollet))
M.new=exp(1.44-0.982*log(A.new))



  #Fecundity pars used by the 3 demographic studies
Fec.Smith=7
Fec.Mollet=8.9


Fec=round(mean(c(Fec.Smith,Fec.Mollet)))



#length of reproductive cycle (in years)
Reprod_cycle=2         


#pups sex ratio
sexratio=0.5          



#---PROCEDURE SECTION---



#---MAIN SECTION---

scenarios=list()
scenarios[[1]]=c(A.Mollet,age.mat.Mollet,M.Mollet)
scenarios[[2]]=c(A.Smith,age.mat.Smith,M.Smith)
scenarios[[3]]=c(A.new,age.mat.new,M.new)
names(scenarios)=NAMES



#Leslie function
Leslie.fn=function(scen)
{
  #1. Set scenarios
    #biological scenario
  A.sim=scen[1]
  age=first.age:A.sim


  #reproduction
  Reprod_cycle.sim=Reprod_cycle
  Pmax.sim=1/Reprod_cycle.sim
  Meanfec.sim=Fec
  age.mat.sim=scen[2]
  StartRep.sim=age.mat.sim+1  #key parameter (Smith et al 1998), set at 1 year after mature

  #natural mortality
  m.sim=scen[3]
                   
  #survivorship
  S=exp(-m.sim)         
  
  #proportion surviving
  lx=rep(NA,length(age))
  lx[1]=1.0
  for (i in 2:(length(age)))lx[i]=lx[i-1]*S
  
  #reproductive schedules
  fec_age=ifelse (age>=StartRep.sim, Meanfec.sim,0)                       
  mx=fec_age*sexratio*Pmax.sim

  #probability of surviving (for birth-pulse, post-breeding census)
  px=vector(length=length(lx))
  for(i in 2:length(lx)) px[i-1]=(lx[i])/(lx[i-1])
  
  #fertility  (birth-pulse, post-breeding census)
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
  
    
  #3. Solve projection matrix
  LAMBDA=lambda(Proyec.matrix)
  r.sim=log(LAMBDA)         
  t2.sim=log(2)/r.sim           #pop doubling time
  v=reproductive.value(Proyec.matrix)    #reproductive value
  Elasticities=elasticity(Proyec.matrix) #elasticities

  
  #4. Extract average elasticities per life stage
  juvenile=(first.age):(StartRep.sim-1)
  adult1=StartRep.sim:A.sim
  Elast.stage=1:n.elas.stages
  names(Elast.stage)=names.elast.stage
  Elast.stage[1]=sum(Elasticities[2:nrow(Elasticities),juvenile])
  Elast.stage[2]=sum(Elasticities[2:nrow(Elasticities),adult1])
  Elast.stage[3]=sum(Elasticities[1,juvenile])
  Elast.stage[4]=sum(Elasticities[1,adult1])
  
  
  return(list(r=r.sim,t2=t2.sim,Elast.stage=Elast.stage,LAMBDA=LAMBDA,v=v,Proyec.matrix=
    Proyec.matrix))
}


#run Leslie function for each life history scenario 
output=vector('list',length=length(scenarios))
names(output)=names(scenarios)

system.time(for(l in 1:length(scenarios))
{
  output.by.year=Leslie.fn(scenarios[[l]])
  output[[l]]=output.by.year
})


 


#---REPORT SECTION---

  # Table 1
#Population parameters
table1=function(data,data1)
{
  Amax=data1[1]
  Amat=data1[2]
  NatMort=data1[3]
  R=data$r
  T2=data$t2
  data=rbind(Amax,Amat,NatMort,R,T2)
  return(data)
}
TEST=matrix(nrow=5,ncol=length(NAMES))
for(i in 1:length(output)) TEST[,i]=table1(output[[i]],scenarios[[i]])
colnames(TEST)=names(output)
TEST=as.data.frame(TEST)
TEST$Variables=c("Max age", "Age mat","M","r","T2")
write.table(TEST,"Table1.csv",sep=",",row.names=F)




#Plot population projections
n.Yr.prol=50
scenario.N0=.2
U.frac=0.1
MeanM=mean(c(scenarios[[1]][3],scenarios[[2]][3],scenarios[[3]][3]))
scenario.U=rep(U.frac*MeanM,n.Yr.prol)


Risk.fn=function(scenario,Proyec.matrix)
{
  Rel.abundance=matrix(NA,ncol=n.Yr.prol,nrow = 1)
  A.sim=scenario[1]
  Proj.Mat=Proyec.matrix
    
  H=diag(nrow(Proj.Mat))
  diag(H)=1-(scenario.U[1]) 
  MH=Proj.Mat%*%H
  
  
  #Project population into future
  No=matrix(stable.stage(Proj.Mat)) 
  Total.No=sum(No)
  n.No=length(No)
  
    #add density-depence in survival
  if(dens.dep=="YES")
  {
    #1. calculated the density-dependent factor   
    fun.dens=function(dens)
    {
      #choose projection matrix
      Mat=Proj.Mat
      
      #add survival density-dependence to projection matrix
      
        #over all age classes, effect of all age classes
      if (WhatClass=="ALL")for (nn in 2:nrow(Mat)) Mat[nn,(nn-1)]=Mat[nn,(nn-1)]*exp(-dens*Total.No)
      
        #over all age classes, effect of older age classes
      if (WhatClass=="ALL.bigger")for (nn in 2:(nrow(Mat)-1))Mat[nn,(nn-1)]=Mat[nn,(nn-1)]*exp(-dens*sum(No[(nn+1):n.No]))  
          
        #over first age class, effect of all age classes
      if (WhatClass=="ONE")for (nn in 2:2) Mat[nn,(nn-1)]=Mat[nn,(nn-1)]*exp(-dens*Total.No)
      
      
      #calculate lambda
      LAMBDA=lambda(Mat)
      
      #obj function
      epsilon=(1-LAMBDA)^2
      
      return(epsilon)
    }
    
    fit.dens=optimize(fun.dens,lower=0.0001,upper=100)
    Dens.factor=fit.dens$minimum
  }
  
   
  #exploited population size
  n.vec=vector("list",length = n.Yr.prol)
  n.vec[[1]]=(No)*scenario.N0
  for (p in 2:n.Yr.prol)
  {
    #no density dependence
    if(dens.dep=="NO") n.vec[[p]]=MH%*%matrix(n.vec[[p-1]])
    
    #Density-depence in survival
    if(dens.dep=="YES")
    {
      dummy=MH          
      if (WhatClass=="ALL")for (nn in 2:nrow(dummy))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
      if (WhatClass=="ALL.bigger")for (nn in 2:(nrow(dummy)-1))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]][(nn+1):n.No]))             #NEW 
      if (WhatClass=="ONE")for (nn in 2:2)dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
      n.vec[[p]]=dummy%*%matrix(n.vec[[p-1]])
    }
  }
  
  Pop.size=rep(0,n.Yr.prol)
  for(y in 1:n.Yr.prol) Pop.size[y]=sum(n.vec[[y]])
    
    
    #6. Store relative abundance 
    Rel.abundance=Pop.size
  
  
  return(Rel.abundance=Rel.abundance)
}  


#Control density dependence
dens.dep="YES"
#dens.dep="NO"

WhatClass="ALL"
#WhatClass="ALL.bigger"
#WhatClass="ONE"

Store.abundance=vector('list',length=length(scenarios))
for(i in 1:length(scenarios)) Store.abundance[[i]]=Risk.fn(scenarios[[i]],output[[i]]$Proyec.matrix)
MAX.Y=1
NN=length(scenarios)
line.type=c(1,2,1)
#line.type=rep(1,NN)
line.col=c("black","grey20","grey60")
#line.col=2:4

tiff(file="Figure1.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
plot(1:n.Yr.prol,seq(0,MAX.Y,length.out=n.Yr.prol),col='transparent',xaxt='n',yaxt='n',xlab="",ylab="",ylim=c(0,MAX.Y))
for(i in 1:length(scenarios))lines(1:n.Yr.prol,Store.abundance[[i]],type='l',
                                   lty=line.type[i],col=line.col[i],lwd=4)
axis(1,at=seq(1,n.Yr.prol,by=1),labels=F,cex.axis=1.3,las=1,tck=-0.015)
axis(2,at=seq(0,MAX.Y,by=MAX.Y/5),labels=seq(0,MAX.Y,by=MAX.Y/5),cex.axis=1.3,las=1,tck=-0.03)
axis(1,at=seq(5,n.Yr.prol,by=5),labels=seq(5,n.Yr.prol,by=5),cex.axis=1.3,las=1,tck=-0.03)

mtext("Projected year",side=1,outer=T,line=-3,font=1,las=0,cex=1.75)
mtext("Relative population size",side=2,outer=T,line=-1.45,font=1,las=0,cex=1.75)
legend("topleft",c("Scenario 1","Scenario 2","Scenario 3"),bty='n',lty=line.type,
       col=line.col,lwd=3,cex=1.5)
dev.off()







# #___________NOT USED___________#  
#   #Figure 1
# COL=c("black","gray55","black")
# LINE=c(1,1,2)
# # COL=c("black","gray55","black","gray55")
# # LINE=c(1,1,2,2)
# 
#   #Plot population parameters
# tiff(file="Figure2.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
# par(mfcol=c(1,1),las=1,mai=c(.6,.6,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
# #reproductive value
# plot(first.age:max(A),first.age:max(A), type='l',col="white",xlab="",ylim=c(0,14),xlim=c(first.age,max(A)),
#      ylab="",cex.axis=1.5)
# for (i in 1:length(A))
# {
#   agE=first.age:A[i]
#   lines(agE,output[[i]]$v,lwd=2,col=COL[i],lty=LINE[i])
# }
# 
# mtext("Age",side=1,line=2,font=1,las=0,cex=1.75)
# mtext("Reproductive value",side=2,line=1.8,font=1,las=0,cex=1.75)
# legend("topleft",NAMES,col=COL,lty=LINE,bty="n",cex=1.5,lwd=2)
# dev.off()
# 
# 
# #elasticities
# tiff(file="Figure3.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
# par(mfcol=c(3,1),las=1,mai=c(.6,.6,.05,.01),omi=c(.001,.01,.001,.01),mgp=c(1.5,1,0))
# for(i in 1:length(A))
# {
#   barplot(output[[i]]$Elast.stage,ylim=c(0,1),xlab="",ylab="",names.arg= "",yaxt='n')
#   text(4.35,.2,NAMES[i],cex=1.5)
#   abline(v=2.6, col =1,lty=2)
#   axis(2,at=seq(.2:1,by=.2),F)
#   axis(1,at=c(0.7, 1.9,3.2, 4.3),F)
#   if(i %in%c(1))
#     {
#       text(1.35,.9,'Survival',cex=1.75)
#       text(3.75,.9,'Fecundity',cex=1.75)
#      }
#   if(i %in%c(3)) axis(1,at=c(0.7, 1.9,3.2, 4.3),c("Juvenile","Adult","Juvenile","Adult"),cex.axis=1.75)
#   axis(2,at=seq(.2:1,by=.2),seq(.2:1,by=.2),cex.axis=1.5)
#   box()
# }
# 
# mtext("Age group",side=1,line=-1.6,font=1,las=0,cex=1.5,outer=T)
# mtext("Elasticities",side=2,line=-1.6,font=1,las=0,cex=1.5,outer=T)
# dev.off()
# 
# 
# 
