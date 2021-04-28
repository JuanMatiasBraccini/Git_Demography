
#------------------White Shark Demographic Analysis-------------------------

#note: this script runs birth-pulse, post-breeding census Leslie Matrix analyses for white shark,
      # sampling from biological pars distributions and projecting population into future based 
      # on different scenarios of harvest rate, selectivities and population sizes 

#scenarios: two scenarios on biol pars; two scenarios of selectivity; four reconstructed catch scenarios; 
#            discrete and distribution scenarios on No (within calculated No bounds)

#INDEX       #---1.  MODELLED SCENARIOS SECTION---
             #---2.  DATA SECTION---
             #---3.  PARAMETERS SECTION---             
             #---4.  FUNCTIONS SECTION---  
             #---5.  LESLIE MATRICES SECTION---  
             #---6.  CATCH RECONSTRUCTIONS SECTION---  
             #---7.  RISK ANALYSES (Population projections) SECTION---  
             #---8.  SENSITIVITY TESTS---
             #---9.  DISPLAY RESULTS SECTION---  
             #---10. FIGURES FOR STEVE SECTION---
             #---11. FUTURE PROJECTIONS STEVE SECTION---
             #---12. PREVIOUS CODE SECTION---

rm(list=ls(all=TRUE))

set.seed(999)  #for reproducibility


#---1.  MODELLED SCENARIOS SECTION---
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/MS.Office.outputs.R"))

#Life history scenarios
Life.hist.scenarios=list(LH.1=1, LH.2=2)

#Initial female population sizes
N.0=list(c(1500,3000),c(3000,5000),c(5000,7500),c(7500,10000),c(10000,20000))    #lower and upper bounds for uniform distribution derived from Ne (see "Get_No.R")
names(N.0)=paste(N.0)
N.0.Rick=c(1500,3000,5000,7500,10000)  #Rick Scenarios

N1975.scenarios=list(N1=.9)  #Population size in 1975


#Catch scenarios
KTCH.scen=paste("C",c(1:4),sep="")

#Post capture survival
PCS.scenarios=c("NO","YES","YES.100")


#---2.  DATA SECTION---

library(triangle)
library(mvtnorm)  		#for multivariate normal pdf
library(MASS)       #for fitting distributions and bivariate kernel densities
library(vioplot)    #for violin plots
library(plotrix)      		#needed for graph legends
library(popbio)     #for solving matrices
library(matrixStats) #for colQuantiles function
library(TeachingDemos)      #for grey scale
#library(xlsReadWrite)
library(grDevices)    #grey scale
library(gplots)

source(handl_OneDrive("Analyses/Demography/White shark/3. Derive.var.covar.Mat.R"))
setwd(handl_OneDrive("Analyses/Demography/White shark/Outputs"))

#Effort

  # 1. GHATF effort (Gummy gillnet Effort 1975:2012, in 1000 km lift)
Effort.GHATF= read.csv(handl_OneDrive("Data/Population dynamics/GHATS_gillnet_effort.csv")) #1971:2014 financial yrs 
id=which(Effort.GHATF$Year==1975)
Effort.GHATF$SA[1:2]=0
Effort.GHATF$Tas[1:2]=0

    #1.1 Fill in missing effort by region with proportional SA effort for 2006
#Effort.GHATF$SA[37:42]=0.49*Effort.GHATF$Total[37:42] 
#Effort.GHATF$Vic[37:42]=0.37*Effort.GHATF$Total[37:42]
#Effort.GHATF$Tas[37:42]=0.14*Effort.GHATF$Total[37:42]

    #1.2 Use SA and West BS effort only (i.e. effort proportin west 147 E)
#Effort.GHATF.prop.BS=read.csv("C:/Matias/Data/Population dynamics/West.BS.eff.prop.csv") 

    #1.3. Total effort used in model
Effort.GHATF$Total=Effort.GHATF$SA
#Effort.GHATF$Total=Effort.GHATF$Vic*Effort.GHATF.prop.BS$West.Eff.prop+Effort.GHATF$Tas/2+Effort.GHATF$SA

Effort.GHATF.prev=Effort.GHATF[1:(id-1),c(1,5)]   
Effort.GHATF=Effort.GHATF[id:nrow(Effort.GHATF),c(1,5)]

    #1.4 Add effort previous to available data
Recons.Effort.GHATF=data.frame(Year=1940:1970,Total=0)
Effort.GHATF=rbind(Recons.Effort.GHATF,Effort.GHATF.prev,Effort.GHATF)


  #2. SA Marine gillnet scalefish (in 1000 km gillnet days)
SA.scale= read.csv(handl_OneDrive("Data/Population dynamics/MSFeffortbyfinancialyear.csv")) #financial yr
SA.scale$Year=as.numeric(substr(SA.scale$Financial.year,1,4))
names(SA.scale)[2]="Effort"

SA.scale.period= read.csv(handl_OneDrive("Data/Population dynamics/MSFeffortbysurveytimeperiod.csv")) 
SA.scale.period=SA.scale.period[2:4,]
names(SA.scale.period)=c("Period","Agg.Eff")   


  #3. WA Shark gillnet         (in km gillnet days)
Effort.WA= read.csv(handl_OneDrive("Analyses/Demography/White shark/Effort.white.csv"))  
Effort.WA[,2]=Effort.WA[,2]/1000  #convert to 1000 km gillent days



  #add effort previous to available data
Effort.WA.1940=data.frame(Year=1940:1953,Total=0)    #0 before 1955
Effort.WA.1955=data.frame(Year=1954:1974,Total=0)    # interpolation to 1955-1974
Effort.WA.1955$Total=Effort.WA.1955$Year*0.2139-417.99   #from excel equation
Effort.WA.1955$Total[1]=0
Effort.WA=rbind(Effort.WA.1940,Effort.WA.1955,Effort.WA)

#Combine efforts
Ef.net=data.frame(Year=Effort.GHATF[,1],Total=Effort.GHATF[,2]+Effort.WA[,2])

#check efforts
par(mfcol=c(1,1),las=1)
plot(Effort.WA$Year,Effort.WA$Total,type='l',ylim=c(0,48),xlab="Year",ylab="Effort (1000 km gn days)",lwd=2,
     cex.lab=1.25,col=3)
#lines(Ef.net[,1],Ef.net[,2],col=1,lwd=2)
lines(Effort.GHATF[,1],Effort.GHATF[,2],col=2,lwd=2)
lines(Effort.GHATF[,1],c(rep(0,length(1940:1982)),SA.scale$Effort),col=4,lwd=2)
legend("topleft",c("TDGDLF","GHATF (SA only)","SA marine"),col=c(3,2,4),lty=1,bty="n",lwd=2,cex=1.25)
#legend("topleft",c("Total","TDGDLF","GHATF (SA only)","SA marine"),col=c(1,3,2,4),lty=1,bty="n",lwd=2,cex=1.25)



#4. Reconstructed TDGDLF catch
# Recons.obs.catch=read.csv("C:/Matias/Analyses/White shark catch reconstruction/Final.catch.obs.csv")
# Recons.pred.catch=read.csv("C:/Matias/Analyses/White shark catch reconstruction/Final.catch.missing.csv")
# Recons.catch=rbind(Recons.obs.catch,Recons.pred.catch)
# Recons.catch=aggregate(cbind(Mean.Catch.App1,Low95.Catch.App1,Up95.Catch.App1)~Time,Recons.catch,sum)

Recons.catch=read.csv(handl_OneDrive("Analyses/White shark catch reconstruction/Total.catch.csv"))
Recons.catch2=read.csv(handl_OneDrive("Analyses/White shark catch reconstruction/Total.catch.Method2.csv"))

names(Recons.catch)[4:6]=c("Mean.Catch.App1","Low95.Catch.App1","Up95.Catch.App1")
names(Recons.catch2)[4:6]=names(Recons.catch)[4:6]



#Reconstructed TDGDLF cpue for calculating Commonwealth catch
Recons.cpue=read.csv(handl_OneDrive("Analyses/White shark catch reconstruction/CPUE.Folly.csv"))
Recons.cpue2=read.csv(handl_OneDrive("Analyses/White shark catch reconstruction/CPUE.Folly.Method2.csv"))

names(Recons.cpue)[c(4,6:7)]=c("Mean.cpue","Low95.cpue","Up95.cpue")
names(Recons.cpue2)[c(4,6:7)]=names(Recons.cpue)[c(4,6:7)]

Recons.cpue=subset(Recons.cpue,Fishing.region=="SS2")
Recons.cpue2=subset(Recons.cpue2,Fishing.region=="SS2")


#5. Reconstructed Lobster and TDGDLF droplines catches for WA
Droplines= read.csv(handl_OneDrive("Data/Population dynamics/Dropline.csv")) 
Droplines$Year=1988:2001


  #sex ratio of catch
Ktch.sex.ratio=0.5


#Scenarios
  #gillnet selectivity
    # Figure 12, Malcolm et al 2001, this is only gillnet fisheries
select.size=seq(1,7,by=.5)
num.rel.alive=c(1,10,10,4,3,7,9,2,1,2,0,0,0)
num.dead=c(0,7,16,6,6,2,6,7,11,4,1,1,1)
select.num=num.rel.alive+num.dead

    # Table 15, Malcolm et al 2001, this is all sources of mortality
select.Robin=data.frame(prop=c(0.2882,0.2620,0.2140,0.1616,0.0655,0.0087),
                        length=c("0-2","2-3","3-4","4-5","5-6","6+"))
prop.sel=select.Robin[,1]

  #gillnet selectivity Steve
prop.sel.Steve=c(0.267,0.267,0.233,0.189,0.044,0)   #from Steve's catch reconstruction          


  #check selectivities
par(mfcol=c(1,1))
barplot(rbind(prop.sel,prop.sel.Steve),beside=T,names.arg=select.Robin$length,
        legend.text = c("Scen1", "Scen2"))


Sel.scenarios=list(sel.1=1,sel.2=2)     

mean.prop.sel=colMeans(rbind(prop.sel,prop.sel.Steve)) 


prop.sel.2=c(0.0725,0.1075,0.13,0.16,0.23,0.30)  #assume wider selectivity for scenario 2



#Maturity data (Francis 1996)
Fem.maturity=data.frame(Size.class=c("300-325","325-350","350-375","375-400","400-425","425-450"
                                     ,"450-475","475-500","500-525","525-550","550-575"),
                        mid.point=c(313,338,363,388,413,438,463,488,513,538,563),
                        n.imm=c(1,0,3,5,3,0,1,1,0,0,0),n.mat=c(1,0,0,0,3,2,3,1,4,1,2))
Fem.maturity$prop.mat=with(Fem.maturity,n.mat/(n.imm +n.mat))
Fem.maturity$prop.mat=ifelse(is.na(Fem.maturity$prop.mat),0,Fem.maturity$prop.mat)



#---3.  PARAMETERS SECTION---


#3.1. Life history

#source: mostly Bruce 2008
first.age=0   #first age class. Assuming post-breeding census (this is the 0-1 age class actually)

A.Smith=36  #Maximum age (60 years unverified, 40-50 years more reasonable, Bruce 2008)
A.Mollet=60
A.Hamady=40
Longevity.Hamady=70
A.1=A.Hamady		       
A.2=Longevity.Hamady+Longevity.Hamady*.3  #Cortes 2002
RangeA.1=c(A.Hamady,A.Mollet)
RangeA.2=c(A.Hamady,A.2)
max.size=c(6,7)    #estimated maximum total length (in metres)
bwt=7.5763e-06;awt=3.0848 #length-weight parameters (Kohler et al 1996)
Dry.w=0.2   #convertion wet to dry weight (Cortes 2002)

  #mortality pars
Aver.T=14   #average water temperature (Reynolds sea surface temperature for 2009)
sd.Aver.T=1.5

	#growth pars
#Use O'Connor 2011 (2 VBG par model) only

#LH1
k=0.056
Linf=7.19
Linf=Linf/100  #put in same order of magnitude as K for random samples
Lo=1.4 #m, TL

# n_age=c(21,114,51,79)		    #sample size age studies (corresponding to Cailliet et al 2005, Wintner and Cliff 1999,
#                           # Malcolm et al 2001, cited in Bruce 2008; and O'Connor 2011)
# Linf=c(7.6737,6.86,6.598,7.19)	#mean Linf of corresponding study (in metres)
# LinfSd=0.1*Linf		        #SD Linf  set to 10% (no SD available in Bruce 2008),10% chosen for having 
#                           # reasonable variability in log space
# k=c(0.058,0.065,0.071,0.056)    #mean k of corresponding study
# kSd=0.1*k		              #SD k set to 30% (no SD available in Bruce 2008)
# to=c(-3.53,-4.4,-2.33,-3.8)    #to of corresponding study
# toSd=0.1*abs(to)          #SD to 
# 
# #note: Wintner and Cliff reports PCL, convert to TL using TL= 1.251*PCL+5.207

#LH2
k.2= k                       
Linf.2= Linf 
Lo.2=Lo

  #reproductive pars 
    #reproductive cycle
Reprod_cycle=3        #length of reproductive cycle (in years) 
Range.reprod_cycle=1:Reprod_cycle   #range of possible reproductive cycles

    #litter size
Rangefec=2:10         #scenario 1
Rangefec2=2:10        #scenario 2

    #age maturity 1
age.mat=c(9,17)      #range age at maturity of females (Bruce 2008 and Smith et al 1998 combined). These studies do not 
                     # specified whether 50% or what

    #age maturity 2
Prop.age.mat=0.38  #from meta analysis of mackerel sharks (see excel spreadsheet)
age.mat2=round(Prop.age.mat*RangeA.2)

    #range size at maturity of females (in metres)
size.mat=c(4.5,5)     

    #pups sex ratio
sexratio=0.5          

  #classes for elasticity analyses              
#pup=first.age
juvenile.LH1=(first.age):(age.mat[1]-1)
juvenile.LH2=(first.age):(age.mat2[1]-1)
adult1.LH1=age.mat[1]:age.mat[2]
adult1.LH2=age.mat2[1]:age.mat2[2]
adult2.LH1=(age.mat[2]+1):RangeA.1[2]
adult2.LH2=(age.mat2[2]+1):RangeA.2[2]

names.elast.stage=c("juv.surv","adult1.surv","adult2.surv",
                    "juv.fec","adult1.fec","adult2.fec") 
n.elas.stages=length(names.elast.stage)

#   #maternity pars
# Max.Size=runif(1,min=max.size[1],max=max.size[2])
# L50=mean(size.mat)              #Arbitrary length at 50% maternal condition
# L50Sd=0.1*L50		                      #Arbitrary SD set to 10% of mean value
# L95=mean(c(size.mat,Max.Size))	#Arbitrary length at 95% maternal condition
# L95Sd=0.1*L95		                      #Arbitrary SD set to 10% of mean value

  #calculate SD for lognormal distributions
# logkSd=sqrt(log(1+(kSd/k)^2))  	          #SD in log space
# logLinfSd=sqrt(log(1+(LinfSd/Linf)^2))
# #logtoSd=sqrt(log(1+(toSd/to)^2))
# logLoSd=sqrt(log(1+(LoSd/Lo)^2))
#logL50Sd=sqrt(log(1+(L50Sd/L50)^2))		    
#logL95Sd=sqrt(log(1+(L95Sd/L95)^2))		    


  #Set up pars for creating distributions
     #growth
#LH1
k.mean=k
Linf.mean=Linf
Lo.mean=Lo


# growth.weights=n_age/sum(n_age)
# k.mean=weighted.mean(k,growth.weights)
# #logkSd.mean=weighted.mean(logkSd,growth.weights)
# Linf.mean=weighted.mean(Linf,growth.weights)
# #logLinfSd.mean=weighted.mean(logLinfSd,growth.weights)
# to.mean=weighted.mean(to,growth.weights)
# #toSd.mean=weighted.mean(toSd,growth.weights)

#LH2
k.mean.2=k.2
Linf.mean.2=Linf.2
Lo.mean.2=Lo.2

    #growth variance covariance matrix from approximating O'Connor 2011 variability
sigma=SIGMA

    #reproduction  
prob.rep.cycle=c(0.05,0.25,0.7)   #here tune the probs for 1-, 2- or 3-year reprod cycle
#sigma.mat=matrix(c((L50Sd/3)^2,0.025,0.025,(L95Sd/3)^2),ncol=2)   #variance covariance matrix for maternal ogive

    #.selectivity
cv=0.3  #assumed CV of vonB function
high.length=c(2,3,4,5,6,7) #upper and lower bounds of length classes
low.length=c(0,2,3,4,5,6)

#Reference Points
Biom.ref.point.1=.3      #population size in Yr.end > Blim (Braccini et al 2014)
Biom.ref.point.2=2*Biom.ref.point.1      #2 times the limit (conservative)
#Biom.ref.point.1=.5      #population size in Yr.end >= 50% No
#Biom.ref.point.2=.6     #population size in Yr.end =  60% No
#Biom.ref.point.3=.7     #population size in Yr.end =  70% No


#3.2 Catch reconstruction

#Proportion of vessels carrying droplines
TDGDLF.prop.hook=0.5
Lobster.prop.hook=0.1
PCS.white=0.5   #arbitrary alternative for risk assesment
#PCS.mako=0.25   #mako shark PCS gillnets (source Braccini et al 2012)


#3.3  Risk analyses

#Add density dependence or not
dens.dep="YES"
#dens.dep="NO"

#where is density dependence added
#WhatClass="ALL"   #applied to all age classes
WhatClass="ONE"  #applied to class one survival only (simil stock-recruitment function)

#Set N1998 to EPBC criteria or not
#Apply.EPBC="YES"
Apply.EPBC="NO"

#Select type of projection
#note: Stochastic projections use different projection matrices with dens. dep. compensation
#       calculated on the most productive matrix. This should avoid N2012 to exceed No
Proy.type="Stochas"   
#Proy.type="Determ"



#---4.  FUNCTIONS SECTION--- 

# Function for calculating Natural mortality
Include.Then_2015="YES"
if(Include.Then_2015=="NO")M.fun=function(A,k,Linf,Aver.T,age.mat,bwt,awt,W.G,W.noG)
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
  #Peterson and Wroblewski 1984 (dry weight in grams)
  wet.weight=(7.5763e-06)*fork^3.0848 
  m.PetWro=1.92*((wet.weight*Dry.w)*1000)^-0.25
  
  #Lorenzen 1996 (weight in grams)
  m.Lorenzen=3*(wet.weight*1000)^-0.288
  
  #Gislason et al (2010) (lengths in grams)
  m.Gislason=1.73*(mid.FL.fem^-1.61)*((Linf*100)^1.44)*k
  if(m.Gislason[1]>1)m.Gislason=rep(NA,length(age))
  
  
  #STEP 2. get mean at age
  nat.mort=data.frame(m.Jensen.2,m.Pauly,m.Hoenig,m.PetWro,m.Lorenzen,m.Gislason)
  N=ncol(nat.mort)
  nn=length(age)
  
  Weights=data.frame(W.J=rep(W.noG,nn),W.P=rep(W.G,nn),W.H=rep(W.noG,nn),
                     W.P=rep(W.noG*1,nn),W.L=rep(W.noG*1,nn),W.G=rep(W.G*1,nn))
  
  nat.mort=data.frame(nat.mort,Weights)
  nat.mort=apply(nat.mort, 1, function(x) weighted.mean(x[1:N], x[(N+1):ncol(nat.mort)],na.rm = T))
  
  
  return(nat.mort)
}
if(Include.Then_2015=="YES")M.fun=function(A,k,Linf,Aver.T,age.mat,bwt,awt,W.G,W.noG)
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
  
  #Then et al (2015)
  m.Then.1=4.899*A^(-0.916)
  m.Then.1=rep(m.Then.1,length(age))
  #m.Then.2=4.118*(k^0.73)*(Linf^(-0.33))
  
  #.age dependent
  #Peterson and Wroblewski 1984 (dry weight in grams)
  wet.weight=(7.5763e-06)*fork^3.0848 
  m.PetWro=1.92*((wet.weight*Dry.w)*1000)^-0.25
  
  #Lorenzen 1996 (weight in grams)
  m.Lorenzen=3*(wet.weight*1000)^-0.288
  
  #Gislason et al (2010) (lengths in grams)
  m.Gislason=1.73*(mid.FL.fem^-1.61)*((Linf*100)^1.44)*k
  if(m.Gislason[1]>1)m.Gislason=rep(NA,length(age))
  
  
  #STEP 2. get mean at age
  nat.mort=data.frame(m.Jensen.2,m.Pauly,m.Hoenig,m.Then.1,m.PetWro,m.Lorenzen,m.Gislason)
  N=ncol(nat.mort)
  nn=length(age)
  nat.mort.all=nat.mort
  
  Weights=data.frame(W.J=rep(W.noG,nn),W.P=rep(W.G,nn),W.H=rep(W.noG,nn),W.T=rep(W.noG,nn),
                     W.P=rep(W.noG*1,nn),W.L=rep(W.noG*1,nn),W.G=rep(W.G*1,nn))
  
  nat.mort=data.frame(nat.mort,Weights)
  nat.mort=apply(nat.mort, 1, function(x) weighted.mean(x[1:N], x[(N+1):ncol(nat.mort)],na.rm = T))
  
  
  return(list(nat.mort=nat.mort,nat.mort.all=nat.mort.all))
}

# Function for calculating Selectivity
Sel.fn=function(A,Linf,k)
{
  #Lower bound of age interval
  age.low=first.age:(A-1)
  
  #Upper bound of age interval
  age.up=(first.age+1):A
  
  #von B function
  mean.length=Lo+(Linf-Lo)*(1-exp(-k*age.low))
  
  sd=mean.length*cv
  logSd=sqrt(log(1+(cv)^2))
  
  #prob of age given length (lognormal)
  prob=function(data)
  {
    #pnorm(data,mean.length,sd)
    plnorm(data,log(mean.length),logSd)
    
  }
  
  
  high.bound=sapply(high.length,prob)
  low.bound=sapply(low.length,prob)
  mat.prob=high.bound-low.bound
  
  #normalised matrix
  st.mat.prob=t(t(mat.prob)/rowSums(t(mat.prob)))
  
  #calculate normalised selectivity
  if(scenario.sel==1)
  {
    Sel=rowSums(st.mat.prob*(matrix(rep(mean.prop.sel,nrow(st.mat.prob)),ncol=(ncol(st.mat.prob)),byrow=T)))
    Sel=Sel/max(Sel)
  }
  
  if(scenario.sel==2)
  {
    Sel=rowSums(st.mat.prob*(matrix(rep(prop.sel.2,nrow(st.mat.prob)),ncol=(ncol(st.mat.prob)),byrow=T)))
    Sel=Sel/max(Sel)
  }
  
  return(Sel=Sel)
}

# Function for filling in ages to max longevity for all iterations
add.missing.age=function(data,A)
{
  test=list()
  for(i in 1:length(data))
  {
    datos=data[[i]]
    if(length(datos)<A)
    {
      extra=A-length(datos)
      datos=c(datos,rep(NA,extra))
    }
    test[[i]]=datos
    if(length(datos)>A)test[[i]]=datos[1:A]
  }
  dummy=do.call(rbind,test)
  return(dummy)
}


# Function for drawing random sample
SamplE=function(Mean,SD)rlnorm(1,log(Mean),SD)

# Function for reconstructing catch
fn.catch=function(DAT,CPUE.zn2,ADD.PCS)
{
  #1. Resample catch within confidence limits
  #1.1. WA catch
  SE.period=(DAT$Up95.Catch.App1-DAT$Low95.Catch.App1)/3.92
  SD.Total.catch=SE.period/Recons.catch$Mean.Catch.App1  
  log.TC.Sd=sqrt(log(1+(SD.Total.catch)^2)) #sd in logspace  
  DAT$Mean.Catch.App1=mapply(SamplE,DAT$Mean.Catch.App1,log.TC.Sd)
  DAT=aggregate(Mean.Catch.App1~Time,DAT,sum)
  
  #1.2 Zn2 cpue
  SE.period=(CPUE.zn2$Up95.cpue-CPUE.zn2$Low95.cpue)/3.92
  SD.cpue=SE.period/CPUE.zn2$Mean.cpue  
  log.cpue.Sd=sqrt(log(1+(SD.cpue)^2)) #sd in logspace  
  CPUE.zn2$Mean.cpue=mapply(SamplE,CPUE.zn2$Mean.cpue,log.cpue.Sd)
  
  
  #2. Calculate annual catch for TDGDLF
  WA.Eff.per=merge(WA.Eff.Per,DAT[,c(1:2)],by.x="Period",by.y="Time")
  
  Eff.WA=merge(Effort.WA,WA.Eff.per,by="Period",all.x=T)
  Eff.WA=Eff.WA[order(Eff.WA$Year),]
  Eff.WA$Eff.prop=with(Eff.WA,Total/Agg.Eff)
  
  #catchability assumption
  Eff.WA$AnnualKtch=with(Eff.WA,Eff.prop*Mean.Catch.App1)             
  
  #fill in early period
  Eff.WA[1:48,]$AnnualKtch=Eff.WA[49,]$AnnualKtch*(Eff.WA[1:48,]$Total/Eff.WA[49,]$Total)
  TOT.Ktch=Eff.WA
  
  
  #3. Calculate annual catch for GHATF
  GHATF.Eff.Per$Mean.Catch.App1=CPUE.zn2$Mean.cpue*GHATF.Eff.Per$Agg.Eff
  Eff.GHATF=merge(Effort.GHATF,GHATF.Eff.Per,by="Period",all.x=T)
  Eff.GHATF=Eff.GHATF[order(Eff.GHATF$Year),]
  Eff.GHATF$Eff.prop=with(Eff.GHATF,Total/Agg.Eff)
  
  #catchability assumption
  Eff.GHATF$AnnualKtch=with(Eff.GHATF,Eff.prop*Mean.Catch.App1)             
  
  #fill in early period
  Eff.GHATF[1:48,]$AnnualKtch=Eff.GHATF[49,]$AnnualKtch*(Eff.GHATF[1:48,]$Total/Eff.GHATF[49,]$Total)
  TOT.Ktch$AnnualKtch.GHATF=Eff.GHATF$AnnualKtch
  
  #add 1938:1939
  TOT.Ktch=rbind(TOT.Ktch[1:2,],TOT.Ktch) 
  TOT.Ktch$Year[1:2]=1938:1939
  
  
  #4. SA marine scale fishery  
  SA.scale.period$Mean.Catch.App1=CPUE.zn2$Mean.cpue*SA.scale.period$Agg.Eff
  SA.scale.eff=merge(SA.scale,SA.scale.period,by="Period",all.x=T)
  SA.scale.eff=SA.scale.eff[order(SA.scale.eff$Year),]
  SA.scale.eff$Eff.prop=with(SA.scale.eff,Effort/Agg.Eff)
  
  #catchability assumption
  SA.scale.eff$AnnualKtch=with(SA.scale.eff,Eff.prop*Mean.Catch.App1)             
  
  #fill in early period
  SA.scale.eff[1:5,]$AnnualKtch=SA.scale.eff[6,]$AnnualKtch*(SA.scale.eff[1:5,]$Effort/SA.scale.eff[6,]$Effort)
  Catch.SA.marine.scale=SA.scale.eff$AnnualKtch
  
  #add missing years
  Catch.SA.marine.scale=c(rep(0,length(1938:1982)),Catch.SA.marine.scale)
  names(Catch.SA.marine.scale)=Yr.proj
  
  
  #5. Albany whaling station
  Catch.Albany.whale=c(rep(0,length(Yr.proj[1]:1971)),rep(10,length(1972:1978)),rep(0,length(1979:Yr.end)))
  
  
  #6. SA game fishers
  Catch.SA.game.fish=c(rep(3.2,length(Yr.proj[1]:1990)),rep(0,length(1991:Yr.end)))
  names(Catch.Albany.whale)=names(Catch.SA.game.fish)=Yr.proj
  
  
  #7. Add droplines from TDGDLF and Lobster fishery 
  Unif.sample=runif(1,0,1)
  DroplineKtch=Droplines$vesselTDGDLF*Unif.sample*TDGDLF.prop.hook+Droplines$vesselLobster*Unif.sample*Lobster.prop.hook
  Catch.Dropline=c(rep(0,length(Yr.proj[1]:1987)),DroplineKtch,rep(0,length(2002:Yr.end)))
  
  names(Catch.Dropline)=Yr.proj
  
  
  #8. PCS from gillnets
  if(ADD.PCS=="YES")
  {
    TOT.Ktch$AnnualKtch[match(1998:Yr.end,TOT.Ktch$Year)]=TOT.Ktch$AnnualKtch[match(1998:Yr.end,TOT.Ktch$Year)]*(1-PCS.white)
    TOT.Ktch$AnnualKtch.GHATF[match(1998:Yr.end,TOT.Ktch$Year)]=TOT.Ktch$AnnualKtch.GHATF[match(1998:Yr.end,TOT.Ktch$Year)]*(1-PCS.white)
    Catch.Dropline[match(1998:Yr.end,names(Catch.Dropline))]=Catch.Dropline[match(1998:Yr.end,names(Catch.Dropline))]*(1-PCS.white)
    Catch.SA.marine.scale[match(1998:Yr.end,names(Catch.SA.marine.scale))]=Catch.SA.marine.scale[match(1998:Yr.end,names(Catch.SA.marine.scale))]*(1-PCS.white)
    
  }
  
  if(ADD.PCS=="YES.100")
  {
    TOT.Ktch$AnnualKtch[match(1998:Yr.end,TOT.Ktch$Year)]=TOT.Ktch$AnnualKtch[match(1998:Yr.end,TOT.Ktch$Year)]*(1-1)
    TOT.Ktch$AnnualKtch.GHATF[match(1998:Yr.end,TOT.Ktch$Year)]=TOT.Ktch$AnnualKtch.GHATF[match(1998:Yr.end,TOT.Ktch$Year)]*(1-1)
    Catch.Dropline[match(1998:Yr.end,names(Catch.Dropline))]=Catch.Dropline[match(1998:Yr.end,names(Catch.Dropline))]*(1-1)
    Catch.SA.marine.scale[match(1998:Yr.end,names(Catch.SA.marine.scale))]=Catch.SA.marine.scale[match(1998:Yr.end,names(Catch.SA.marine.scale))]*(1-1)
  }
  
  #9. Total catch
  TOT.Ktch$TOT.CATCH=round(TOT.Ktch$AnnualKtch+TOT.Ktch$AnnualKtch.GHATF+Catch.Albany.whale+Catch.SA.game.fish+
                             Catch.SA.marine.scale+Catch.Dropline)
  
  return(list(TOT.Ktch=TOT.Ktch,
              Catch.Dropline=Catch.Dropline,
              TDGDLF=TOT.Ktch$AnnualKtch,
              GHATF=TOT.Ktch$AnnualKtch.GHATF,
              SA.marine.scale=Catch.SA.marine.scale))
}

# Function for Plotting reconstructed catch
fn.plot.ktch=function(DAT,XLIM)
{
  #get mean and ci
  Mean.Sim.Ktch=colMeans(DAT)
  Low.CI.Sim.Ktch=apply(DAT, 2, function(x)quantile(x,probs=0.025))
  UP.CI.Sim.Ktch=apply(DAT, 2, function(x)quantile(x,probs=0.975))
  
  #construct polygons
  Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
  Biom.Vec <- c(Low.CI.Sim.Ktch, tail(UP.CI.Sim.Ktch, 1), rev(UP.CI.Sim.Ktch), Low.CI.Sim.Ktch[1])
  
  #plot
  plot(YR,Mean.Sim.Ktch,xlim=XLIM,ylim=c(0,YMAX),type="l",lwd=3,ylab="",xlab="",cex.axis=1.25,
       cex.lab=1.75,xaxt='n',yaxt='n')
  polygon(Year.Vec, Biom.Vec, col = "grey80", border = "grey20")
  lines(YR,Mean.Sim.Ktch,col=1,lwd=3)
  axis(1,at=YR,labels=F,tck=-0.015)
  axis(1,at=seq(1940,2010,10),labels=F,tck=-0.03)
  axis(1,at=seq(1940,2010,20),labels=F,tck=-0.045)
  axis(2,at=c(seq(0,YMAX,250)),labels=F,tck=-0.03)
}
fn.plot.ktch.percentile=function(DAT,XLIM)
{
  #get percentiles
  Low.percentile=function(Nper) apply(DAT, 2, function(x) quantile(x, (0+Nper)/100))
  High.percentile=function(Nper) apply(DAT, 2, function(x) quantile(x, (100-Nper)/100))
  
  #50% of data
  Nper=(100-50)/2
  LOW.50=Low.percentile(Nper)
  UP.50=High.percentile(Nper)
  
  #75% of data
  Nper=(100-75)/2
  LOW.75=Low.percentile(Nper)
  UP.75=High.percentile(Nper)
  
  #95% of data
  Nper=(100-95)/2
  LOW.95=Low.percentile(Nper)
  UP.95=High.percentile(Nper)
  
  #100% of data
  Nper=(100-100)/2
  LOW.100=Low.percentile(Nper)
  UP.100=High.percentile(Nper)
  
  
  #construct polygons
  Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
  Biom.Vec.50 <- c(LOW.50, tail(UP.50, 1), rev(UP.50), LOW.50[1])
  Biom.Vec.75 <- c(LOW.75, tail(UP.75, 1), rev(UP.75), LOW.75[1])
  Biom.Vec.95 <- c(LOW.95, tail(UP.95, 1), rev(UP.95), LOW.95[1])
  Biom.Vec.100 <- c(LOW.100, tail(UP.100, 1), rev(UP.100), LOW.100[1])
  
  
  #plot
  plot(YR,UP.95,xlim=XLIM,ylim=c(0,YMAX),type="l",ylab="",xlab="",cex.axis=1.25,
       cex.lab=1.75,xaxt='n',yaxt='n',col='transparent')
  
  
  colfunc <- colorRampPalette(c("grey40","white"))
  COLS=colfunc(4)
  #   polygon(Year.Vec, Biom.Vec.100, col = COLS[4], border = "grey20")
  polygon(Year.Vec, Biom.Vec.95, col = COLS[3], border = "grey20")
  polygon(Year.Vec, Biom.Vec.75, col = COLS[2], border = "grey20")
  polygon(Year.Vec, Biom.Vec.50, col = COLS[1], border = "grey20")
  axis(1,at=YR,labels=F,tck=-0.015)
  axis(1,at=seq(1940,2010,10),labels=F,tck=-0.03)
  axis(1,at=seq(1940,2010,20),labels=F,tck=-0.045)
  axis(2,at=c(seq(0,YMAX,250)),labels=F,tck=-0.03)
}

#Function for incorporating density dependence to risk analysis
fun.dens=function(dens,Matrix,Total.No)
{
  #add survival density-dependence to projection matrix
  if (WhatClass=="ALL")for (nn in 2:nrow(Matrix)) Matrix[nn,(nn-1)]=Matrix[nn,(nn-1)]*exp(-dens*Total.No)
  if (WhatClass=="ONE")for (nn in 2:2) Matrix[nn,(nn-1)]=Matrix[nn,(nn-1)]*exp(-dens*Total.No)
  
  
  #calculate lambda
  LAMBDA=lambda(Matrix)
  
  #obj function
  epsilon=(1-LAMBDA)^2
  if(LAMBDA<1e-2)epsilon=1e5  #penalty to avoid nonsense dens dep factor
  
  return(epsilon)
}

#Function for Risk analysis 
Risk.fun=function(scenario,Tot.Ktch,N.0.samp)
{
  #1. Create elements to fill in
  Pop.size.ratio=N0.boot=No.vec=rep(NA,length = Ns)
  Rel.abundance=Mature.rel.abundance=A3M.rel.abundance=F_vec=matrix(NA,ncol=N.mdld,nrow = Ns)
  MATRIX=vector('list',length=Ns)
  #2. Monte Carlo loop
  for (s in iterations)
  {
    #2.1 Sample catch to add catch uncertainty                                    
    TC=round(Tot.Ktch[s,]*Ktch.sex.ratio)
    N.init=N.0.samp[s]
    
    #2.2. Vary vital rates for projection matrices 
    
    #draw max age sample
    if(scenario==1) A.sim=ASim.1[s]+1 #use same A.sim for all projections to keep same size matrix
    if(scenario==2) A.sim=ASim[s]+1
    
    #draw  age at maturity sample
    if(scenario==1) A.mat=AgeMatSim.1[s] 
    if(scenario==2) A.mat=AgeMatSim[s]
    if(scenario==1) A.3M=round(Age.3MSim.1[s]) 
    if(scenario==2) A.3M=round(Age.3MSim[s])
    
    #take matrix sample of same dimension
    #select matrices of A.sim dimensions
    Proj.Mat=r.numb=NULL
    if(scenario==1) 
    {
      condition=lapply(Proyec.matrix.1, function(x) nrow(x) == A.sim) 
      DATA=Proyec.matrix.1[unlist(condition)]
    }
    
    if(scenario==2) 
    {
      condition=lapply(Proyec.matrix, function(x) nrow(x) == A.sim) 
      DATA=Proyec.matrix[unlist(condition)]
    }
    
    #resample if there are <N.mdld
    if(Proy.type=="Stochas") r.numb=sample(1:length(DATA),N.mdld,replace=T)   #stochastic projection
    if(Proy.type=="Determ") r.numb=rep(1,N.mdld)                            #deterministic projection
    
    #keep only n.Yr.prol proj matrices
    Proj.Mat=DATA[r.numb]
    
    
    #2.2. Calculate selectivity
    Selectivity.sim=NULL
    if(scenario==1) Selectivity.sim=SelSim.1[s,]
    if(scenario==2) Selectivity.sim=SelSim[s,]
    
    Selectivity.sim=subset(Selectivity.sim,!is.na(Selectivity.sim)) #remove NAs
    Selectivity.sim=c(Selectivity.sim,Selectivity.sim[length(Selectivity.sim)])
    
    #2.3. Add harvesting                   
    harvest.matrix=function(matrix) 
    {
      
      H=diag(nrow(matrix))
      diag(H)=1-(U*Selectivity.sim) #apply U and selectivity
      MH=matrix%*%H
      return(MH)
    }
    Harvest.Proyec.mat=vector("list",length = n.Yr.tot)
    
    
    #pre-exploitation population size
    if(Proy.type=="Determ") Matrix=Proj.Mat[[1]] #previous
    if(Proy.type=="Stochas")
    {
      ind=unlist(lapply(Proj.Mat,lambda))
      ind1=sort(ind)
      ind2=round(length(ind1)/2)
      ind3=which(ind==ind1[ind2])[1]            #use mode lambda to calculate dens factor    
      Matrix=Proj.Mat[[ind3]]
      
      ind4=which(ind==ind1[ind2+1])[1]      
      ind5=which(ind==ind1[ind2-1])[1]      
      Mat.dummy=list(Proj.Mat[[ind3]],Proj.Mat[[ind4]],Proj.Mat[[ind5]])  
      n.dummy=sample(1:length(Mat.dummy),n.Yr.tot,replace=T)   
      Proj.Mat=Mat.dummy[n.dummy]
    }
    
    No=matrix(stable.stage(Matrix)) 
    Total.No=sum(No)
    Mature.No=sum(No[A.mat:A.sim])
    A3M.No=sum(No[A.3M:A.sim])
    
    #2.4. Add density-depence in survival
    if(dens.dep=="YES")   
    {
      fn_obj=function(dens)fun.dens(dens,Matrix,Total.No) #objfun to MAXIMIZE       
      if (WhatClass=="ALL")fit.dens=optimize(fn_obj,lower=1e-5,upper=1e-1)
      if (WhatClass=="ONE")fit.dens=optimize(fn_obj,lower=1e-2,upper=100)
      Dens.factor=fit.dens$minimum
    }
    
    #Burn-in for stabilising population
    n.vec=vector("list",length =N.mdld)
    n.vec[[1]]=No
    
    for (p in 2:n.burn.in)
    {
      #Density-depence in survival
      dummy=Matrix                     
      if (WhatClass=="ALL")for (nn in 2:nrow(dummy))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
      if (WhatClass=="ONE")for (nn in 2:2)dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
      n.vec[[p]]=dummy%*%matrix(n.vec[[p-1]])          
    }
    
    
    #2.5. Population trajectories with exploitation
    #exploited population size
    Fvec=rep(NA,n.Yr.tot)
    for (p in 1:n.Yr.tot)
    {
      pp=p+n.burn.in  #index for population that considers burn in
      
      Project.M=Proj.Mat[[p]]
      
      #Density-depence in survival
      if(dens.dep=="YES")
      {
        if(TC[p]<sum(c(n.vec[[pp-1]]))*N.init)
        {
          U=TC[p]/(sum(c(n.vec[[pp-1]]))*N.init)      
          U=-log(1-U)             #convert annual exploitation rate to instaneous F 
        }
        if(TC[p]>=sum(c(n.vec[[pp-1]]))*N.init) U=-log(1-0.9999999999)
        
        if(TC[p]>0)Harvest.Proyec.mat[[p]]=harvest.matrix(Project.M)  
        if(TC[p]==0)Harvest.Proyec.mat[[p]]=Project.M   
        
        dummy=Harvest.Proyec.mat[[p]]          
        if (WhatClass=="ALL")for (nn in 2:nrow(dummy))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[pp-1]]))
        if (WhatClass=="ONE")for (nn in 2:2)dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[pp-1]]))
        
        n.vec[[pp]]=dummy%*%matrix(n.vec[[pp-1]])          
        if(sum(n.vec[[pp]])<0) n.vec[[pp]]=rep(0,length(n.vec[[pp-1]]))
        Fvec[p]=U
      }
    }
    
    Pop.size=Mature.Pop.size=A3M.Pop.size=rep(0,(N.mdld))
    for(y in 1:(N.mdld))
    {
      Pop.size[y]=sum(n.vec[[y]])*N.init
      Mature.Pop.size[y]=sum(n.vec[[y]][A.mat:A.sim])*N.init
      A3M.Pop.size[y]=sum(n.vec[[y]][A.3M:A.sim])*N.init
    }
    
    
    #2.6. Store outputs
    N0.boot[s]=N.init
    Rel.abundance[s,]=Pop.size
    Mature.rel.abundance[s,]=Mature.Pop.size
    A3M.rel.abundance[s,]=A3M.Pop.size
    No.vec[s]=Mature.No
    F_vec[s,]=c(rep(0,n.burn.in),Fvec)
    MATRIX[[s]]=Matrix
  }
  
  return(list(Rel.abundance=Rel.abundance,N0=N0.boot,Mature.rel.abundance=Mature.rel.abundance,
              No.vec=No.vec,F_series=F_vec,Matrix=MATRIX,A3M.rel.abundance=A3M.rel.abundance))
}

#Function for plotting distributions
BarPlot=function(Data,DataName,Names.Args)
{
  a=table(Data)/sum(table(Data))
  barplot(a,ylim=c(0,max(a)*1.1),ylab="",xlab=DataName,names.arg=Names.Args,
          axes=T,col=Basic,axis.lty=1,cex.names=Ax.cx,cex.lab=Lab.xc,yaxt="n")
  box()
}
fn.barplo=function(data1,data2,range1,range2,TITLE)
{
  table.1=table(data1)/sum(table(data1))
  table.=table(data2)/sum(table(data2))
  
  largo=length(table.)
  nm.1=as.numeric(names(table.1))
  nm.=as.numeric(names(table.))
  
  min.=nm.[1]-nm.1[1]  
  max.=nm.[length(nm.)]-nm.1[length(nm.1)]
  
  table.1[(length(table.1)+1):(length(table.1)+max.)]=0
  if(!min.==0)
  {
    table.[(largo+1):(largo+1+min.)]=0
    table.=table.[c((largo+1):(largo+min.),1:largo)]
  }
  
  combined=rbind(table.1,table.)
  
  barplot(combined, beside = TRUE,ylim=c(0,max(combined)*1.1),mgp = c(2.5, 0.6, 0),
          names.arg= range1[1]:range2[2],xlab=TITLE,yaxt="n",
          ylab="", axis.lty=1, axes=T,col=greyscale,cex.names=Ax.cx,las=1,cex.lab=Lab.xc,space=c(.1,.1))
  box()
  
  
}

#Function for plotting population trajectories
plot.mean.CI.Report=function(LH1,LH2,CEXX,Fem_3M,Whr.leg) 
{
  DATA=LH1[,These.Yrs]
  DATA.LH2=LH2[,These.Yrs]
  YRANGE=c(0,1.05*max(c(DATA[,1],DATA.LH2[,2])))
  if(Fem_3M=="YES") YRANGE=c(0,YLM3M)
  Low.percentile=function(Nper,DAT) apply(DAT, 2, function(x) quantile(x, (0+Nper)/100))
  High.percentile=function(Nper,DAT) apply(DAT, 2, function(x) quantile(x, (100-Nper)/100))
  
  plot(1:n.Yr.prol,type='l',lty=1,col="transparent",lwd=2,ylim=YRANGE,xaxt='n',yaxt='n',xlab="",ylab="")
  
  #LH1
  #Median
  Median=apply(DATA, 2, median)
  
  #95% of data
  Nper=(100-95)/2
  LOW.95=Low.percentile(Nper,DATA)
  UP.95=High.percentile(Nper,DATA)
  YR=1:n.Yr.prol
  Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
  Biom.Vec.95 <- c(LOW.95, tail(UP.95, 1), rev(UP.95), LOW.95[1])
  
  #LH2
  #Median
  Median.LH2=apply(DATA.LH2, 2, median)
  
  #95% of data
  Nper=(100-95)/2
  LOW.95.LH2=Low.percentile(Nper,DATA.LH2)
  UP.95.LH2=High.percentile(Nper,DATA.LH2)
  Biom.Vec.95.LH2 <- c(LOW.95.LH2, tail(UP.95.LH2, 1), rev(UP.95.LH2), LOW.95.LH2[1])
  
  COLS=colfunc(2)
  polygon(Year.Vec, Biom.Vec.95, col = COLS[2], border = Brdr)
  polygon(Year.Vec, Biom.Vec.95.LH2, col = COLS[1], border = Brdr)
  lines(Median,col=LIN,lwd=2)
  lines(Median.LH2,col=LIN,lwd=2,lty=2)
  
  axis(1,at=seq(1,n.Yr.prol,by=1),labels=F,las=1,tck=-0.02)
  axis(1,at=seq(match(1940,Yr.proj),n.Yr.prol,by=10),labels=F,las=1,tck=-0.04)
  axis(1,at=seq(match(1940,Yr.proj),n.Yr.prol,by=20),labels=F,tck=-0.06)
  if(z==5)axis(1,at=seq(match(1940,Yr.proj),n.Yr.prol,by=20),labels=c("1940","1960","1980","2000"),tck=-0.06,cex.axis=CEXX)
  #  if(z==5)axis(1,at=seq(match(1940,Yr.proj),n.Yr.prol,by=20),labels=c("1940-41","1960-61","1980-81","2000-01"),tck=-0.06,cex.axis=CEXX)
  
  YAxis=seq(0,N.0.Rick[z],length.out=5)
  axis(2,at=YAxis,labels=F,las=1,tck=-0.06)
  if(y==1)axis(2,at=YAxis,labels=YAxis,las=1,tck=-0.06,cex.axis=CEXX)
  if(z==1 & y==1) legend(Whr.leg,c("LH1","LH2"),fill=rev(COLS),bty='n',cex=1.75)
  box()
  return(list(Median.LH1=Median,Median.LH2=Median.LH2))
}
plot.mean.CI.Paper=function(LH1,LH2,CEXX) 
{
  DATA=LH1[,These.Yrs]
  DATA=DATA/DATA[,1]
  DATA.LH2=LH2[,These.Yrs]
  DATA.LH2=DATA.LH2/DATA.LH2[,1]
  YRANGE=c(0,1.05*max(c(DATA[,1],DATA.LH2[,2])))
  Low.percentile=function(Nper,DAT) apply(DAT, 2, function(x) quantile(x, (0+Nper)/100))
  High.percentile=function(Nper,DAT) apply(DAT, 2, function(x) quantile(x, (100-Nper)/100))
  
  plot(1:n.Yr.prol,type='l',lty=1,col="transparent",lwd=2,ylim=YRANGE,xaxt='n',yaxt='n',xlab="",ylab="")
  
  #LH1
  #Median
  Median=apply(DATA, 2, median)
  
  #95% of data
  Nper=(100-95)/2
  LOW.95=Low.percentile(Nper,DATA)
  UP.95=High.percentile(Nper,DATA)
  YR=1:n.Yr.prol
  Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
  Biom.Vec.95 <- c(LOW.95, tail(UP.95, 1), rev(UP.95), LOW.95[1])
  
  #LH2
  #Median
  Median.LH2=apply(DATA.LH2, 2, median)
  
  #95% of data
  Nper=(100-95)/2
  LOW.95.LH2=Low.percentile(Nper,DATA.LH2)
  UP.95.LH2=High.percentile(Nper,DATA.LH2)
  Biom.Vec.95.LH2 <- c(LOW.95.LH2, tail(UP.95.LH2, 1), rev(UP.95.LH2), LOW.95.LH2[1])
  
  #   #50% of data
  #   Nper=(100-50)/2
  #   LOW.50=Low.percentile(Nper)
  #   UP.50=High.percentile(Nper)
  #   Biom.Vec.50 <- c(LOW.50, tail(UP.50, 1), rev(UP.50), LOW.50[1])
  #   
  #   #75% of data
  #   Nper=(100-75)/2
  #   LOW.75=Low.percentile(Nper)
  #   UP.75=High.percentile(Nper)
  #   Biom.Vec.75 <- c(LOW.75, tail(UP.75, 1), rev(UP.75), LOW.75[1])
  #   
  #   #100% of data
  #   Nper=(100-100)/2
  #   LOW.100=Low.percentile(Nper)
  #   UP.100=High.percentile(Nper)
  
  COLS=colfunc(2)
  polygon(Year.Vec, Biom.Vec.95, col = COLS[2], border = Brdr)
  polygon(Year.Vec, Biom.Vec.95.LH2, col = COLS[1], border = Brdr)
  lines(Median,col=LIN,lwd=2)
  lines(Median.LH2,col=LIN,lwd=2,lty=2)
  
  axis(1,at=seq(1,n.Yr.prol,by=1),labels=F,las=1,tck=-0.02)
  axis(1,at=seq(match(1940,Yr.proj),n.Yr.prol,by=10),labels=F,las=1,tck=-0.04)
  axis(1,at=seq(match(1940,Yr.proj),n.Yr.prol,by=20),labels=F,tck=-0.06)
  if(z==5)
  {
    Nnn=seq(match(1940,Yr.proj),n.Yr.prol,by=20)
    axis(1,at=Nnn,labels=Yr.proj[Nnn],tck=-0.06,cex.axis=CEXX)    
  }
  
  YAxis=seq(0,1,length.out=5)
  axis(2,at=YAxis,labels=F,las=1,tck=-0.06)
  if(y==1)axis(2,at=YAxis,labels=YAxis,las=1,tck=-0.06,cex.axis=CEXX)
  if(z==5 & y==1) legend("bottomleft",c("LH1","LH2"),fill=rev(COLS),bty='n',cex=1.75)
  
  box()
  return(list(Median.LH1=Median,Median.LH2=Median.LH2))
}
fn.yr.axis=function(CEx)
{
  axis(1,at=seq(match(1940,Yr.proj),n.Yr.prol,by=20),
       #labels=seq(Yr.proj[match(1940,Yr.proj)],Yr.proj[length(Yr.proj)],by=20),
       labels=c("1940/41","1960/61","1980/81","2000/01"),
       cex.axis=CEx,las=1,tck=-0.06)
}

#Function for getting depletion and increase rate since 1998
fn.inc.rate.depl=function(LH1,LH2)
{
  DATA=LH1[,These.Yrs]
  DATA.LH2=LH2[,These.Yrs]
  
  #Depletion level
  Depletn=DATA[,ncol(DATA)]/DATA[,1]
  Depletn.LH2=DATA.LH2[,ncol(DATA.LH2)]/DATA.LH2[,1]
  
  #Increase since 1998
  Rat.inc=DATA[,ncol(DATA)]/LH1[,id.1998]
  Rat.inc[is.na(Rat.inc)]=0
  Rat.inc=100*(Rat.inc-1)
  Rat.inc=ifelse(Rat.inc<0,0,Rat.inc)
  
  Rat.inc.LH2=DATA.LH2[,ncol(DATA.LH2)]/LH2[,id.1998]
  Rat.inc.LH2[is.na(Rat.inc.LH2)]=0
  Rat.inc.LH2=100*(Rat.inc.LH2-1)
  Rat.inc.LH2=ifelse(Rat.inc.LH2<0,0,Rat.inc.LH2)
  
  
  return=list(Depletion=c(LH1=median(Depletn),LH2=median(Depletn.LH2)),
              Rate.inc=c(LH1=median(Rat.inc),LH2=median(Rat.inc.LH2)))
}

#Function for Steve's graphs
fn.barplot=function(DAT,NAMES)
{
  #Mean
  a=reshape(DAT[,c(1:2,4)],v.names = "Mean.cpue.folly", timevar = "Time", idvar= "Fishing.region", direction = "wide")
  
  #Up
  b=reshape(DAT[,c(1:2,7)],v.names = "Up95.cpue.folly", timevar = "Time",idvar = "Fishing.region", direction = "wide")
  
  #low
  d=reshape(DAT[,c(1:2,6)],v.names = "Low95.cpue.folly", timevar = "Time",idvar = "Fishing.region", direction = "wide")
  
  if(NAMES=="YES")
  {
    barplot2(as.matrix(a[c(3,1:2),2:4]), plot.ci=TRUE, ci.l=as.matrix(d[c(3,1:2),2:4]), 
             ci.u=as.matrix(b[c(3,1:2),2:4]),beside=T,ylim=c(0,max(b[,2:4])*1.01), xpd=FALSE,main="",
             names.arg=c("1988/89-1996/97","1997/98-2004/05","2005/06-2014/15"),cex.axis=1.25,
             ylab="",xlab="",col=c("grey80","white","grey50"),cex.names=1.25,axis.lty=1)
  }else
  {
    barplot2(as.matrix(a[c(3,1:2),2:4]), plot.ci=TRUE, ci.l=as.matrix(d[c(3,1:2),2:4]), 
             ci.u=as.matrix(b[c(3,1:2),2:4]),beside=T,ylim=c(0,max(b[,2:4])*1.01), xpd=FALSE,main="",
             names.arg=c("","",""),ylab="",xlab="",col=c("grey80","white","grey50"),cex.axis=1.25,axis.lty=1)
    
  }
  
  box()
}


#Functions for sensitivity analyses
fn.get.mean.min.max=function(d) return(data.frame(Mean=mean(d,na.rm=T),Min=min(d,na.rm=T),Max=max(d,na.rm=T)))
fun.sensitivity=function(SenSc,N.init)
{
  #1. get values of relevant quantities
  
  # Catch
  Full.ktch=do.call(rbind,TOTAL.CATCH.PCS$NO)
  if(SenSc$Par=="Catch")
  {
    if(SenSc$value=="Min") Tot.Ktch=apply(Full.ktch,2,min) 
    if(SenSc$value=="Max") Tot.Ktch=apply(Full.ktch,2,max)      
  } else  Tot.Ktch=colMeans(Full.ktch)
  
  # PCM
  if(SenSc$Par=="PCM")
  {
    if(SenSc$value=="0") PCM_s=0
    if(SenSc$value=="100") PCM_s=1     
  } else  PCM_s=0.5
  
  #Apply PCM to Catch
  Tot.Ktch[match(1998:Yr.end,1938:Yr.end)]=Tot.Ktch[match(1998:Yr.end,1938:Yr.end)]*PCM_s
  
  
  # Selectivity
  if(SenSc$Par=="Selectivity")
  {
    if(SenSc$value=="Min") 
    {
      SEL1=apply(SelSim.1,2,min,na.rm=T)
      SEL2=apply(SelSim,2,min,na.rm=T)
      SEL1[(length(SEL1)+1):length(SEL2)]=mean(SEL1[(length(SEL1)-10):length(SEL1)])
      Selectivity.sim=apply(rbind(SEL1,SEL2),2,min,na.rm=T)        
    }
    if(SenSc$value=="Max") 
    {
      SEL1=apply(SelSim.1,2,max,na.rm=T)
      SEL2=apply(SelSim,2,max,na.rm=T)
      SEL1[(length(SEL1)+1):length(SEL2)]=mean(SEL1[(length(SEL1)-10):length(SEL1)])
      Selectivity.sim=apply(rbind(SEL1,SEL2),2,max,na.rm=T)        
    }       
  } else  
  {
    SEL1=apply(SelSim.1,2,mean,na.rm=T)
    SEL2=apply(SelSim,2,mean,na.rm=T)
    SEL1[(length(SEL1)+1):length(SEL2)]=mean(SEL1[(length(SEL1)-10):length(SEL1)])
    Selectivity.sim=apply(rbind(SEL1,SEL2),2,mean,na.rm=T)        
  }
  
  # Max.Age
  A.SIM=c(ASim,ASim.1)
  if(SenSc$Par=="Max.Age")
  {
    if(SenSc$value=="Min")  A.sim=min(A.SIM,na.rm=T)
    if(SenSc$value=="Max")  A.sim=max(A.SIM,na.rm=T)-1  
  } else  A.sim=round(mean(A.SIM,na.rm=T))
  
  #age classes
  age=first.age:A.sim
  
  
  # Mat.Age
  A.mat.SIM=c(AgeMatSim.1,AgeMatSim)
  if(SenSc$Par=="Mat.Age")
  {
    if(SenSc$value=="Min") 
    {
      age.mat.sim=min(A.mat.SIM,na.rm=T)  #reset max age due to correlation age mat and max age
      A.sim=min(A.SIM,na.rm=T)
    }
    
    if(SenSc$value=="Max")
    {
      age.mat.sim=max(A.mat.SIM,na.rm=T) 
      A.sim=max(A.SIM,na.rm=T)-1  
    }
    
  } else  age.mat.sim=round(mean(A.mat.SIM,na.rm=T))
  
  # Rep.cycl
  if(SenSc$Par=="Rep.cycl")
  {
    if(SenSc$value=="Min")  Reprod_cycle.sim=min(Range.reprod_cycle)
    if(SenSc$value=="Max")  Reprod_cycle.sim=max(Range.reprod_cycle)          
  } else  Reprod_cycle.sim=mean(Range.reprod_cycle)
  Pmax.sim=1/Reprod_cycle.sim
  
  # Pups
  if(SenSc$Par=="Pups")
  {
    if(SenSc$value=="Min")  Meanfec.sim=min(Rangefec)
    if(SenSc$value=="Max")  Meanfec.sim=max(Rangefec)        
  } else  Meanfec.sim=round(mean(Rangefec))
  
  # M
  if(SenSc$Par=="M")
  {
    if(SenSc$value=="Min")
    {
      M1=apply(mSim.1,2,min,na.rm=T)
      M2=apply(mSim,2,min,na.rm=T)
      M1[(length(M1)+1):length(M2)]=mean(M1[(length(M1)-10):length(M1)])
      m.sim=apply(rbind(M1,M2),2,min,na.rm=T)
    }
    
    if(SenSc$value=="Max")
    {
      M1=apply(mSim.1,2,max,na.rm=T)
      M2=apply(mSim,2,max,na.rm=T)
      M1[(length(M1)+1):length(M2)]=mean(M1[(length(M1)-10):length(M1)])
      m.sim=apply(rbind(M1,M2),2,max,na.rm=T)
    }
    
  } else
  {
    M1=apply(mSim.1,2,mean,na.rm=T)
    M2=apply(mSim,2,mean,na.rm=T)
    M1[(length(M1)+1):length(M2)]=mean(M1[(length(M1)-10):length(M1)])
    m.sim=apply(rbind(M1,M2),2,mean,na.rm=T)     
  }
  
  #Recalculate m based on life history pars tested and mean Lo, K, Linf for max Age.max to allow equilibrium condition with no fishing
  #   total.length=Lo+(Linf*100-Lo)*(1-exp(-k*age))
  #   fork=(0.9442*total.length*100)-5.7441     
  #   mid.FL.fem=fork
  #   m.sim=M.fun(A=A.sim,k=k,Linf=Linf,Aver.T=Aver.T,age.mat=age.mat.sim,bwt,awt,W.G=1,W.noG=1)    
  
  
  
  #2. get projection matrix
  
  #survivorship
  S=exp(-m.sim[1:length(age)])         
  
  #proportion surviving
  lx=rep(0,length(age))
  lx[1]=1.0
  for (i in 2:(length(age)))lx[i]=lx[i-1]*S[i]
  
  #reproductive schedules   
  MF=c(rep(0,(age.mat.sim)),rep(Meanfec.sim,length(age)-age.mat.sim))
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
  DATA=matrix(0,nrow=n,ncol=n)
  diag(DATA)[-nrow(DATA)]=PX
  DATA=rbind(matrix(0,nrow=1,ncol=n),DATA)
  DATA=DATA[-(n+1),]
  DATA[1,]=BX
  rownames(DATA)=colnames(DATA)=(first.age+1):n
  
  
  
  #3. get population trajectories
  
  #Sample catch to add catch uncertainty                                    
  TC=round(Tot.Ktch*Ktch.sex.ratio)
  A.mat=age.mat.sim
  Proj.Mat=DATA
  Selectivity.sim=Selectivity.sim[1:length(age)]
  
  #Add harvesting                   
  harvest.matrix=function(matrix) 
  {
    
    H=diag(nrow(matrix))
    diag(H)=1-(U*Selectivity.sim) #apply U and selectivity
    MH=matrix%*%H
    return(MH)
  }
  Harvest.Proyec.mat=vector("list",length = n.Yr.tot)
  Matrix=Proj.Mat  
  No=matrix(stable.stage(Matrix)) 
  Total.No=sum(No)
  
  #Add density-depence in survival
  if(dens.dep=="YES")   
  {
    fn_obj=function(dens)fun.dens(dens,Matrix,Total.No) #objfun to MAXIMIZE       
    if (WhatClass=="ALL")fit.dens=optimize(fn_obj,lower=1e-5,upper=1e-1)
    if (WhatClass=="ONE")fit.dens=optimize(fn_obj,lower=1e-2,upper=100)
    Dens.factor=fit.dens$minimum
  }
  
  #Burn-in for stabilising population
  n.vec=vector("list",length =N.mdld)
  n.vec[[1]]=No
  
  for (p in 2:n.burn.in)
  {
    #Density-depence in survival
    dummy=Matrix                     
    if (WhatClass=="ALL")for (nn in 2:nrow(dummy))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
    if (WhatClass=="ONE")for (nn in 2:2)dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
    n.vec[[p]]=dummy%*%matrix(n.vec[[p-1]])          
  }
  
  
  #Population trajectories with exploitation
  #exploited population size
  Fvec=rep(NA,n.Yr.tot)
  for (p in 1:n.Yr.tot)
  {
    pp=p+n.burn.in  #index for population that considers burn in
    
    Project.M=Proj.Mat
    
    #Density-depence in survival
    if(dens.dep=="YES")
    {
      if(TC[p]<sum(c(n.vec[[pp-1]]))*N.init)
      {
        U=TC[p]/(sum(c(n.vec[[pp-1]]))*N.init)      
        U=-log(1-U)             #convert annual exploitation rate to instaneous F 
      }
      if(TC[p]>=sum(c(n.vec[[pp-1]]))*N.init) U=-log(1-0.9999999999)
      
      if(TC[p]>0)Harvest.Proyec.mat[[p]]=harvest.matrix(Project.M)  
      if(TC[p]==0)Harvest.Proyec.mat[[p]]=Project.M   
      
      dummy=Harvest.Proyec.mat[[p]]          
      if (WhatClass=="ALL")for (nn in 2:nrow(dummy))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[pp-1]]))
      if (WhatClass=="ONE")for (nn in 2:2)dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[pp-1]]))
      
      n.vec[[pp]]=dummy%*%matrix(n.vec[[pp-1]])          
      if(sum(n.vec[[pp]])<0) n.vec[[pp]]=rep(0,length(n.vec[[pp-1]]))
      Fvec[p]=U
    }
  }
  
  Pop.size=rep(0,(N.mdld))
  for(y in 1:(N.mdld))Pop.size[y]=sum(n.vec[[y]])*N.init
  
  rm(Tot.Ktch,PCM_s,Selectivity.sim,A.sim,age.mat.sim,Pmax.sim,Meanfec.sim,m.sim)
  
  return(list(Rel.abundance=Pop.size))  
}
fun.plot.Sens.Ref=function(REF,CE)
{
  plot(1:length(These.Yrs),REF$Rel.abundance[These.Yrs],xaxt='n',ylab="",cex.axis=CE,
       xlab="",ylim=c(0,max(REF$Rel.abundance[These.Yrs])),type='l',lwd=2)
  axis(1,1:length(These.Yrs),F,tck=-0.03)  
  axis(1,at=AT,labels=F,tck=-0.06,cex.axis=CE)  
}
fun.plot.Sens=function(REF,CL,LTY) lines(1:length(These.Yrs),REF$Rel.abundance[These.Yrs],lwd=3,col=CL,lty=LTY)
fun.sensitivity.get.input.values=function(SenSc,N.init)
{
  #1. get values of relevant quantities
  
  # Catch
  Full.ktch=do.call(rbind,TOTAL.CATCH.PCS$NO)
  if(SenSc$Par=="Catch")
  {
    if(SenSc$value=="Min") Tot.Ktch=apply(Full.ktch,2,min) 
    if(SenSc$value=="Max") Tot.Ktch=apply(Full.ktch,2,max)      
  } else  Tot.Ktch=colMeans(Full.ktch)
  
  # PCM
  if(SenSc$Par=="PCM")
  {
    if(SenSc$value=="0") PCM_s=0
    if(SenSc$value=="100") PCM_s=1     
  } else  PCM_s=0.5
  
  #Apply PCM to Catch
  Tot.Ktch[match(1998:Yr.end,1938:Yr.end)]=Tot.Ktch[match(1998:Yr.end,1938:Yr.end)]*PCM_s
  Tot.Ktch=round(Tot.Ktch*Ktch.sex.ratio)
  
  # Selectivity
  if(SenSc$Par=="Selectivity")
  {
    if(SenSc$value=="Min") 
    {
      SEL1=apply(SelSim.1,2,min,na.rm=T)
      SEL2=apply(SelSim,2,min,na.rm=T)
      SEL1[(length(SEL1)+1):length(SEL2)]=mean(SEL1[(length(SEL1)-10):length(SEL1)])
      Selectivity.sim=apply(rbind(SEL1,SEL2),2,min,na.rm=T)        
    }
    if(SenSc$value=="Max") 
    {
      SEL1=apply(SelSim.1,2,max,na.rm=T)
      SEL2=apply(SelSim,2,max,na.rm=T)
      SEL1[(length(SEL1)+1):length(SEL2)]=mean(SEL1[(length(SEL1)-10):length(SEL1)])
      Selectivity.sim=apply(rbind(SEL1,SEL2),2,max,na.rm=T)        
    }       
  } else  
  {
    SEL1=apply(SelSim.1,2,mean,na.rm=T)
    SEL2=apply(SelSim,2,mean,na.rm=T)
    SEL1[(length(SEL1)+1):length(SEL2)]=mean(SEL1[(length(SEL1)-10):length(SEL1)])
    Selectivity.sim=apply(rbind(SEL1,SEL2),2,mean,na.rm=T)        
  }
  
  # Max.Age
  A.SIM=c(ASim,ASim.1)
  if(SenSc$Par=="Max.Age")
  {
    if(SenSc$value=="Min")  A.sim=min(A.SIM,na.rm=T)
    if(SenSc$value=="Max")  A.sim=max(A.SIM,na.rm=T)-1  
  } else  A.sim=round(mean(A.SIM,na.rm=T))
  
  #age classes
  age=first.age:A.sim
  
  
  # Mat.Age
  A.mat.SIM=c(AgeMatSim.1,AgeMatSim)
  if(SenSc$Par=="Mat.Age")
  {
    if(SenSc$value=="Min") 
    {
      age.mat.sim=min(A.mat.SIM,na.rm=T)  #reset max age due to correlation age mat and max age
      A.sim=min(A.SIM,na.rm=T)
    }
    
    if(SenSc$value=="Max")
    {
      age.mat.sim=max(A.mat.SIM,na.rm=T) 
      A.sim=max(A.SIM,na.rm=T)-1  
    }
    
  } else  age.mat.sim=round(mean(A.mat.SIM,na.rm=T))
  
  # Rep.cycl
  if(SenSc$Par=="Rep.cycl")
  {
    if(SenSc$value=="Min")  Reprod_cycle.sim=min(Range.reprod_cycle)
    if(SenSc$value=="Max")  Reprod_cycle.sim=max(Range.reprod_cycle)          
  } else  Reprod_cycle.sim=mean(Range.reprod_cycle)
  Pmax.sim=1/Reprod_cycle.sim
  
  # Pups
  if(SenSc$Par=="Pups")
  {
    if(SenSc$value=="Min")  Meanfec.sim=min(Rangefec)
    if(SenSc$value=="Max")  Meanfec.sim=max(Rangefec)        
  } else  Meanfec.sim=round(mean(Rangefec))
  
  # M
  if(SenSc$Par=="M")
  {
    if(SenSc$value=="Min")
    {
      M1=apply(mSim.1,2,min,na.rm=T)
      M2=apply(mSim,2,min,na.rm=T)
      M1[(length(M1)+1):length(M2)]=mean(M1[(length(M1)-10):length(M1)])
      m.sim=apply(rbind(M1,M2),2,min,na.rm=T)
    }
    
    if(SenSc$value=="Max")
    {
      M1=apply(mSim.1,2,max,na.rm=T)
      M2=apply(mSim,2,max,na.rm=T)
      M1[(length(M1)+1):length(M2)]=mean(M1[(length(M1)-10):length(M1)])
      m.sim=apply(rbind(M1,M2),2,max,na.rm=T)
    }
    
  } else
  {
    M1=apply(mSim.1,2,mean,na.rm=T)
    M2=apply(mSim,2,mean,na.rm=T)
    M1[(length(M1)+1):length(M2)]=mean(M1[(length(M1)-10):length(M1)])
    m.sim=apply(rbind(M1,M2),2,mean,na.rm=T)     
  }
  
  return(list(Tot.Ktch=Tot.Ktch,PCM=PCM_s,Selectivity=Selectivity.sim,
              Age_max=A.sim,Age_mat=age.mat.sim,Cycle=Reprod_cycle.sim,Fec=Meanfec.sim,M=m.sim))  
}
show.sens.par=function(Y,a,b,d)
{
  YY=c(a,b,d)
  YLIM=c(0,max(YY))
  plot(1:length(Y),a,xaxt='n',ylab="",type='l',xlab="",cex.axis=1.5,
       lty=LTY.ppr[3],col=CL.ppr[3],lwd=3,ylim=YLIM)
  axis(1,1:length(Y),F,tck=-0.03)  
  lines(1:length(Y),b,lty=LTY.ppr[1],col=CL.ppr[1],lwd=3)
  lines(1:length(Y),d,lty=LTY.ppr[2],col=CL.ppr[2],lwd=3)  
}



#---5.  LESLIE MATRICES SECTION---  


# Calculate maturity ogive
      #logistic glm approach
glm1 <-glm(cbind(n.mat,n.imm)~mid.point,data=Fem.maturity,family=binomial)
pred.glm=predict(glm1,type ="response")

      #get L50 and L95
lrPerc <-function(cf,p) (log(p/(1-p))-cf[1])/cf[2]
L50 <-lrPerc(coef(glm1),0.5)
L95 <-lrPerc(coef(glm1),0.95)
plot(Fem.maturity$mid.point,Fem.maturity$prop.mat)
lines(Fem.maturity$mid.point,pred.glm,col=3)



#---Monte Carlo simulations of Leslie Matrix----

Ns <- 12500	  #number of simulations
iterations=1:Ns

type=paste("iter.",iterations,sep="")
output=vector('list',length=length(Life.hist.scenarios))

   #1.1. run simulations 
system.time(for(l in 1:length(Life.hist.scenarios))      #takes 0.03 secs per iteration
 {    
    #1. Set scenarios
    #biological scenario
    scenario=Life.hist.scenarios[[l]]
  
    
    #2. Draw random samples of parameters and calculate population quantities
    A.SIM=k.SIM=Linf.SIM=Reprod_cycle.SIM=Pmax.SIM=age.mat.SIM=Age.3M.SIM=Aver.T.SIM=vector(length = Ns)
    Compare.M=Meanfec.SIM=m.SIM=px.SIM=bx.SIM=Proyec.matrix=vector("list",length = Ns)
    
    for (s in iterations)
    {
  
      #1. Draw parameter samples and build functions
      #maximum age
      if(scenario==1)A.sim=floor(rtriangle(1, a=RangeA.1[1], b=RangeA.1[2]+1, c=RangeA.1[1]))    #triangular dist
      if(scenario==2)A.sim=floor(rtriangle(1, a=RangeA.2[1], b=RangeA.2[2]+1, c=Longevity.Hamady))
      A.SIM[s]=A.sim
      
      #growth
      if(scenario==1) 
      {    
        growth.pars.sim=rmvnorm(1,mean=c(k.mean,Linf.mean),sigma=sigma)    #multivariate normal dist
        if(growth.pars.sim[,2]<.06 | growth.pars.sim[,2]>.075 | growth.pars.sim[,1]<=0)    #repeat until sensible pars obtained
        { repeat 
        {
          growth.pars.sim=rmvnorm(1,mean=c(k.mean,Linf.mean),sigma=sigma)
          if(growth.pars.sim[,2]>.06 & growth.pars.sim[,2]<.075)break
        }
        }
        k.sim=growth.pars.sim[1]
        Linf.sim=growth.pars.sim[2]*100
        Lo.sim=Lo
      }
      
      if(scenario==2) 
      {    
        growth.pars.sim=rmvnorm(1,mean=c(k.mean.2,Linf.mean.2),sigma=sigma)    #multivariate normal dist
        if(growth.pars.sim[,2]<.06 | growth.pars.sim[,2]>.075 | growth.pars.sim[,1]<=0)    #repeat until sensible pars obtained
        { repeat 
        {
          growth.pars.sim=rmvnorm(1,mean=c(k.mean.2,Linf.mean.2),sigma=sigma)
          if(growth.pars.sim[,2]>.06 & growth.pars.sim[,2]<.075)break
        }
        }
        k.sim=growth.pars.sim[1]
        Linf.sim=growth.pars.sim[2]*100
        Lo.sim=Lo.2
      }
      k.SIM[s]=k.sim;Linf.SIM[s]=Linf.sim;

      #number of age classes
      age=first.age:A.sim
      
      #reproductive cycle
      Reprod_cycle.sim=sample(Range.reprod_cycle,1,replace=T,prob = prob.rep.cycle) #triangular dist
      Reprod_cycle.SIM[s]=Reprod_cycle.sim
      Pmax.sim=1/Reprod_cycle.sim
      Pmax.SIM[s]=Pmax.sim
      
      #fecundity and age at maturity
      if(scenario==1) 
      {
        Meanfec.sim=rep(ceiling(rtriangle(1,a=Rangefec[1]-1,b=Rangefec[length(Rangefec)],
                                        c=Rangefec[length(Rangefec)])),length(age)) #triangular dist
        age.mat.sim=round(rtriangle(1,a=age.mat[1],b=age.mat[2],c=floor(mean(age.mat))))    #triangular dist
      }
      if(scenario==2) 
      {
        Meanfec.sim=rep(ceiling(runif(1,min=Rangefec2[1]-1,max=Rangefec2[length(Rangefec2)])),length(age)) #uniform dist
        age.mat.sim=round(rtriangle(1,a=age.mat2[1],b=age.mat2[2],c=floor(mean(age.mat2))))
        #age.mat.sim=floor(rtriangle(1,a=age.mat2[1],b=age.mat2[2]+1,c=age.mat2[1]))
      }
      Meanfec.SIM[[s]]=Meanfec.sim
      age.mat.SIM[s]=age.mat.sim
      
      #find age of 3m female shark 
      total.length=Lo.sim+(Linf.sim-Lo.sim)*(1-exp(-k.sim*age))
      Age.3M.sim=(log((1-((3-Lo.sim)/(Linf.sim-Lo.sim)))))/(-k.sim)
      Age.3M.SIM[s]=Age.3M.sim
      
      #temperature
      Aver.T.sim=rnorm(1,Aver.T,sd.Aver.T)    #normal dist
      Aver.T.SIM[s]=Aver.T.sim
      
       #total length
      total.length=Lo.sim+(Linf.sim-Lo.sim)*(1-exp(-k.sim*age))
      
      #fork length
      fork=(0.9442*total.length*100)-5.7441 # (Kohler et al 1996, length in cm)
      mid.FL.fem=fork
      
      #predicted size at maturity
      Pred.mat=exp(coef(glm1)[1]+coef(glm1)[2]*(total.length*100))/(1+exp(coef(glm1)[1]+coef(glm1)[2]*(total.length*100)))
      
      #natural mortality
      m.sim=M.fun(A=A.sim,k=k.sim,Linf=Linf.sim,Aver.T=Aver.T.sim,age.mat=age.mat.sim,
                  bwt,awt,W.G=1,W.noG=1)
      Compare.M[[s]]=m.sim$nat.mort.all
      m.sim=m.sim$nat.mort
      m.SIM[[s]]=m.sim
      
      #survivorship
      S=exp(-m.sim)         
      
      #proportion surviving
      lx=rep(NA,length(age))
      lx[1]=1.0
      for (i in 2:(length(age)))lx[i]=lx[i-1]*S[i]
      
      #reproductive schedules   
      #mx=Meanfec.sim*sexratio*Pmax.sim *Pred.mat
      MF=c(rep(0,(age.mat.sim-1)),Meanfec.sim[age.mat.sim:length(Meanfec.sim)])
      mx=MF*sexratio*Pmax.sim
      
      #probability of surviving (for birth-pulse, post-breeding census)
      px=vector(length=length(lx))
      for(i in 2:length(lx)) px[i-1]=(lx[i])/(lx[i-1])
      px.SIM[[s]]=px
      
      #fertility  (for birth-pulse, post-breeding census)
      bx=mx*px
      bx.SIM[[s]]=bx
      
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
      Proyec.matrix[[s]]=Data
    }
    
    #3. Solve projection matrix
    LAMBDA=sapply(Proyec.matrix,lambda)
    r.sim=log(LAMBDA)         
    
    t2.sim=log(2)/r.sim           #pop doubling time
    v=sapply(Proyec.matrix,reproductive.value)    #reproductive value
    Elasticities=sapply(Proyec.matrix,elasticity) #elasticities
    
      
    #4. Extract average elasticities per life stage
    Elast.stage=matrix(,ncol=n.elas.stages,nrow=length(Elasticities))
    nn=1:n.elas.stages
    if(scenario==1)
    {
      juv=juvenile.LH1
      adul1=adult1.LH1
      adul2=adult2.LH1
    }
    if(scenario==2)
    {
      juv=juvenile.LH2
      adul1=adult1.LH2
      adul2=adult2.LH2
      
    }
    
    
    for (e in 1:length(Elasticities))
    {
      test=Elasticities[[e]]
      Elast.stage[e,nn[1]]=sum(test[2:nrow(test),juv])
      Elast.stage[e,nn[2]]=sum(test[2:nrow(test),adul1])
      Elast.stage[e,nn[3]]=sum(test[2:nrow(test),adul2[adul2<=A.SIM[e]]])
      
      Elast.stage[e,nn[4]]=sum(test[1,juv])
      Elast.stage[e,nn[5]]=sum(test[1,adul1])
      Elast.stage[e,nn[6]]=sum(test[1,adul2[adul2<=A.SIM[e]]])
    }
    colnames(Elast.stage)=names.elast.stage 
    
    
    #5 Fill in missing age classes 
    if(scenario==1)
    {
      Meanfec.SIM=add.missing.age(Meanfec.SIM,RangeA.1[2])
      m.SIM=add.missing.age(m.SIM,RangeA.1[2])
      v.SIM=add.missing.age(v,RangeA.1[2])
    }
   
    
    if(scenario==2)
    {
      Meanfec.SIM=add.missing.age(Meanfec.SIM,RangeA.2[2])
      m.SIM=add.missing.age(m.SIM,RangeA.2[2])
      v.SIM=add.missing.age(v,RangeA.2[2])
    }
 
    #6. Remove  nonsense 
    #note: remove param combos that yield negative growth
    LAMBDA.reject=LAMBDA[LAMBDA<1]  #store negative growth separetely
    ID.lam.rej=which(LAMBDA<1)
    if(length(ID.lam.rej)>0)
    {
      LAMBDA=LAMBDA[-ID.lam.rej]        
      r.sim=r.sim[-ID.lam.rej]
      t2.sim=t2.sim[-ID.lam.rej]
      k.SIM=k.SIM[-ID.lam.rej]
      Linf.SIM=Linf.SIM[-ID.lam.rej]
      Reprod_cycle.SIM=Reprod_cycle.SIM[-ID.lam.rej]
      Pmax.SIM=Pmax.SIM[-ID.lam.rej]
      age.mat.SIM=age.mat.SIM[-ID.lam.rej]
      Age.3M.SIM=Age.3M.SIM[-ID.lam.rej]
      A.SIM=A.SIM[-ID.lam.rej]
      Aver.T.SIM=Aver.T.SIM[-ID.lam.rej]
      
      
      m.SIM=m.SIM[-ID.lam.rej,]
      Meanfec.SIM=Meanfec.SIM[-ID.lam.rej,]
      v.SIM=v.SIM[-ID.lam.rej,]
        
      Elasticities=Elasticities[-ID.lam.rej]
      Elast.stage=Elast.stage[-ID.lam.rej,]
      Proyec.matrix=Proyec.matrix[-ID.lam.rej]
        
    }
    
    #7. Store stuff
    output[[l]]=list(r=r.sim,t2=t2.sim,k=k.SIM,Linf=Linf.SIM,Rep.cycle=Reprod_cycle.SIM,
            Pmax=Pmax.SIM,Age.3M=Age.3M.SIM,age.mat=age.mat.SIM,A=A.SIM,Temp=Aver.T.SIM,
            m=m.SIM,meanfec=Meanfec.SIM,v.SIM=v.SIM,
            Elast.stage=Elast.stage,LAMBDA=LAMBDA,LAMBDA.reject=LAMBDA.reject,Proyec.matrix=
            Proyec.matrix,Compare.M=Compare.M)
    
  })

   #1.2. store scenarios 
output.by.year.1.1=output[[1]]   #scenario 1 of vital rates 
output.by.year.2.1=output[[2]]   #scenario 2 of vital rates 


  #1.3. extract values from list

    #1.3.1 check rejected par combos
Rejected.lambdas=c(Lambda.rej.1=length(output.by.year.1.1[LAMBDA.reject]),
                  Lambda.rej.2=length(output.by.year.2.1[LAMBDA.reject]))
print(Rejected.lambdas)

    #1.3.2 life history pars
ASim=output.by.year.2.1$A
ASim.1=output.by.year.1.1$A

AgeMatSim.1=output.by.year.1.1$age.mat
AgeMatSim=output.by.year.2.1$age.mat

Age.3MSim.1=output.by.year.1.1$Age.3M
Age.3MSim=output.by.year.2.1$Age.3M

MeanFecSim.1=output.by.year.1.1$meanfec
MeanFecSim=output.by.year.2.1$meanfec

PmaxSim=output.by.year.2.1$Pmax

kSim=output.by.year.2.1$k
LinfSim=output.by.year.2.1$Linf

kSim.1=output.by.year.1.1$k
LinfSim.1=output.by.year.1.1$Linf

TempSim=output.by.year.2.1$Temp

mSim=output.by.year.2.1$m
mSim.1=output.by.year.1.1$m

    #1.3.3 Leslie outputs
LAmBDA=output.by.year.2.1$LAMBDA
LAmBDA.1=output.by.year.1.1$LAMBDA
LAmBDA.reject=output.by.year.2.1$LAMBDA.reject
LAmBDA.reject.1=output.by.year.1.1$LAMBDA.reject
rSim=output.by.year.2.1$r
rSim.1=output.by.year.1.1$r
t2Sim=output.by.year.2.1$t2
t2Sim.1=output.by.year.1.1$t2
vSim=output.by.year.2.1$v.SIM
vSim.1=output.by.year.1.1$v.SIM
Elast.stage=output.by.year.2.1$Elast.stage
Elast.stage.1=output.by.year.1.1$Elast.stage
Proyec.matrix=output.by.year.2.1$Proyec.matrix
Proyec.matrix.1=output.by.year.1.1$Proyec.matrix

#free up memory
rm(list=c('output','output.by.year.1.1','output.by.year.2.1'))


#Set selectivity scenarios
Ns <- length(LAmBDA)    
iterations=1:Ns
Selectivity.SIM=Selectivity.SIM.1=vector("list",length = Ns)
scenario.sel=1
for (s in iterations) Selectivity.SIM.1[[s]]=Sel.fn(ASim.1[s],LinfSim.1[s],kSim.1[s])
scenario.sel=2
for (s in iterations) Selectivity.SIM[[s]]=Sel.fn(ASim[s],LinfSim[s],kSim[s])
SelSim.1=add.missing.age(Selectivity.SIM.1,RangeA.1[2])
SelSim=add.missing.age(Selectivity.SIM,RangeA.2[2])



#---6.  CATCH RECONSTRUCTIONS SECTION--- 
#note: use reconstructed catches for time 1 (1988-1996) to infer 1975 catches

# Years modeled
Yr.start=1838
Yr.end=Effort.WA$Year[nrow(Effort.WA)]

N.modl.yrs=Yr.start:Yr.end  #total years modelled
N.mdld=length(N.modl.yrs)

Yr.proj=1938:Yr.end    #years with catch                           
n.Yr.prol=length(Yr.proj)
Yr.burning=subset(N.modl.yrs,!(N.modl.yrs%in%(Yr.proj)))  #years burn in
n.burn.in=length(Yr.burning)  
These.Yrs=match(Yr.proj,N.modl.yrs)
id.1998=match(1998,N.modl.yrs)
id.2012=match(Yr.end,N.modl.yrs)

#Re set iterations, no need to do 10,000 (24 hours syst time) as 1,000 yields same result
Ns <- 1000    #number of simulations
iterations=1:Ns

#reconstructed WA catch periods
TIME.1=1988:1996
TIME.2=1997:2004
TIME.3=2005:Yr.end

NN.1=length(TIME.1)
NN.2=length(TIME.2)
NN.3=length(TIME.3)

#Catch history from WA TDGDLF and lobster, SA GHATF, marine scalefish, and gamefishers, Albany  
Effort.WA$Period=c(rep(NA,48),rep(1,NN.1),rep(2,NN.2),rep(3,NN.3))
Effort.GHATF$Period=c(rep(NA,48),rep(1,NN.1),rep(2,NN.2),rep(3,NN.3))
SA.scale$Period=c(rep(NA,5),rep(1,NN.1),rep(2,NN.2),rep(3,NN.3))
SA.scale.period$Period=1:3

WA.Eff.Per=aggregate(Total~Period,subset(Effort.WA,Year>1987),sum)
names(WA.Eff.Per)[2]="Agg.Eff"
GHATF.Eff.Per=aggregate(Total~Period,subset(Effort.GHATF,Year>1987),sum)
names(GHATF.Eff.Per)[2]="Agg.Eff"

  #Create lists for storing catches
TOT.Ktch_1.1.0=TOT.Ktch_1.2.0=TOT.Ktch_2.1.0=TOT.Ktch_2.2.0=vector('list',length=Ns)
Dropline_1.1.0=Dropline_1.2.0=Dropline_2.1.0=Dropline_2.2.0=TOT.Ktch_1.1.0
ByFishery_1.1.0=ByFishery_1.2.0=ByFishery_2.1.0=ByFishery_2.2.0=TOT.Ktch_1.1.0

  #calculate double the cpue for zone2
Recons.cpueX2=Recons.cpue
Recons.cpueX2$Mean.cpue=Recons.cpueX2$Mean.cpue*2
Recons.cpueX2$Low95.cpue=Recons.cpueX2$Low95.cpue*2
Recons.cpueX2$Up95.cpue=Recons.cpueX2$Up95.cpue*2

Recons.cpue2X2=Recons.cpue2
Recons.cpue2X2$Mean.cpue=Recons.cpue2X2$Mean.cpue*2
Recons.cpue2X2$Low95.cpue=Recons.cpue2X2$Low95.cpue*2
Recons.cpue2X2$Up95.cpue=Recons.cpue2X2$Up95.cpue*2


Export.Catches="NO"  #select if catches are exported or not

TOTAL.CATCH.PCS=vector('list',length(PCS.scenarios))
names(TOTAL.CATCH.PCS)=PCS.scenarios

system.time(for(pcs in 1:length(PCS.scenarios))      #(takes 0.086 seconds per iteration)
  {
    PCS=PCS.scenarios[pcs]
    for(i in 1:Ns)   
    {
      #Method 1, cpue Zn2 X 1
      a=fn.catch(Recons.catch,Recons.cpue,PCS)   
      TOT.Ktch_1.1.0[[i]]=a$TOT.Ktch
      Dropline_1.1.0[[i]]=a$Catch.Dropline
      ByFishery_1.1.0[[i]]=a
      
      #Method 1, cpue Zn2 X 2
      a=fn.catch(Recons.catch,Recons.cpueX2,PCS)
      TOT.Ktch_1.2.0[[i]]=a$TOT.Ktch
      Dropline_1.2.0[[i]]=a$Catch.Dropline
      ByFishery_1.2.0[[i]]=a
      
      #Method 2, cpue Zn2 X 1
      a=fn.catch(Recons.catch2,Recons.cpue2,PCS)
      TOT.Ktch_2.1.0[[i]]=a$TOT.Ktch
      Dropline_2.1.0[[i]]=a$Catch.Dropline
      ByFishery_2.1.0[[i]]=a
      
      #Method 2, cpue Zn2 X 2
      a=fn.catch(Recons.catch2,Recons.cpue2X2,PCS)
      TOT.Ktch_2.2.0[[i]]=a$TOT.Ktch
      Dropline_2.2.0[[i]]=a$Catch.Dropline
      ByFishery_2.2.0[[i]]=a
    }
    
    #Store and plot reconstructed total catch
    TOT.KTCH=list(TOT.Ktch_1.1.0=TOT.Ktch_1.1.0,TOT.Ktch_1.2.0=TOT.Ktch_1.2.0,
                  TOT.Ktch_2.1.0=TOT.Ktch_2.1.0,TOT.Ktch_2.2.0=TOT.Ktch_2.2.0)  
    
    TOTAL.CATCH=vector('list',length(TOT.KTCH))
    names(TOTAL.CATCH)=names(TOT.KTCH)
    
    for (j in 1:length(TOTAL.CATCH))
    {
      dum=TOT.KTCH[[j]]
      Tc=matrix(nrow=Ns,ncol=nrow(dum[[1]]))
      for(i in 1:Ns) Tc[i,]=dum[[i]]$TOT.CATCH
      TOTAL.CATCH[[j]]=Tc
    }
    
    #Extract dropline catch
    if(PCS=="NO")      #PCS is not applicable to dropline so only run it once
    {
      Dropline.KTCH=list(Dropline_1.1.0=Dropline_1.1.0,Dropline_1.2.0=Dropline_1.2.0,
                         Dropline_2.1.0=Dropline_2.1.0,Dropline_2.2.0=Dropline_2.2.0)
      Dropline.CATCH=vector('list',length(Dropline.KTCH))
      names(Dropline.CATCH)=names(Dropline.KTCH)    
      for (j in 1:length(Dropline.CATCH))
      {
        dum=Dropline.KTCH[[j]]
        Tc=matrix(nrow=Ns,ncol=length(dum[[1]]))
        for(i in 1:Ns) Tc[i,]=dum[[i]]
        Dropline.CATCH[[j]]=Tc
      }
    }
    
    YR=TOT.Ktch_1.1.0[[1]]$Year
    
    if(Export.Catches=="YES")
    {
      #export total catches
      for(i in 1:length(TOTAL.CATCH))write.csv(TOTAL.CATCH[[i]],
            paste("Catch scenarios/",KTCH.scen[i],"_PCS_",PCS.scenarios[pcs],".csv",sep=""))
      #export catch by fishery
      byFishery.KTCH=list(ByFishery_1.1.0=ByFishery_1.1.0,ByFishery_1.2.0=ByFishery_1.2.0,
                          ByFishery_2.1.0=ByFishery_2.1.0,ByFishery_2.2.0=ByFishery_2.2.0)    
      for (j in 1:length(byFishery.KTCH))
      {
        dum=byFishery.KTCH[[j]]
        
        MAT=matrix(nrow=Ns,ncol=length(dum[[1]]$TDGDLF))
        TDGDLF.ktch=GHATF.ktch=SA.marine.scale.ktch=MAT
        for(x in 1:Ns)
        {
          TDGDLF.ktch[x,]=dum[[x]]$TDGDLF
          GHATF.ktch[x,]=dum[[x]]$GHATF
          SA.marine.scale.ktch[x,]=dum[[x]]$SA.marine.scale
        }
        
        write.csv(TDGDLF.ktch,paste("Catch_for_Steve/TDGDLF.ktch.",names(byFishery.KTCH)[j],"_PCS_",PCS.scenarios[pcs],".csv",sep=""),row.names=F)
        write.csv(GHATF.ktch,paste("Catch_for_Steve/GHATF.ktch.",names(byFishery.KTCH)[j],"_PCS_",PCS.scenarios[pcs],".csv",sep=""),row.names=F)
        write.csv(SA.marine.scale.ktch,paste("Catch_for_Steve/SA.marine.scale.ktch.",names(byFishery.KTCH)[j],"_PCS_",PCS.scenarios[pcs],".csv",sep=""),row.names=F)
      }      
      rm(byFishery.KTCH)
    }
        
    TOTAL.CATCH.PCS[[pcs]]=TOTAL.CATCH
  })    



#---7.  RISK ANALYSES (Population projections) SECTION--- 
#Set Scenarios on initial and 1975 population sizes
NN=1
n.Yr.tot=length(YR)     
equil.yrs=1:16
  
#Free up memory
rm(Elasticities)
rm(bx.SIM)
rm(px.SIM)
rm(v.SIM)



#------Run Risk Analysis for each scenario------------


#1. REPORT  

#note: discrete scenarios on No
#Store outputs
LH1.N1.Sel1.Us.ab=vector('list',length=length(PCS.scenarios))
names(LH1.N1.Sel1.Us.ab)=PCS.scenarios
LH2.N1.Sel2.Us.ab=LH1.N1.Sel1.Us.ab

system.time(for(pcs in 1:length(PCS.scenarios))
{
  TOTAL.CATCH=TOTAL.CATCH.PCS[[pcs]]  
  STore.sim1=vector('list',length(TOTAL.CATCH))
  names(STore.sim1)=names(TOTAL.CATCH)
  STore.sim2=STore.sim1
  
  for(i in 1:length(TOTAL.CATCH))
  {      
    STR.sim1=vector('list',length(N.0.Rick))
    names(STR.sim1)=N.0.Rick
    STR.sim2=STR.sim1
    
    for (no in 1:length(N.0.Rick))
    {
      #Sample initial population size 
      N.0.sim=rep(N.0.Rick[no],Ns)
      
      #Run risk analysis
      STR.sim1[[no]]=Risk.fun(Life.hist.scenarios[[1]],TOTAL.CATCH[[i]],N.0.samp=N.0.sim)
      STR.sim2[[no]]=Risk.fun(Life.hist.scenarios[[2]],TOTAL.CATCH[[i]],N.0.samp=N.0.sim)  
    }
    #Store outputs
    STore.sim1[[i]]=STR.sim1       
    STore.sim2[[i]]=STR.sim2
  }
  LH1.N1.Sel1.Us.ab[[pcs]]=STore.sim1       
  LH2.N1.Sel2.Us.ab[[pcs]]=STore.sim2
  
})     #(takes 105 seconds per iteration)  
LH1.N1.Sel1.Us.ab.Report=LH1.N1.Sel1.Us.ab
LH2.N1.Sel2.Us.ab.Report=LH2.N1.Sel2.Us.ab

rm(LH1.N1.Sel1.Us.ab,LH2.N1.Sel2.Us.ab,STore.sim1,STore.sim2,STR.sim1,STR.sim2)



#2. PAPER 

#note: No sampled from distributions

Do.PAPER="NO"

  #Store outputs
LH1.N1.Sel1.Us.ab=vector('list',length=length(PCS.scenarios))
names(LH1.N1.Sel1.Us.ab)=PCS.scenarios
LH2.N1.Sel2.Us.ab=LH1.N1.Sel1.Us.ab

if(Do.PAPER=="YES")
{
  system.time(for(pcs in 1:length(PCS.scenarios))
  {
    TOTAL.CATCH=TOTAL.CATCH.PCS[[pcs]]  
    STore.sim1=vector('list',length(TOTAL.CATCH))
    names(STore.sim1)=names(TOTAL.CATCH)
    STore.sim2=STore.sim1
    
    for(i in 1:length(TOTAL.CATCH))
    {      
      STR.sim1=vector('list',length(N.0))
      names(STR.sim1)=N.0
      STR.sim2=STR.sim1
      
      for (no in 1:length(N.0))
      {
        #Sample initial population size 
        N.0.sim=round(runif(Ns,N.0[[no]][1],N.0[[no]][2]))
        
        #Run risk analysis
        STR.sim1[[no]]=Risk.fun(Life.hist.scenarios[[1]],TOTAL.CATCH[[i]],N.0.samp=N.0.sim)
        STR.sim2[[no]]=Risk.fun(Life.hist.scenarios[[2]],TOTAL.CATCH[[i]],N.0.samp=N.0.sim)  
      }
      #Store outputs
      STore.sim1[[i]]=STR.sim1       
      STore.sim2[[i]]=STR.sim2
    }
    LH1.N1.Sel1.Us.ab[[pcs]]=STore.sim1       
    LH2.N1.Sel2.Us.ab[[pcs]]=STore.sim2
    
  })       #(takes 101 seconds per iteration)  
  
}

#export list of model runs
# for(x in 1:length(PCS.scenarios))
# {
#   for(y in 1:length(KTCH.scen))  
#   {
#     for (z in 1:length(N.0.Rick))   
#     {
#       dput(LH1.N1.Sel1.Us.ab.Report[[x]][[y]][[z]], file = paste("Store.runs/LH1.N1.Sel1.Us.ab.Report_",x,"_",y,"_",z,".txt",sep=""))
#       dput(LH2.N1.Sel2.Us.ab.Report[[x]][[y]][[z]], file = paste("Store.runs/LH2.N1.Sel2.Us.ab.Report_",x,"_",y,"_",z,".txt",sep=""))
#     }
#     for (z in 1:length(N.0))   
#     {
#       dput(LH1.N1.Sel1.Us.ab[[x]][[y]][[z]], file = paste("Store.runs/LH1.N1.Sel1.Us.ab_",x,"_",y,"_",z,".txt",sep=""))
#       dput(LH2.N1.Sel2.Us.ab[[x]][[y]][[z]], file = paste("Store.runs/LH2.N1.Sel2.Us.ab_",x,"_",y,"_",z,".txt",sep=""))
#     }
#   }
# }


#read in list of model runs
#LH1.N1.Sel1.Us.ab <- dget(file = 'LH1.N1.Sel1.Us.ab.txt') 
#LH2.N1.Sel2.Us.ab <- dget(file = 'LH2.N1.Sel2.Us.ab.txt')



#---8.  SENSITIVITY TESTS---
#note: Sensitivity tests for informing future research
#     Max Age and Age maturity highly correlated so text Max Age only

Sens.pars=c("Catch","PCM","Selectivity","Max.Age","Rep.cycl","Pups","M")
Sens.Matrix=data.frame(Par=rep(Sens.pars,each=2))
Sens.Matrix$value=c("Min","Max","0","100",rep(c("Min","Max"),5))

N.0.Sens.Report=N.0.Rick
STORE.SENS.reprt=vector('list',length(N.0.Sens.Report))
names(STORE.SENS.reprt)=N.0.Sens.Report
Reference.report=STORE.SENS.reprt

for (no in 1:length(N.0.Sens.Report))
{
  STORE=vector('list',nrow(Sens.Matrix))
  names(STORE)=paste(Sens.Matrix$Par,Sens.Matrix$value)
  
  #Sensitivity tests
  for(s in 1:nrow(Sens.Matrix)) STORE[[s]]=fun.sensitivity(Sens.Matrix[s,],N.init=N.0.Sens.Report[no])
  STORE.SENS.reprt[[no]]=STORE
  
  #Base case
  Reference.report[[no]]=fun.sensitivity(data.frame(Par="NIL",value="NIL"),N.init=N.0.Sens.Report[no])
}

  #store values of tested quantities
STORE.SENS.inputs=Reference.inputs=STORE.SENS.reprt
for (no in 1:length(N.0.Sens.Report))
{
  STORE=vector('list',nrow(Sens.Matrix))
  names(STORE)=paste(Sens.Matrix$Par,Sens.Matrix$value)
  
  #Sensitivity tests
  for(s in 1:nrow(Sens.Matrix)) STORE[[s]]=fun.sensitivity.get.input.values(SenSc=Sens.Matrix[s,],N.init=N.0.Sens.Report[no])
  STORE.SENS.inputs[[no]]=STORE
  
  #Base case
  Reference.inputs[[no]]=fun.sensitivity.get.input.values(data.frame(Par="NIL",value="NIL"),N.init=N.0.Sens.Report[no])
}




#---9.  DISPLAY RESULTS SECTION---  

plot.power.point="NO"
#plot.power.point="YES"

Ax.cx=1.35
Lab.xc=2
LWD=3

  #colors for paper
LH1.col="grey70"
LH2.col="grey20"

if(plot.power.point=="NO")
{
  setwd(handl_OneDrive("Analyses/Demography/White shark/Outputs"))
  greyscale=c(LH1.col,LH2.col)
  Basic="black"
  line.col=gray.colors(NN,start=0,end=0.8)
  COLS.prob=rep(Basic,3)
}

#colors for power point presentation
if(plot.power.point=="YES")
{
  setwd(handl_OneDrive("Analyses/Demography/White shark/Outputs/Power.point"))
  greyscale=c("red","chartreuse4")  
  Basic="blue"
  #line.col=c("black","red","chartreuse4","blue","mediumorchid1")
  line.col=c(heat.colors(length(N.0)/2),rev(topo.colors(length(N.0)/2)))
  fun.back=function()rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(.2, .2, .2,0.15))
  COLS.prob=c(greyscale,Basic)
}


#9.1 DEMOGRAPHIC ANALYSIS OUTPUTS

  #Figure 1. plot priors (age invariant)
tiff(file="Figure1.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(3,2),las=1,mai=c(.55,.3,.1,0),omi=c(.1,.2,.1,.2))


  #1 scenarios for maximum age
fn.barplo(ASim.1,ASim,c(RangeA.1[1],RangeA.1[2]),c(RangeA.2[1]+1,RangeA.2[2]),"Maximum age (years)")
if(plot.power.point=="YES") legend("topright",c("LH1","LH2"),bty='n',cex=2,fill=greyscale)
if(plot.power.point=="YES") fun.back()  
  
  #2 scenarios for maturity           
fn.barplo(AgeMatSim.1,AgeMatSim,age.mat,age.mat2,"Age at maturity (years)")
if(plot.power.point=="YES") fun.back()

  #3 scenarios for fecundity
table.MeanFec.1=table(MeanFecSim.1[,1])/sum(table(MeanFecSim.1[,1]))
table.MeanFec=table(MeanFecSim[,1])/sum(table(MeanFecSim[,1]))
if(length(table.MeanFec)<length(table.MeanFec.1))
{
  add.names=names(table.MeanFec.1)[-match(names(table.MeanFec),names(table.MeanFec.1))]
  start.add=(length(table.MeanFec)+1)
  end.add=length(table.MeanFec.1)
  table.MeanFec[start.add:end.add]=0
  names(table.MeanFec)[start.add:end.add]=add.names
  table.MeanFec=table.MeanFec[match(sort(as.numeric(names(table.MeanFec))),
                                    as.numeric(names(table.MeanFec)))]
  
}

combined.Fec=rbind(table.MeanFec.1,table.MeanFec)
barplot(combined.Fec, beside = TRUE,ylim=c(0,max(combined.Fec)*1.1),mgp = c(2.5, 0.6, 0),
        names.arg= Rangefec,xlab="Number of pups",yaxt="n",
        ylab="", axis.lty=1, axes=T,col=greyscale,cex.names=Ax.cx,las=1,cex.lab=Lab.xc,space=c(.1,.1))
box()
if(plot.power.point=="YES") fun.back()

  #4 scenarios for reproductive cycle
BarPlot(PmaxSim,"Reproductive cycle length (years)",c("3","2","1"))
if(plot.power.point=="YES") fun.back()


  #5 scenarios for k
LinesPlot=function(Data,LineType,Cols)lines(density(Data, adjust=2), lty=LineType,col=Cols,lwd=2.5)

plot(density(kSim, adjust=2), type="l",lty=1,col=Basic,xlab=expression("k " (years^-1)),ylab="",main="",
       yaxt="n",cex.axis=Ax.cx,cex.lab=Lab.xc,lwd=LWD)
#LinesPlot(kSim,1,greyscale[2])
#LinesPlot(kSim.1,1,greyscale[1])
if(plot.power.point=="YES") fun.back()

  #6 scenarios for Linf
plot(density(LinfSim, adjust=2), type="l",lty=1,col=Basic,xlab="Linf (m)",ylab="",main="",
     yaxt="n",cex.axis=Ax.cx,cex.lab=Lab.xc,lwd=LWD)
if(plot.power.point=="YES") fun.back()
# plot(density(toSim, adjust=2), type="l",lty=1,col="white",xlab="to (years)",ylab="",main="",
#        xlim=c(-6,-1),ylim=c(0,0.8),yaxt="n",cex.axis=1.1,cex.lab=1.4)
# LinesPlot(toSim,1,greyscale[2])
# LinesPlot(toSim.1,1,greyscale[1])

# plot(density(L50Sim, adjust=2), type="l",lty=2,lwd=1.75,col="white",xlab="Total length (m)",ylab="",main="",
#        xlim=c(4,6),ylim=c(0,2.5),yaxt="n",cex.axis=1.1,cex.lab=1.4)
# LinesPlot(L50Sim,2,1)
# LinesPlot(L95Sim,3,1)
#LinesPlot(LinfSim.1,1,1)
#LinesPlot(LinfSim,1,greyscale[2])
# 
# box()
# legend("topright", c("L50","L95"), lty=c(2,3),col=c(1,1),bty="n",lwd=1.75,cex=1.2)

mtext("Distribution",side=2,outer=T,line=-0.75,font=1,las=0,cex=2)

dev.off()

  #Scenarios for Temperature
plot(density(TempSim, adjust=2), type="l",lty=1,col="white",xlab="Temperature (?C)",ylab="",main="",
       xlim=c(5,25),ylim=c(0,0.285),yaxt="n",cex.axis=1,cex.lab=1.4)
LinesPlot(TempSim,1,1)
#expression(Temperature~(~degree~C))




  #Figure 2. plot priors (at-age schedules)
M.median.range.LH2=round(c(median(mSim[,1]),median(mSim[,ncol(mSim)],na.rm=T)),3)
M.median.range.LH1=round(c(median(mSim.1[,1]),median(mSim.1[,ncol(mSim.1)],na.rm=T)),3)

    #create age vector
AGE.1=factor(rep((first.age+1):RangeA.1[2],each=nrow(mSim.1)),levels=first.age:RangeA.1[2])
AGE.2=factor(rep((first.age+1):RangeA.2[2],each=nrow(mSim)),levels=first.age:RangeA.2[2])
AGE.sel.1=factor(rep((first.age+1):RangeA.1[2],each=nrow(SelSim.1)),levels=first.age:RangeA.1[2])
AGE.sel.2=factor(rep((first.age+1):RangeA.2[2],each=nrow(SelSim)),levels=first.age:RangeA.2[2])


tiff(file="Figure2.1.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(2,1),las=1,mai=c(.2,.9,.3,.05),omi=c(.5,.1,.1,.05),mgp=c(5,0.75,0))

#M    
boxplot(c(mSim)~AGE.2,col="white",ylab="",xaxt="n",cex.axis=1.5,cex.lab=2,outline=F,
        border=greyscale[2],varwidth=T,ylim=c(0,0.21))
boxplot(c(mSim.1)~AGE.1,col="white",ylab="",xaxt="n",cex.axis=1.5,cex.lab=2,outline=F,
        border=greyscale[1],add=T,varwidth=T)

axis(1,at=c(seq(2,A.2,2)),labels=F,tck=-0.03)
axis(1,at=c(seq(1,A.2,2)),labels=F,tck=-0.015)
mtext(expression("Natural mortality  " (year^-1)),side=2,line=3.25,font=1,las=0,cex=1.65)
if(plot.power.point=="YES") legend("topright",c("LH1","LH2"),bty='n',cex=2,fill=greyscale)
if(plot.power.point=="YES") fun.back()

#Selectivity
boxplot(c(SelSim)~AGE.sel.2,col="white",xaxt="n",cex.axis=1.5,cex.lab=2,border=greyscale[2],outline=F,
        varwidth=T,ylim=c(0,1))
 boxplot(c(SelSim.1)~AGE.sel.1,col="white",ylab="",xaxt="n",cex.axis=1.5,cex.lab=2,border=greyscale[1],
         add=T,outline=F,varwidth=T)
axis(1,at=c(seq(2,A.2,2)),labels=F,tck=-0.03)
axis(1,at=c(seq(1,A.2,2)),labels=F,tck=-0.015)
mtext("Relative selectivity",side=2,line=3.75,font=1,las=0,cex=1.65)
mtext("Age (years)",side=1,line=2,font=1,las=0,cex=1.65)
axis(1,at=seq(2,A.2,2),labels=seq(2,A.2,2),cex.axis=1.35,las=1,tck=-0.03)
if(plot.power.point=="YES") fun.back()
dev.off()


#Table 1. Stats (95% CI are mean + and - CI)
este.r=LAmBDA
este.r=este.r[!(is.na(este.r))]
rStats=c(mean=mean(este.r),SD=sd(este.r),CV=sd(este.r)/mean(este.r),
         CI95=sort(este.r)[c(floor(0.025*length(este.r)),
                             ceiling(0.975*length(este.r)))],median=median(este.r))

este.r.1=LAmBDA.1
este.r.1=este.r.1[!(is.na(este.r.1))]
rStats.1=c(mean=mean(este.r.1),SD=sd(este.r.1),CV=sd(este.r.1)/mean(este.r.1),
           CI95=sort(este.r.1)[c(floor(0.025*length(este.r.1)),
                                 ceiling(0.975*length(este.r.1)))],median=median(este.r.1))

este.t2=t2Sim
este.t2=este.t2[!(is.na(este.t2))]
t2Stats=c(mean=mean(este.t2),SD=sd(este.t2),CI95=sort(este.t2)[c(floor(0.025*length(este.t2)),
                                   ceiling(0.975*length(este.t2)))],median=median(este.t2))

este.t2.1=t2Sim.1
este.t2.1=este.t2.1[!(is.na(este.t2.1))]
t2Stats.1=c(mean=mean(este.t2.1),SD=sd(este.t2.1),CI95=sort(este.t2.1)[c(floor(0.025*length(este.t2.1)),
                                     ceiling(0.975*length(este.t2.1)))],median=median(este.t2.1))


#Export population params summaries
write.table(rStats,"LambdaStats.scen2.csv",sep=",");write.table(rStats.1,"LambdaStats.scen1.csv",sep=",")
write.table(t2Stats,"t2Stats.scen2.csv",sep=",");write.table(t2Stats.1,"t2Stats.scen1.csv",sep=",")


#Figure 3
#Plot population parameters
tiff(file="Figure3.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
if(plot.power.point=="NO")
{
  par(mfcol=c(3,1),las=1,mai=c(.45,.5,.075,.5),omi=c(.1,.4,.1,.5),mgp=c(1.5,0.6,0))
  #  par(mfcol=c(2,2),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
  
  #Lambda
  plot(density(LAmBDA.1, adjust=2,na.rm =T), type="l",lty=1,col=greyscale[1],xlab="",lwd=LWD,
       ylab="",main="",xlim=c(1,1.2),ylim=c(0,25),yaxt="n",cex.axis=1.75,cex.lab=1.65)
  lines(density(LAmBDA, adjust=2,na.rm =T), lty=1,col=greyscale[2],lwd=LWD)
  #legend("topright",c("what",""),cex=1.2,bty="n",col=c(1,1))
  mtext(expression(paste(lambda," " (year^-1),sep="")),
        #mtext(expression(paste("Population growth rate, ",lambda," " (years^-1),sep="")),
        side=1,line=2.5,font=1,las=0,cex=1.5)
  #t2
  plot(density(t2Sim.1, adjust=3,na.rm =T,from=0,to=60), type="l",lty=1,col=greyscale[1],xlab="",lwd=LWD,
       ylab="",main="",xlim=c(0,60),yaxt="n",cex.axis=1.75,cex.lab=1.65,ylim=c(0,.105))
  lines(density(t2Sim, adjust=3,na.rm =T,from=0,to=60), lty=1,col=greyscale[2],lwd=LWD)
  mtext("TD (years)",side=1,line=2.5,font=1,las=0,cex=1.5)
  mtext("                                 Distribution",side=2,line=3,font=1,las=0,cex=1.8)
  
  #reproductive value
  vSim.mean=colMeans(vSim, na.rm = T)
  vSim.upCI=colQuantiles(vSim,probs=0.975, na.rm = T)
  vSim.lowCI=colQuantiles(vSim,probs=0.025, na.rm = T)
  vSim.1.mean=colMeans(vSim.1, na.rm = T)
  vSim.1.upCI=colQuantiles(vSim.1,probs=0.975, na.rm = T)
  vSim.1.lowCI=colQuantiles(vSim.1,probs=0.025, na.rm = T)
  
  agE=1:RangeA.2[2]
  plot(agE,vSim.mean, type='l',col=greyscale[2],xlab="", lwd=LWD,ylim=c(0,14),
       ylab="",cex.axis=1.75,cex.lab=1.65)
  # segments(agE,vSim.mean,agE,vSim.upCI,col="gray55", lwd=2)
  # segments(agE,vSim.mean,agE,vSim.lowCI,col="gray55", lwd=2)
  lines(agE[1:RangeA.1[2]],vSim.1.mean, pch=16,col=greyscale[1],lwd=LWD)
  # segments(agE,vSim.1.mean,agE,vSim.1.upCI,col=1, lwd=2)
  # segments(agE,vSim.1.mean,agE,vSim.1.lowCI,col=1, lwd=2)
  mtext("Age (years)",side=1,line=2.5,font=1,las=0,cex=1.5)
  mtext("Reproductive value",side=2,line=3,font=1,las=0,cex=1.8)
  
  
  #elasticities
  #   table.Elast=colMeans(Elast.stage, na.rm = T)
  #   table.Elast=table.Elast/sum(table.Elast)  #standardise
  #   table.Elast.1=colMeans(Elast.stage.1, na.rm = T)
  #   table.Elast.1=table.Elast.1/sum(table.Elast.1)
  #   RangeAge.group=rep(c("Juv","Ad1","Ad2"),2)
  #   
  #   combined.Elast=rbind(table.Elast.1,table.Elast)
  #   barplot(combined.Elast, beside = TRUE,ylim=c(0,max(combined.Elast)*1.1),mgp = c(2.5, 0.6, 0),
  #           names.arg= RangeAge.group,xlab="",ylab="",
  #           axis.lty=1, axes=T,col=greyscale,cex.names=1.05,las=1,cex.axis=1.2,cex.lab=1.65,space=c(.0,.99))
  #   box()
  #   abline(v=9.46, col =1,lty=2)
  #   legend('topleft','survival',bty='n',cex=1.25)
  #   legend('topright','fecundity',bty='n',cex=1.25)
  #   mtext("Age group",side=1,line=2,font=1,las=0,cex=1.3)
  #   mtext("Elasticities",side=2,line=2.2,font=1,las=0,cex=1.3)
  
}
if(plot.power.point=="YES")
{
  par(mfcol=c(1,1),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
  
  #Lambda
  plot(density(LAmBDA.1, adjust=2,na.rm =T), type="l",lty=1,col=greyscale[1],xlab="",lwd=LWD,
       ylab="",main="",xlim=c(1,1.2),ylim=c(0,35),yaxt="n",cex.axis=1.3,cex.lab=1.65)
  lines(density(LAmBDA, adjust=2,na.rm =T), lty=1,col=greyscale[2],lwd=LWD)
  mtext(expression(paste(lambda," " (year^-1),sep="")),
        side=1,line=2.1,font=1,las=0,cex=1.5)
  fun.back()
  text(rStats[1]*1.05,34,
       paste(round(rStats[6],3),"= ",round((rStats[6]-1)*100,1),"%",
             " (pop. doubling=",round(t2Stats[5])," years)", sep=""),
       col=greyscale[2],cex=1.5)
  text(rStats.1[1]*1.05,19,
       paste(round(rStats.1[6],3),"= ",round((rStats.1[6]-1)*100,1),"%",
             " (pop. doubling=",round(t2Stats.1[5])," years)", sep=""),
       col=greyscale[1],cex=1.5)
  legend("topright",c("LH1","LH2"),bty='n',lty=1,col=greyscale,cex=1.5,lwd=LWD)
  
  mtext("Distribution",side=2,line=-1,font=1,las=0,cex=2,outer=T)
}
dev.off()



#9.2 RISK ANALYSIS OUTPUTS

#Figure 2. Catch time series
#note: show full catches (i.e. with no PCS)
tiff(file="Recons.Catch.CI.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
YMAX=1500
par(mfcol=c(2,2),las=1,mai=c(.1,.2,.075,.15),omi=c(.5,.7,.1,.05),mgp=c(1.5,0.9,0))
for (j in 1:length(TOTAL.CATCH.PCS$NO))
{
  #fn.plot.ktch(TOTAL.CATCH[[j]],XLIM=c(1938,Yr.end))
  fn.plot.ktch.percentile(TOTAL.CATCH.PCS$NO[[j]],XLIM=c(1938,Yr.end))
  if(j%in%c(2,4))axis(1,at=seq(1940,2010,20),labels=c("1940","1960","1980","2000"),tck=-0.06,cex.axis=1.25)
  if(j%in%1:2) axis(2,at=c(seq(0,YMAX,250)),labels=c(seq(0,YMAX,250)),tck=-0.03,cex.axis=1.5)
  #if(j==1)   legend("topleft",c("50th percentile","75th percentile","95th percentile"),bty='n',cex=1.5,fill=COLS)
  legend("topright",KTCH.scen[j],bty='n',cex=1.85)
  #if(j==1)mtext("Method 1",3,0,cex=1.5)
  #if(j==5)mtext("Method 2",3,0,cex=1.5)
}
mtext("Total catch (in numbers)",2,2.5,outer=T,cex=1.5,las=3)
mtext("Financial year",1,1.75,outer=T,cex=1.5)
dev.off()


#get values from catch series
get.ktch.per=function(DAT)
{
  #get percentiles
  Low.percentile=function(Nper) apply(DAT, 2, function(x) quantile(x, (0+Nper)/100))
  High.percentile=function(Nper) apply(DAT, 2, function(x) quantile(x, (100-Nper)/100))
  
  #50% of data
  Nper=(100-50)/2
  LOW.50=Low.percentile(Nper)
  UP.50=High.percentile(Nper)
  
  #75% of data
  Nper=(100-75)/2
  LOW.75=Low.percentile(Nper)
  UP.75=High.percentile(Nper)
  
  #95% of data
  Nper=(100-95)/2
  LOW.95=Low.percentile(Nper)
  UP.95=High.percentile(Nper)
  
  #100% of data
  Nper=(100-100)/2
  LOW.100=Low.percentile(Nper)
  UP.100=High.percentile(Nper)
  
  return(list(MAX.range=c(max(LOW.50),max(UP.50)),MAXUncertainty=c(max(LOW.95),max(UP.95)),
         Current.50=c(LOW.50[length(LOW.50)],UP.50[length(UP.50)])))
  
}
get.ktch.per(TOTAL.CATCH.PCS$NO[[2]])
get.ktch.per(TOTAL.CATCH.PCS$NO[[3]])


#Figure 4. Population projections

#--- 4.1 REPORT ----
LH1.col="grey80"
LH2.col="grey50"
plot.cols="NO"
if(plot.cols=="NO")
{
  colfunc <- colorRampPalette(c(LH2.col,LH1.col))
  LIN="black"
  Brdr ="grey20"
}  
if(plot.cols=="YES")
{
  colfunc <- colorRampPalette(c("lightblue1","deepskyblue2"))  
  Brdr="dodgerblue4"
  LIN="dodgerblue4"
}

#4.1.1. Total numbers    
Store.Report.median.all=vector('list',length(PCS.scenarios))
names(Store.Report.median.all)=PCS.scenarios
names(Store.Report.median.all)=ifelse(names(Store.Report.median.all)=="YES","PCM_50",
                    ifelse(names(Store.Report.median.all)=="YES.100","PCM_0","PCM_100"))
for(x in 1:length(PCS.scenarios))
{
  tiff(file=paste("Figure4_Report_total_",names(Store.Report.median.all)[x],".tiff",sep=""),
       width = 2400, height = 2000,units = "px", res = 300, compression = "lzw")   
  
  par(mfcol=c(5,4),las=1,mai=c(.1,.1,.1,0.01),omi=c(.45,.55,.2,.01),mgp=c(1.5,0.6,0),xpd=F)
  STR.LH.y=vector('list',length(KTCH.scen))
  names(STR.LH.y)=KTCH.scen
  
  for(y in 1:length(KTCH.scen))  
  {
    STR.LH.z=vector('list',length(N.0.Rick))
    names(STR.LH.z)=N.0.Rick
    for (z in 1:length(N.0.Rick))   
    {
      STR.LH.z[[z]]=plot.mean.CI.Report(LH1.N1.Sel1.Us.ab.Report[[x]][[y]][[z]]$Rel.abundance,
                    LH2.N1.Sel2.Us.ab.Report[[x]][[y]][[z]]$Rel.abundance,CEXX=1.15,Fem_3M="NO","bottomleft")
      if(z==1)mtext(KTCH.scen[y],side=3,line=0,font=1,las=0,cex=1.75)      
    }
    
    STR.LH.y[[y]]=STR.LH.z
    
  }
  
  Store.Report.median.all[[x]]=STR.LH.y
  mtext("Financial year",side=1,outer=T,line=1.95,font=1,las=0,cex=1.75)
  mtext("Number of females",side=2,outer=T,line=2.2,font=1,las=0,cex=1.75)
  dev.off()
}

#Table of rate of increase and depletion level
Store.Report.Depltn.Inc=Store.Report.median.all
for(x in 1:length(PCS.scenarios))
{
  STR.LH.y=vector('list',length(KTCH.scen))
  names(STR.LH.y)=KTCH.scen
  
  for(y in 1:length(KTCH.scen))  
  {
    STR.LH.z=vector('list',length(N.0.Rick))
    names(STR.LH.z)=N.0.Rick
    for (z in 1:length(N.0.Rick))   
    {
      STR.LH.z[[z]]=fn.inc.rate.depl(LH1.N1.Sel1.Us.ab.Report[[x]][[y]][[z]]$Rel.abundance,
                                        LH2.N1.Sel2.Us.ab.Report[[x]][[y]][[z]]$Rel.abundance)
    }
    
    STR.LH.y[[y]]=STR.LH.z
    
  }
  
  Store.Report.Depltn.Inc[[x]]=STR.LH.y

}

TBL.Rate.inc.1998=as.data.frame(matrix(nrow=length(KTCH.scen)*length(PCS.scenarios),ncol=12))
colnames(TBL.Rate.inc.1998)=c("PCM","Catch",paste(rep(N.0.Rick,each=2),c("LH1","LH2")))
TBL.Rate.inc.1998$PCM=rep(names(Store.Report.Depltn.Inc),each=4)
TBL.Rate.inc.1998$Catch=rep(KTCH.scen,length(PCS.scenarios))
TBL.Depltn=TBL.Rate.inc.1998
for(x in 1:length(PCS.scenarios))
{
 for(y in 1:length(KTCH.scen))
  {
   a=round(c(sapply(Store.Report.Depltn.Inc[[x]][[y]], "[[", 2)),1)
   b=round(c(sapply(Store.Report.Depltn.Inc[[x]][[y]], "[[", 1)),2)
   q=which(TBL.Rate.inc.1998$PCM==names(Store.Report.Depltn.Inc)[x] & TBL.Rate.inc.1998$Catch==KTCH.scen[y])
   TBL.Rate.inc.1998[q,3:12]=a
   TBL.Depltn[q,3:12]=b
  }
}

  #export rate of increase since protection
fn.word.table(WD=getwd(),TBL=TBL.Rate.inc.1998,Doc.nm="Table.Increase.since.1998.Report",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")

  #export depletion level
fn.word.table(WD=getwd(),TBL=TBL.Depltn,Doc.nm="Table.Depletion.Report",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")


#4.1.2 Females larger than 3 m
  #plot trajectories
Store.Report.median.3m=Store.Report.median.all
for(x in 1:length(PCS.scenarios))
{
  tiff(file=paste("Figure4_Report_3m_",names(Store.Report.median.3m)[x],".tiff",sep=""),
       width = 2400, height = 2000,units = "px", res = 300, compression = "lzw")   
  
  par(mfcol=c(5,4),las=1,mai=c(.1,.1,.1,0.01),omi=c(.45,.55,.2,.01),mgp=c(1.5,0.6,0),xpd=F)
  STR.LH.y=vector('list',length(KTCH.scen))
  names(STR.LH.y)=KTCH.scen
  
  for(y in 1:length(KTCH.scen))  
  {
    STR.LH.z=vector('list',length(N.0.Rick))
    names(STR.LH.z)=N.0.Rick
    for (z in 1:length(N.0.Rick))   
    {
      YLM3M=N.0.Rick[z]
      STR.LH.z[[z]]=plot.mean.CI.Report(LH1.N1.Sel1.Us.ab.Report[[x]][[y]][[z]]$A3M.rel.abundance,
                                        LH2.N1.Sel2.Us.ab.Report[[x]][[y]][[z]]$A3M.rel.abundance,
                                        CEXX=1.25,Fem_3M="YES","topright")
      if(z==1)mtext(KTCH.scen[y],side=3,line=0,font=1,las=0,cex=1.75)      
    }
    
    STR.LH.y[[y]]=STR.LH.z
    
  }
  
  Store.Report.median.3m[[x]]=STR.LH.y
  mtext("Financial year",side=1,outer=T,line=1.95,font=1,las=0,cex=1.75)
  mtext("Number of females > 3m TL",side=2,outer=T,line=2.2,font=1,las=0,cex=1.75)
  dev.off()
}

  #extract depletion and increase rate
for(x in 1:length(PCS.scenarios))
{
  STR.LH.y=vector('list',length(KTCH.scen))
  names(STR.LH.y)=KTCH.scen
  
  for(y in 1:length(KTCH.scen))  
  {
    STR.LH.z=vector('list',length(N.0.Rick))
    names(STR.LH.z)=N.0.Rick
    for (z in 1:length(N.0.Rick))   
    {
      STR.LH.z[[z]]=fn.inc.rate.depl(LH1.N1.Sel1.Us.ab.Report[[x]][[y]][[z]]$A3M.rel.abundance,
                                     LH2.N1.Sel2.Us.ab.Report[[x]][[y]][[z]]$A3M.rel.abundance)
    }
    
    STR.LH.y[[y]]=STR.LH.z
    
  }
  
  Store.Report.Depltn.Inc[[x]]=STR.LH.y
  
}


TBL.Depltn_3m=TBL.Depltn
TBL.Rate.inc.1998_3m=TBL.Rate.inc.1998
for(x in 1:length(PCS.scenarios))
{
  for(y in 1:length(KTCH.scen))
  {
    a=round(c(sapply(Store.Report.Depltn.Inc[[x]][[y]], "[[", 2)),1)
    b=round(c(sapply(Store.Report.Depltn.Inc[[x]][[y]], "[[", 1)),2)
    q=which(TBL.Rate.inc.1998$PCM==names(Store.Report.Depltn.Inc)[x] & TBL.Rate.inc.1998$Catch==KTCH.scen[y])  
    TBL.Rate.inc.1998_3m[q,3:12]=a
    TBL.Depltn_3m[q,3:12]=b
  }
}
fn.word.table(WD=getwd(),TBL=TBL.Rate.inc.1998_3m,Doc.nm="Table.Increase.since.1998.Report_3m",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")

fn.word.table(WD=getwd(),TBL=TBL.Depltn_3m,Doc.nm="Table.Depletion.Report_3m",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")




#--- 4.2 PAPER ----
if(Do.PAPER=="YES")
{
  plot.cols="NO"
  if(plot.cols=="NO")
  {
    colfunc <- colorRampPalette(c(LH2.col,LH1.col))
    LIN="black"
    Brdr ="grey20"
  }  
  if(plot.cols=="YES")
  {
    colfunc <- colorRampPalette(c("lightblue1","deepskyblue2"))  
    Brdr="dodgerblue4"
    LIN="dodgerblue4"
  }
  #Total biomass
  Store.Paper.median.all=Store.Report.median.all
  for(x in 1:length(PCS.scenarios))
  {
    tiff(file=paste("Figure4_Paper_total_",names(Store.Paper.median.all)[x],".tiff",sep=""),
         width = 2600, height = 2000,units = "px", res = 300, compression = "lzw")   
    
    par(mfcol=c(5,4),las=1,mai=c(.1,.1,.1,0.01),omi=c(.45,.8,.2,.01),mgp=c(1.5,0.6,0),xpd=F)
    STR.LH.y=vector('list',length(KTCH.scen))
    names(STR.LH.y)=KTCH.scen
    
    for(y in 1:length(KTCH.scen))  
    {
      STR.LH.z=vector('list',length(N.0))
      names(STR.LH.z)=N.0
      for (z in 1:length(N.0))   
      {
        STR.LH.z[[z]]=plot.mean.CI.Paper(LH1.N1.Sel1.Us.ab[[x]][[y]][[z]]$Rel.abundance,
                                         LH2.N1.Sel2.Us.ab[[x]][[y]][[z]]$Rel.abundance,CEXX=1.25)
        if(z==1)mtext(KTCH.scen[y],side=3,line=0,font=1,las=0,cex=1.75)      
        if(y==1)mtext(paste("[",N.0[z][[1]][1],",",N.0[z][[1]][2],"]",sep=""),side=2,line=5,font=1,las=0,cex=1)
      }
      
      STR.LH.y[[y]]=STR.LH.z
      
    }
    
    Store.Paper.median.all[[x]]=STR.LH.y
    mtext("Financial year",side=1,outer=T,line=1.95,font=1,las=0,cex=1.75)
    mtext("Relative number of females",side=2,outer=T,line=2,font=1,las=0,cex=1.75)
    dev.off()
  }
  
  #Table of rate of increase and depletion level  
  Store.Report.Depltn.Inc.Paper=Store.Report.median.all
  for(x in 1:length(PCS.scenarios))
  {
    STR.LH.y=vector('list',length(KTCH.scen))
    names(STR.LH.y)=KTCH.scen
    
    for(y in 1:length(KTCH.scen))  
    {
      STR.LH.z=vector('list',length(N.0))
      names(STR.LH.z)=N.0
      for (z in 1:length(N.0))   
      {
        STR.LH.z[[z]]=fn.inc.rate.depl(LH1.N1.Sel1.Us.ab[[x]][[y]][[z]]$Rel.abundance,
                                       LH2.N1.Sel2.Us.ab[[x]][[y]][[z]]$Rel.abundance)
      }
      
      STR.LH.y[[y]]=STR.LH.z
      
    }
    
    Store.Report.Depltn.Inc.Paper[[x]]=STR.LH.y
    
  }
  
  TBL.Rate.inc.1998=as.data.frame(matrix(nrow=length(KTCH.scen)*length(PCS.scenarios),ncol=12))
  colnames(TBL.Rate.inc.1998)=c("PCM","Catch",paste(rep(N.0,each=2),c("LH1","LH2")))
  TBL.Rate.inc.1998$PCM=rep(names(Store.Report.Depltn.Inc.Paper),each=4)
  TBL.Rate.inc.1998$Catch=rep(KTCH.scen,length(PCS.scenarios))
  TBL.Depltn=TBL.Rate.inc.1998
  for(x in 1:length(PCS.scenarios))
  {
    for(y in 1:length(KTCH.scen))
    {
      a=round(c(sapply(Store.Report.Depltn.Inc.Paper[[x]][[y]], "[[", 2)),1)
      b=round(c(sapply(Store.Report.Depltn.Inc.Paper[[x]][[y]], "[[", 1)),2)
      q=which(TBL.Rate.inc.1998$PCM==names(Store.Report.Depltn.Inc)[x] & TBL.Rate.inc.1998$Catch==KTCH.scen[y]) 
      TBL.Rate.inc.1998[q,3:12]=a
      TBL.Depltn[q,3:12]=b
    }
  }
  fn.word.table(WD=getwd(),TBL=TBL.Rate.inc.1998,Doc.nm="Table.Increase.since.1998.Paper",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
  fn.word.table(WD=getwd(),TBL=TBL.Depltn,Doc.nm="Table.Depletion.Paper",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
}




#9.3 SENSITIVITY TESTS
  #Show population trajectories
AT=seq(match(1940,Yr.proj),n.Yr.prol,by=20)
CL.ppr=c("black","grey40","grey70")
LTY.ppr=1:3
for(p in 1:length(Reference.report))
{
  Nme=names(Reference.report)[p]
  tiff(file=paste("Sensitivity.N0_",Nme,".tiff",sep=""),width=2400,height=2000,units="px",res=300,compression="lzw")
  par(mfcol=c(3,2),las=1,omi=c(.25,.35,.1,.05),mai=c(.15,.3,.15,.2),mgp=c(1.5,0.8,0))
  
  #catch
  fun.plot.Sens.Ref(Reference.report[[p]],1.25)
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"Catch Min",CL.ppr[2],LTY.ppr[2])
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"Catch Max",CL.ppr[3],LTY.ppr[3])
  mtext("Total catch",3)
  legend("bottomleft",c("Minimum","Average","Maximum"),lty=c(LTY.ppr[c(2,1,3)]),lwd=3,col=CL.ppr[c(2,1,3)],bty='n',cex=1.5)
  

  #Selectivity
  fun.plot.Sens.Ref(Reference.report[[p]],1.25)
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"Selectivity Min",CL.ppr[2],LTY.ppr[2])
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"Selectivity Max",CL.ppr[3],LTY.ppr[3])
  mtext("Selectivity",3)
  
  #Max Age
  fun.plot.Sens.Ref(Reference.report[[p]],1.25)
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"Max.Age Min",CL.ppr[2],LTY.ppr[2])
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"Max.Age Max",CL.ppr[3],LTY.ppr[3])
  axis(1,at=AT,labels=Yr.proj[AT],tck=-0.06,cex.axis=1.25)
  mtext("Maximum age",3)
  
  #Rep cycle
  fun.plot.Sens.Ref(Reference.report[[p]],1.25)
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"Rep.cycl Min",CL.ppr[2],LTY.ppr[2])
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"Rep.cycl Max",CL.ppr[3],LTY.ppr[3])
  mtext("Reproductive cycle length",3)
  
  #Pups
  fun.plot.Sens.Ref(Reference.report[[p]],1.25)
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"Pups Min",CL.ppr[2],LTY.ppr[2])
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"Pups Max",CL.ppr[3],LTY.ppr[3])
  mtext("Number of pups",3)
  
  #Natural mortality
  fun.plot.Sens.Ref(Reference.report[[p]],1.25)
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"M Min",CL.ppr[2],LTY.ppr[2])
  fun.plot.Sens(STORE.SENS.reprt[[p]]$"M Max",CL.ppr[3],LTY.ppr[3])
  axis(1,at=AT,labels=Yr.proj[AT],tck=-0.06,cex.axis=1.25)
  mtext("Natural mortality",3)
  
  mtext("Number of females",2,line=1,outer=T,cex=1.25,las=3)
  mtext("Financial year",1,line=0.75,outer=T,cex=1.25)
  dev.off()
}

  #display min, average and max values used in sensitivity analysis
#Catch
tiff(file=paste("Sensitivity.inputs.tiff",sep=""),width=2000,height=2400,units="px",res=300,compression="lzw")
par(mfcol=c(3,1),las=1,mai=c(.45,.5,.075,.15),omi=c(.1,.4,.1,.15),mgp=c(1.5,0.8,0))
show.sens.par(These.Yrs,STORE.SENS.inputs$"3000"$"Catch Min"$"Tot.Ktch",
              Reference.inputs$"3000"$"Tot.Ktch",
              STORE.SENS.inputs$"3000"$"Catch Max"$"Tot.Ktch")
axis(1,at=AT,labels=Yr.proj[AT],tck=-0.06,cex.axis=1.5)
mtext("Total catch",2,las=3,line=4,cex=1.5)
mtext("Financial year",1,line=2.5,cex=1.5)
legend("topleft",c("Minimum","Average","Maximum"),lty=c(LTY.ppr[c(3,1,2)]),lwd=3,col=CL.ppr[c(3,1,2)],bty='n',
       cex=2)

#Selectivity
show.sens.par(0:90,STORE.SENS.inputs$"3000"$"Selectivity Min"$"Selectivity",
              Reference.inputs$"3000"$"Selectivity",
              STORE.SENS.inputs$"3000"$"Selectivity Max"$"Selectivity")
#axis(1,at=seq(0,90,10),labels=seq(0,90,10),tck=-0.06,cex.axis=1.5)
mtext("Selectivity",2,las=3,line=4,cex=1.5)

#M
show.sens.par(0:90,STORE.SENS.inputs$"3000"$"M Min"$"M",
              Reference.inputs$"3000"$"M",
              STORE.SENS.inputs$"3000"$"M Max"$"M")
axis(1,at=seq(0,90,10),labels=seq(0,90,10),tck=-0.06,cex.axis=1.5)
mtext(expression("Natural mortality  " (year^-1)),2,las=3,line=4,cex=1.5)
mtext("Age",1,line=2.5,cex=1.5)
dev.off()


#Maximum age
Max_age_min=STORE.SENS.inputs$"3000"$"Max.Age Min"$"Age_max"
Max_age_avg=Reference.inputs$"3000"$"Age_max"
Max_age_max=STORE.SENS.inputs$"3000"$"Max.Age Max"$"Age_max"

#Reproductive cycle
Cycle_min=STORE.SENS.inputs$"3000"$"Rep.cycl Min"$"Cycle"
Cycle_avg=Reference.inputs$"3000"$"Cycle"
Cycle_max=STORE.SENS.inputs$"3000"$"Rep.cycl Max"$"Cycle"

#Fecundity
Fec_min=STORE.SENS.inputs$"3000"$"Pups Min"$"Fec"
Fec_avg=Reference.inputs$"3000"$"Fec"
Fec_max=STORE.SENS.inputs$"3000"$"Pups Max"$"Fec"

Out.sens.inputs=data.frame(Max.age=c(Max_age_min,Max_age_avg,Max_age_max),
           Cycle=c(Cycle_min,Cycle_avg,Cycle_max),
           Num.pups=c(Fec_min,Fec_avg,Fec_max))
row.names(Out.sens.inputs)=c("Min","Average","Max")
write.csv(Out.sens.inputs,"Out.sens.inputs.csv")



#---10.  FIGURES FOR STEVE SECTION---
#Dropline catch
tiff(file="Steve.Recons.Dropline.CI.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
YMAX=120
par(mfcol=c(1,1),las=1,mai=c(.1,.1,.075,.15),omi=c(.9,.9,.1,.05),mgp=c(1.5,1.5,0))
for (j in 1:1)
{
  fn.plot.ktch(Dropline.CATCH[[j]],XLIM=c(1988.5,2000.5))
  axis(1,at=seq(1988,2001,2),labels=c("1988/89","1990/91","1992/93","1994/95",
                                      "1996/97","1998/99","2000/01"),tck=-0.04,cex.axis=1.5)
  axis(2,at=c(seq(0,YMAX,20)),labels=c(seq(0,YMAX,20)),tck=-0.04,cex.axis=1.5)
}
mtext("Estimated white shark catch (in numbers)",2,3,outer=T,cex=1.75,las=3)
mtext("Financial year",1,2.5,outer=T,cex=1.75)
dev.off()

#Barplots of cpue
cpue.Folly=read.csv(handl_OneDrive("Data/Population dynamics/CPUE.Folly.csv"))
cpue.Folly2=read.csv(handl_OneDrive("Data/Population dynamics/CPUE.Folly.Method2.csv"))
tiff(file="Steve.CPUEs.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
par(mfcol=c(2,1),las=1,mai=c(.1,.1,.075,.15),omi=c(.8,.8,.1,.05),mgp=c(1.5,.8,0))
fn.barplot(cpue.Folly,"NO")
legend("top",c("WC","Z1","Z2"),bty='n',fill=c("grey80","white","grey50"),cex=1.5)
fn.barplot(cpue.Folly2,"YES")
mtext("Number of white sharks per 1000 km gn.d",2,2,outer=T,cex=1.75,las=3)
mtext("Time period",1,2.5,outer=T,cex=1.75)
dev.off()



#---11. FUTURE PROJECTIONS STEVE SECTION---
#NOTE: Check if population recovers in the future
Future.yrs=100
Check.future="NO"
if(Check.future=="YES")
{
  Ns=100
  iterations=1:Ns
  N.mdld=N.mdld+Future.yrs
  LH1.N1.Sel1.Us.ab.future=vector('list',length=length(PCS.scenarios))
  names(LH1.N1.Sel1.Us.ab.future)=PCS.scenarios
  LH2.N1.Sel2.Us.ab.future=LH1.N1.Sel1.Us.ab.future
  
  system.time(for(pcs in 1:length(PCS.scenarios))
  {
    TOTAL.CATCH=TOTAL.CATCH.PCS[[pcs]]  
    STore.sim1=vector('list',length(TOTAL.CATCH))
    names(STore.sim1)=names(TOTAL.CATCH)
    STore.sim2=STore.sim1
    
    for(i in 1:length(TOTAL.CATCH))
    {      
      STR.sim1=vector('list',length(N.0))
      names(STR.sim1)=N.0
      STR.sim2=STR.sim1
      
      for (no in 1:length(N.0))
      {
        #Sample initial population size 
        N.0.sim=round(runif(Ns,N.0[[no]][1],N.0[[no]][2]))
        
        CTCH=cbind(TOTAL.CATCH[[i]],matrix(rep(0,Future.yrs*nrow(TOTAL.CATCH[[i]])),ncol=Future.yrs))
        n.Yr.tot=ncol(CTCH)
        
        #Run risk analysis
        STR.sim1[[no]]=Risk.fun(Life.hist.scenarios[[1]],CTCH,N.0.samp=N.0.sim)
        STR.sim2[[no]]=Risk.fun(Life.hist.scenarios[[2]],CTCH,N.0.samp=N.0.sim)  
      }
      #Store outputs
      STore.sim1[[i]]=STR.sim1       
      STore.sim2[[i]]=STR.sim2
    }
    LH1.N1.Sel1.Us.ab.future[[pcs]]=STore.sim1       
    LH2.N1.Sel2.Us.ab.future[[pcs]]=STore.sim2
    
  })       #(takes 59 seconds per iteration)  
  
  n.Yr.prol=n.Yr.tot
  Yr.proj=c(Yr.proj,seq((Yr.proj[length(Yr.proj)]+1),(Yr.proj[length(Yr.proj)]+Future.yrs)))
  These.Yrs=c(These.Yrs,seq((These.Yrs[length(These.Yrs)]+1),(These.Yrs[length(These.Yrs)]+Future.yrs)))
  
  Yrs.future=Yrs.future=((n.Yr.prol-Future.yrs+1):n.Yr.prol)
  fn.future.pol=function(Yrs) polygon(c(Yrs[1],Yrs[length(Yrs)],Yrs[length(Yrs)],Yrs[1]),
                             c(-0.5,-0.5,1.5,1.5),col=rgb(.1,.6,.3,alpha=.2))
  
  for(x in 1:length(PCS.scenarios))
  {
    tiff(file=paste("Future_projections.",names(Store.Paper.median.all)[x],".tiff",sep=""),
         width = 2600, height = 2000,units = "px", res = 300, compression = "lzw")   
    
    par(mfcol=c(5,4),las=1,mai=c(.1,.1,.1,0.01),omi=c(.45,.8,.2,.01),mgp=c(1.5,0.6,0),xpd=F)
    STR.LH.y=vector('list',length(KTCH.scen))
    names(STR.LH.y)=KTCH.scen
    
    for(y in 1:length(KTCH.scen))  
    {
      STR.LH.z=vector('list',length(N.0))
      names(STR.LH.z)=N.0
      for (z in 1:length(N.0))   
      {
        STR.LH.z[[z]]=plot.mean.CI.Paper(LH1.N1.Sel1.Us.ab.future[[x]][[y]][[z]]$Rel.abundance,
                                         LH2.N1.Sel2.Us.ab.future[[x]][[y]][[z]]$Rel.abundance,CEXX=0.9)
        if(z==1)mtext(KTCH.scen[y],side=3,line=0,font=1,las=0,cex=1.75)      
        if(y==1)mtext(paste("[",N.0[z][[1]][1],";",N.0[z][[1]][2],"]",sep=""),side=2,line=3,font=1,las=0,cex=1.15)
         fn.future.pol(Yrs.future) 
      }
      
    }
    mtext("Financial year",side=1,outer=T,line=1.95,font=1,las=0,cex=1.75)
    mtext("Relative number of females",side=2,outer=T,line=4,font=1,las=0,cex=1.75)
    dev.off()
  }
  
}




##################################################################################

#---12. PREVIOUS CODE SECTION---
USED.PREVIOUSY="YES"
if(USED.PREVIOUSY=="YES")
{
  #Figure 2.2 Exploitation rates                       
  
  tiff(file="Figure2.2.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  par(mfcol=c(2,1),las=1,omi=c(.65,.85,.1,.01),mai=c(.15,.3,.1,.25),mgp=c(1.5,0.8,0),xpd=F)
  
  plot.mean.CI(LH1.N1.Sel1.Us.ab$F_series,1,leg="YES",RELATIVE="NO",c(0,0.07),PERCENTILES="YES")
  legend("topright","LH1",bty='n',cex=2)  
  plot.mean.CI(LH2.N1.Sel2.Us.ab$F_series,1,leg="NO",RELATIVE="NO",c(0,0.07),PERCENTILES="YES")
  legend("topright","LH2",bty='n',cex=2)
  axis(1,at=seq(match(1940,Yr.proj),n.Yr.prol,by=20),
       labels=seq(Yr.proj[match(1940,Yr.proj)],Yr.proj[length(Yr.proj)],by=20),
       cex.axis=1.5,las=1,tck=-0.06)
  mtext("Year",side=1,outer=T,line=1.5,font=1,las=0,cex=2)
  mtext("Exploitation rate",side=2,outer=T,line=2.5,font=1,las=0,cex=2)
  dev.off()
  
  
  #calculate probabilities of population doubling from 1998 to Yr.end (for Yr.end estimates)
  #note: this shows prob of doubling since protection (i.e. nyears)
  cum.dist.t2.1=ecdf(este.t2.1)
  cum.dist.t2=ecdf(este.t2)
  nyears=length(1998:Yr.end)
  Prob.1998_2012.t2.1=cum.dist.t2.1(nyears)  #scenario 1
  Prob.1998_2012.t2=cum.dist.t2(nyears)      #scenario 2
  
  
  #Mature biomass in current year    NEW
  Nmat.LH1=matrix(NA,nrow=NN,ncol=3)
  Nmat.LH1=as.data.frame(Nmat.LH1)
  colnames(Nmat.LH1)=c("Low95CI","Median","Up95CI")
  Nmat.LH2=Nmat.LH1
  tiff(file="Nmature.total.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  par(mfcol=c(5,2),las=1,omi=c(.4,.5,.1,.01),mai=c(.15,.3,.1,.25),mgp=c(1.5,0.8,0),xpd=F)
  for(i in 1:NN)
  {
    aa=LH1.N1.Sel1.Us.ab[[i]]$Mature.rel.abundance
    bb=LH2.N1.Sel2.Us.ab[[i]]$Mature.rel.abundance
    Nmat=ncol(aa)
    aa=aa[,Nmat]
    bb=bb[,Nmat]
    
    Nmat.LH1[i,]=quantile(aa,probs=c(.025,.5,.975))
    Nmat.LH2[i,]=quantile(bb,probs=c(.025,.5,.975))
    
    plot(density(bb, adjust=2,na.rm =T), type="l",lty=1,col=greyscale[1],xlab="",lwd=LWD,
         ylab="",main="",yaxt="n",xlim=c(0,max(aa*1.15)),cex.axis=1.3,cex.lab=1.65)
    lines(density(aa, adjust=2,na.rm =T), lty=1,col=greyscale[1],lwd=LWD)
    lines(density(bb, adjust=2,na.rm =T), lty=1,col=greyscale[2],lwd=LWD)
    if(i==6)legend('topright',c('LH1','LH2'),lty=1,col=greyscale,bty='n',cex=1.25,lwd=LWD)
    
  }
  mtext("Distribution",side=2,line=1,outer=T,font=1,las=0,cex=1.75)
  mtext(expression('N'[mature(Yr.end)]),side=1,line=2.2,font=1,outer=T,las=0,cex=1.75)
  dev.off()
  
  write.table(Nmat.LH1,"Nmat.Yr.end.LH1.csv",sep=",",row.names=F)
  write.table(Nmat.LH2,"Nmat.Yr.end.LH2.csv",sep=",",row.names=F)
  
  
  
  # Figure 4 for Rick, combining 0% and 50% survival
  #NO.PCS1=LH1.Sel1.ab  #Store here the run for 0% survival
  #NO.PCS2=LH2.Sel2.ab
  
  
  # tiff(file="Figure4.total.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  # 
  # par(mfcol=c(4,2),las=1,omi=c(.6,.7,.1,.1),mai=c(.15,.3,.15,.25),mgp=c(1.5,0.8,0),xpd=F)
  # for(i in 1:NN) 
  # {
  #   plot.mean.CI(NO.PCS1[[i]]$Rel.abundance,NO.PCS2[[i]]$Rel.abundance)
  #   if(i==1) mtext("0% survival",3)
  #   if(i ==4)
  #   {
  #     axis(1,at=seq(match(1940,Yr.proj),n.Yr.prol,by=20),
  #          labels=seq(Yr.proj[match(1940,Yr.proj)],Yr.proj[length(Yr.proj)],by=20),
  #          cex.axis=1.5,las=1,tck=-0.06)
  #   }
  #   
  # }
  # for(i in 1:NN) 
  # {
  #   plot.mean.CI(LH1.Sel1.ab[[i]]$Rel.abundance,LH2.Sel2.ab[[i]]$Rel.abundance)
  #   if(i==1) mtext("50% survival",3)  
  #   if(i ==4)
  #   {
  #     axis(1,at=seq(match(1940,Yr.proj),n.Yr.prol,by=20),
  #          labels=seq(Yr.proj[match(1940,Yr.proj)],Yr.proj[length(Yr.proj)],by=20),
  #          cex.axis=1.5,las=1,tck=-0.06)
  #   }
  #   
  # }
  # mtext("Year",side=1,outer=T,line=1.5,font=1,las=0,cex=1.6)
  # mtext("Total number of females",side=2,outer=T,line=1.8,font=1,las=0,cex=1.6)
  # 
  # dev.off()
  
  
  
  
  
  #4.2. Semi-log plot
  # handle=function()
  # {
  #   plot(1:n.Yr.prol,seq(MIN.Y,MAX.Y,length.out=n.Yr.prol),col='transparent',xaxt='n',yaxt='n',xlab="",
  #        ylab="",xlim=c(1,(n.Yr.prol+4)))
  #   axis(1,at=1:n.Yr.prol,labels=F,cex.axis=1.35,las=1,tck=-0.03)
  # }
  # 
  # 
  # line.type=1:NN
  # 
  # 
  # plot.what=function(DATA,line.type1,line.col1) 
  # {
  #   MED=log(colMedians(DATA))[These.Yrs]  #select 1956 to Yr.end
  #   MED=ifelse(MED==-Inf,0,MED)
  #   lines(1:n.Yr.prol,MED,type='l',lty=line.type1,col=line.col1,lwd=LWD)
  # }  
  # 
  
  # 
  # poly.prot=function(MAX.Y)polygon(x=c(id.1998,id.2012,id.2012,id.1998),y=c(0,0,MAX.Y,MAX.Y),lwd=2,
  #                                  col=rgb(.2, .2, .2,0.15), border=NA)
  # 
  # #Pop.size.leg=c('N1','N2','N3')
  
  # N.N=N.0 
  # THIS.F=N.N[c(1,4,5,6)]
  # #THIS.F=N.N[c(1,9,19,20)]
  # 
  # AXIS1=function()axis(1,at=seq(1,n.Yr.prol,by=5),labels=seq(Yr.proj[1],Yr.proj[length(Yr.proj)],by=5),cex.axis=1.1,las=1,tck=-0.06)
  # AXIS2=function()axis(2,at=log(THIS.F),labels=THIS.F/1000,cex.axis=1.35,las=1,tck=-0.035)
  # AXIS3=function()axis(1,at=seq(1,n.Yr.prol,by=5),labels=F,cex.axis=1.1,las=1,tck=-0.06)
  # AXIS4=function()axis(2,at=log(N.N),labels=F,cex.axis=1.1,las=1,tck=-0.02)
  # ab.fn=function(n) lines(1:n.Yr.prol,log(rep(n,n.Yr.prol)),lty=5,col='grey',lwd=2.5)
  # ab.text=function(n)text(1,log(n*1.5),n,cex=1.35,col='grey40')
  
  
  # tiff(file="Figure4.log.total.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  #   par(mfcol=c(2,1),las=1,omi=c(.7,1,.1,.3),mai=c(.15,.01,.15,.175),mgp=c(1.5,0.8,0),xpd=F)
  # 
  #   
  #   #LH1, sel 2
  #   MAX.Y=log(max(N.N)*1.15)
  #   MIN.Y=0.08
  #   handle();for(i in 1:NN) plot.what(LH1.Sel1.ab[[i]]$Rel.abundance,1,line.col[i])
  #   ab.fn(100)
  #   AXIS2(); AXIS3();AXIS4()
  #   mtext("LH1",side=4,line=0.4,font=1,las=1,cex=1.3,col='black')
  #   legend("bottomright",F_names[This.F],bty='n',lty=1,col=line.col[This.F],lwd=LWD,cex=0.8)
  #   
  #   if(plot.power.point=="YES")
  #   {
  #     fun.back()
  #     ab.fn(100)
  #     ab.text(100)
  # 
  #   }
  # 
  #   #LH2, sel 2
  #   handle();for(i in 1:NN) plot.what(LH2.Sel2.ab[[i]]$Rel.abundance,1,line.col[i])
  #   AXIS2();AXIS1();AXIS4()
  #   ab.fn(100)
  #   mtext("LH2",side=4,line=0.4,font=1,las=1,cex=1.3,col='black')
  #   if(plot.power.point=="YES")
  #   {
  #     fun.back()
  #     ab.fn(100)
  #     ab.text(100)
  # 
  #   }
  #   mtext("Year",side=1,outer=T,line=1.5,font=1,las=0,cex=1.6)
  #   mtext("Total number of females (in thousands)",side=2,outer=T,line=3,font=1,las=0,cex=1.6)
  #   legend("topright",F_names[This.F2],bty='n',lty=1,col=line.col[This.F2],lwd=LWD,cex=0.8,horiz = F,
  #          inset = 0,box.lwd = 0,box.col = "white",bg = "white")
  #   
  # dev.off()
  
  
  # Mature biomass
  #Mature biomass
  # if(length(N1975.scenarios)==1)
  # {
  #   tiff(file="Figure4.log.mature.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  #   par(mfcol=c(2,1),las=1,omi=c(.7,1,.1,.3),mai=c(.15,.01,.15,.175),mgp=c(1.5,0.8,0),xpd=F)
  #   N.N=N.mat
  #   
  #   #LH1, sel 2
  #   MAX.Y=log(max(N.N)*1.15)
  #   MIN.Y=0.08
  #   handle();for(i in 1:NN) plot.what(LH1.Sel1.ab[[i]]$Mature.rel.abundance,1,line.col[i])
  #   ab.fn(100)
  #   axis(2,at=log(THIS.F/10),labels=(THIS.F/10)/1000,cex.axis=1.35,las=1,tck=-0.035)
  #   AXIS3();AXIS4()
  #   mtext("LH1",side=4,line=0.4,font=1,las=1,cex=1.3,col='black')
  #   legend("bottomright",F_names[This.F],bty='n',lty=1,col=line.col[This.F],lwd=LWD,cex=0.8)
  #   if(plot.power.point=="YES") fun.back()
  #   
  #   
  #   #LH2, sel 2
  #   handle();for(i in 1:NN) plot.what(LH2.Sel2.ab[[i]]$Mature.rel.abundance,1,line.col[i])
  #   axis(2,at=log(THIS.F/10),labels=(THIS.F/10)/1000,cex.axis=1.35,las=1,tck=-0.035)
  #   AXIS1();AXIS4()
  #   ab.fn(100)
  #   mtext("LH2",side=4,line=0.4,font=1,las=1,cex=1.3,col='black')
  #   if(plot.power.point=="YES") fun.back()
  #   legend("topright",F_names[This.F2],bty='n',lty=1,col=line.col[This.F2],lwd=LWD,cex=0.8,horiz = F,
  #          inset = 0,box.lwd = 0,box.col = "white",bg = "white")
  #   
  #   
  #   mtext("Year",side=1,outer=T,line=1.5,font=1,las=0,cex=1.6)
  #   mtext("Total number of females (in thousands)",side=2,outer=T,line=3,font=1,las=0,cex=1.6)
  #   dev.off()
  # }
  
  
  # plausible=c("F1","F2","F3","F4","F5")
  # ID.plaus=match(plausible,names(LH1.Sel1.ab[[1]]))
  # Fishing.leg.pl=names(LH2.Sel2.ab.plaus[[3]])
  
  
  
  #Outputs for risk assessment
  
  #Plot all runs 
  #Total biomass
  # tiff(file="N.traj.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  # par(mfcol=c(5,4),las=1,omi=c(.7,1,.1,.1),mai=c(.15,.3,.15,.25),mgp=c(1.5,0.8,0),xpd=F)
  # for (a in 1:length(N.0))
  # {
  #   plot(LH1.N1.Sel1.Us.ab[[a]][[1]][1,],ylim=c(0,N.0[a]*1.1),main=N.0[a],ylab="",xlab="")
  #   for (i in 1:Ns)
  #   {
  #     lines(LH1.N1.Sel1.Us.ab[[a]]$Rel.abundance[i,])
  #     lines(LH2.N1.Sel2.Us.ab[[a]]$Rel.abundance[i,],col=2)
  #     legend("bottomleft",c("LH1","LH2"),bty="n",lty=1,col=c(1,2))
  #   }
  #   
  # }
  # dev.off()
  
  
  
  #Plot mean and confidence intervals 
  #Total biomass
  plot.mean.CI=function(DATA) 
  {
    DATA=DATA[,These.Yrs]
    MED=colMedians(DATA)
    LOW=apply(DATA, 2, function(x) quantile(x, 0.025))
    UP=apply(DATA, 2, function(x) quantile(x, 0.975))
    
    plot(1:n.Yr.prol,MED,type='l',lty=1,col=1,lwd=LWD,ylim=c(0,MED[1]*1.25),xaxt='n',xlab="",ylab="",
         cex.axis=0.95)
    lines(1:n.Yr.prol,LOW,type='l',lty=1,col=2,lwd=LWD)
    lines(1:n.Yr.prol,UP,type='l',lty=1,col=2,lwd=LWD)
    axis(1,at=seq(1,n.Yr.prol,by=1),labels=F,las=1,tck=-0.02)
    axis(1,at=seq(1,n.Yr.prol,by=5),labels=seq(Yr.proj[1],Yr.proj[length(Yr.proj)],by=5),cex.axis=0.95,las=1,tck=-0.04)
  }  
  tiff(file="For Rick/N.traj.CI.LH1.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  par(mfcol=c(5,2),las=1,omi=c(.7,1,.1,.1),mai=c(.15,.3,.15,.25),mgp=c(1.5,0.8,0),xpd=F)
  for(i in 1:NN) 
  {
    plot.mean.CI(LH1.Sel1.ab[[i]]$Rel.abundance)
    legend("bottomleft",paste("No=",N.0[i]),bty='n')
  }
  mtext("Year",side=1,outer=T,line=1.5,font=1,las=0,cex=1.6)
  mtext("Total number of females",side=2,outer=T,line=1.8,font=1,las=0,cex=1.6)
  dev.off()
  
  tiff(file="For Rick/N.traj.CI.LH2.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  par(mfcol=c(5,2),las=1,omi=c(.7,1,.1,.1),mai=c(.15,.3,.15,.25),mgp=c(1.5,0.8,0),xpd=F)
  for(i in 1:NN) 
  {
    plot.mean.CI(LH2.Sel2.ab[[i]]$Rel.abundance)
    legend("bottomleft",paste("No=",N.0[i]),bty='n')
  }
  mtext("Year",side=1,outer=T,line=1.5,font=1,las=0,cex=1.6)
  mtext("Total number of females",side=2,outer=T,line=1.8,font=1,las=0,cex=1.6)
  dev.off()
  
  
  
  
  
  #get values
  #Total biomass
  n.pop=length(N1975.scenarios)
  
  Dep.level.LH1.N1=vector('list',NN)
  Dep.level.LH1.N2=Dep.level.LH1.N3=Dep.level.LH2.N1=Dep.level.LH2.N2=Dep.level.LH2.N3=Dep.level.LH1.N1
  N2012.LH1.N1=N2012.LH1.N2=N2012.LH1.N3=N2012.LH2.N1=N2012.LH2.N2=N2012.LH2.N3=Dep.level.LH1.N1
  
  
  get.what=function(DATA,scaler) 
  {
    MED=colMedians(DATA)
    Depletion=MED[id.2012]/scaler
    return(list(Depletion=Depletion,N2012=MED[id.2012]))
  }  
  
  
  for(i in 1:NN)
  {
    Dep.level.LH1.N1[[i]]=get.what(LH1.Sel1.ab[[i]]$Rel.abundance,N.0[i])$Depletion   
    Dep.level.LH2.N1[[i]]=get.what(LH2.Sel2.ab[[i]]$Rel.abundance,N.0[i])$Depletion  
    
    
    N2012.LH1.N1[[i]]=get.what(LH1.Sel1.ab[[i]]$Rel.abundance,N.0[i])$N2012 
    N2012.LH2.N1[[i]]=get.what(LH2.Sel2.ab[[i]]$Rel.abundance,N.0[i])$N2012  
    
  }
  Depletion.LH1=round(matrix(unlist(Dep.level.LH1.N1),ncol=NN,nrow=n.pop,byrow=T),3)
  Depletion.LH2=round(matrix(unlist(Dep.level.LH2.N1),ncol=NN,nrow=n.pop,byrow=T),3)
  
  N2012.LH1=round(matrix(unlist(N2012.LH1.N1),ncol=NN,nrow=n.pop,byrow=T),3)
  N2012.LH2=round(matrix(unlist(N2012.LH2.N1),ncol=NN,nrow=n.pop,byrow=T),3)
  
  
  
  colnames(Depletion.LH1)=colnames(Depletion.LH2)=colnames(N2012.LH1)=colnames(N2012.LH2)=N.0
  
  write.table(Depletion.LH1,"Figure4.Depletion.LH1.total.csv",sep=",",row.names=F)
  write.table(Depletion.LH2,"Figure4.Depletion.total.LH2.csv",sep=",",row.names=F)
  write.table(N2012.LH1,"Figure4.N2012.LH1.total.csv",sep=",",row.names=F)
  write.table(N2012.LH2,"Figure4.N2012.LH2.total.csv",sep=",",row.names=F)
  
  
  
  
  
  #Figure 5 Probability of population increase 
  
  fn.prob.inc=function(DATA,factor)
  {
    a=DATA[,id.2012]
    b=DATA[,id.1998]
    Prob.Pop.increase=sum(a>b*factor)/nrow(DATA)
  }
  
  
  #Total biomass
  LH1.N1.Sel1.Us=LH2.N1.Sel2.Us=LH1.N1.Sel1.Us.1=LH2.N1.Sel2.Us.1=
    LH1.N1.Sel1.Us.2=LH2.N1.Sel2.Us.2=LH1.N1.Sel1.Us.3=LH2.N1.Sel2.Us.3=rep(NA,NN)
  for (i in 1:NN)
  {
    LH1.N1.Sel1.Us[i]=fn.prob.inc(LH1.N1.Sel1.Us.ab[[i]]$Rel.abundance,1) 
    LH2.N1.Sel2.Us[i]=fn.prob.inc(LH2.N1.Sel2.Us.ab[[i]]$Rel.abundance,1)
    LH1.N1.Sel1.Us.1[i]=fn.prob.inc(LH1.N1.Sel1.Us.ab[[i]]$Rel.abundance,1.01) 
    LH2.N1.Sel2.Us.1[i]=fn.prob.inc(LH2.N1.Sel2.Us.ab[[i]]$Rel.abundance,1.01)
    LH1.N1.Sel1.Us.2[i]=fn.prob.inc(LH1.N1.Sel1.Us.ab[[i]]$Rel.abundance,1.02) 
    LH2.N1.Sel2.Us.2[i]=fn.prob.inc(LH2.N1.Sel2.Us.ab[[i]]$Rel.abundance,1.02)
    LH1.N1.Sel1.Us.3[i]=fn.prob.inc(LH1.N1.Sel1.Us.ab[[i]]$Rel.abundance,1.03) 
    LH2.N1.Sel2.Us.3[i]=fn.prob.inc(LH2.N1.Sel2.Us.ab[[i]]$Rel.abundance,1.03)
  }
  
  
  DIMNAMES=list(N.0,names(N1975.scenarios))
  
  LH1.Sel1=matrix(LH1.N1.Sel1.Us,ncol=n.pop,nrow=NN,byrow=F,dimnames=DIMNAMES)
  LH2.Sel2=matrix(LH2.N1.Sel2.Us,ncol=n.pop,nrow=NN,byrow=F,dimnames=DIMNAMES)
  LH1.Sel1.1=matrix(LH1.N1.Sel1.Us.1,ncol=n.pop,nrow=NN,byrow=F,dimnames=DIMNAMES)
  LH2.Sel2.1=matrix(LH2.N1.Sel2.Us.1,ncol=n.pop,nrow=NN,byrow=F,dimnames=DIMNAMES)
  LH1.Sel1.2=matrix(LH1.N1.Sel1.Us.2,ncol=n.pop,nrow=NN,byrow=F,dimnames=DIMNAMES)
  LH2.Sel2.2=matrix(LH2.N1.Sel2.Us.2,ncol=n.pop,nrow=NN,byrow=F,dimnames=DIMNAMES)
  LH1.Sel1.3=matrix(LH1.N1.Sel1.Us.3,ncol=n.pop,nrow=NN,byrow=F,dimnames=DIMNAMES)
  LH2.Sel2.3=matrix(LH2.N1.Sel2.Us.3,ncol=n.pop,nrow=NN,byrow=F,dimnames=DIMNAMES)
  
  plausible=rownames(LH2.Sel2)
  N.plaus=NN
  
  
  #create plotting functions
  handle=function()
  {
    plot(1:N.plaus,seq(0,1,length.out=N.plaus),col='transparent',xaxt='n',yaxt='n',xlab="",ylab="",
         panel.first =abline(h = 0.5, lty = 5, lwd=1.5,col = 'grey'))
    axis(1,at=1:N.plaus,labels=F,cex.axis=1.35,las=1,tck=-0.03)  
  }
  if(plot.power.point=="NO")
  {
    CLS.prob=c("black","grey30","grey50","grey70")
    ln.type=1:4
  }
  
  if(plot.power.point=="YES")
  {
    CLS.prob=c("darkgreen","red","blue","orange")
    ln.type=rep(1,4)
  }
  
  plot.what=function(data,ln,lncol) 
  {
    if(Drop=="LH2.F2")data[1,]=c(NA,NA,NA)
    for(i in 1:ncol(data))lines(1:N.plaus,data[,i],type='o',lty=ln,col=lncol,lwd=LWD)
  }  
  
  #LEG.LAB=paste("(",plausible,")",sep="")
  if(length(N1975.scenarios)==1)
  {
    LEG.LAB2=N.0/1000
    #   LEG.LAB2=c(expression(paste("(",N[1975],"=",2250,")",sep="")),
    #              expression(paste("(",N[1975],"=",4500,")",sep="")),
    #              expression(paste("(",N[1975],"=",9000,")",sep="")),
    #              expression(paste("(",N[1975],"=",18000,")",sep="")))
  }
  
  
  tiff(file="Figure5.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  par(mfcol=c(2,1),las=1,omi=c(.9,.8,.1,.1),mai=c(.15,.2,.15,.3),mgp=c(1.5,0.8,0),xpd=F)
  
  Drop="NONE"
  handle()   
  plot.what(LH1.Sel1,ln.type[1],CLS.prob[1])
  plot.what(LH1.Sel1.1,ln.type[2],CLS.prob[2])
  plot.what(LH1.Sel1.2,ln.type[3],CLS.prob[3])
  plot.what(LH1.Sel1.3,ln.type[4],CLS.prob[4])
  if(plot.power.point=="YES") fun.back()
  axis(2,at=seq(0,1,by=.2),labels=seq(0,1,by=.2),cex.axis=1.3,las=1,tck=-0.03)
  legend('topright',c("LH1"),bty='n',cex=1.5)
  
  #Drop="LH2.F2"  #this drops F2 for LH2
  Drop="NONE"
  handle()
  plot.what(LH2.Sel2,ln.type[1],CLS.prob[1])
  plot.what(LH2.Sel2.1,ln.type[2],CLS.prob[2])
  plot.what(LH2.Sel2.2,ln.type[3],CLS.prob[3])
  plot.what(LH2.Sel2.3,ln.type[4],CLS.prob[4])
  axis(2,at=seq(0,1,by=.2),labels=seq(0,1,by=.2),cex.axis=1.3,las=1,tck=-0.03)
  axis(1,at=1:N.plaus,labels=LEG.LAB2,cex.axis=1.3,las=1,tck=-0.03)
  # axis(1, at=1:N.plaus, labels = LEG.LAB, line=1.05,cex.axis=0.625,col="transparent") #second axis 
  legend('topright',c("LH2"),bty='n',cex=1.5)
  legend("topleft",c("Increase","1% increase","2% increase","3% increase"),
         bty='n',lty=ln.type,col=CLS.prob,lwd=LWD,cex=1.25)
  if(plot.power.point=="YES") fun.back()
  #mtext(expression("Fishing mortality  " (year^-1)),side=1,outer=T,line=3.5,font=1,las=0,cex=1.75)
  mtext(expression(paste(N[0]," (in thousands of females)")),side=1,outer=T,line=2,font=1,las=0,cex=1.75)
  mtext("Probability of population increase",side=2,outer=T,line=2,font=1,las=0,cex=1.75)
  dev.off()
  
  
  write.table(LH1.Sel1,"Figure5.LH1.probs.csv",sep=",",row.names=T)
  write.table(LH2.Sel2,"Figure5.LH2.probs.csv",sep=",",row.names=T)
  
  
  #Figure 6 Probability of recovery for plausible scenarios
  
  #calculate prob of N2012 > Biom.ref.point
  dummy=vector(length=Ns)
  probs=vector(length=NN)
  probs1=vector('list',length=n.pop)
  
  fn.prob.RP=function(DATA,Ref.point,scaler)
  {
    for (b in 1:NN)
    {   
      n2012=DATA[[b]][[1]][,id.2012]
      probs[b]=sum(n2012>Ref.point*scaler[b])/length(n2012)
    }
    
    return(probs)  
  }
  #Biom.ref.point.1
  one.one=fn.prob.RP(LH1.Sel1.ab,Biom.ref.point.1,N.0)
  two.two=fn.prob.RP(LH2.Sel2.ab,Biom.ref.point.1,N.0)
  
  LH1.Sel1.brp1=one.one
  LH2.Sel2.brp1=two.two
  
  #LH1.Sel1.brp1=t(do.call(rbind,one.one))
  #LH2.Sel2.brp1=t(do.call(rbind,two.two))
  
  #Biom.ref.point.2
  one.one=fn.prob.RP(LH1.Sel1.ab,Biom.ref.point.2,N.0)
  two.two=fn.prob.RP(LH2.Sel2.ab,Biom.ref.point.2,N.0)
  
  LH1.Sel1.brp2=one.one
  LH2.Sel2.brp2=two.two
  
  
  
  #create plotting functions
  handle=function()
  {
    plot(1:N.plaus,seq(0,1,length.out=N.plaus),col='transparent',xaxt='n',yaxt='n',xlab="",ylab="",
         panel.first =abline(h = 0.5, lty = 5, lwd=1.5,col = 'grey'))
    axis(1,at=1:N.plaus,labels=F,cex.axis=1.35,las=1,tck=-0.03)  
  }
  
  
  plot.what=function(data,ln,lncol) 
  {
    if(Drop=="LH2.F2")data[1,]=c(NA,NA,NA)
    lines(1:N.plaus,data,type='o',lty=ln,col=lncol,lwd=LWD)
  }  
  
  #Total biomass
  
  tiff(file="Figure6.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  par(mfcol=c(2,1),las=1,omi=c(.9,.8,.1,.1),mai=c(.15,.2,.15,.3),mgp=c(1.5,0.8,0),xpd=F)
  
  #LH1
  Drop="NONE"
  handle()   
  plot.what(LH1.Sel1.brp1,ln.type[1],CLS.prob[1])
  plot.what(LH1.Sel1.brp2,ln.type[2],CLS.prob[2])
  axis(2,at=seq(0,1,by=.2),labels=seq(0,1,by=.2),cex.axis=1.3,las=1,tck=-0.03)
  legend("bottomleft",c(expression("N"[paste(30,"%")]),expression("N"[paste(60,"%")])),
         bty='n',lty=ln.type[1:2],col=CLS.prob[1:2],lwd=LWD,cex=1.75)
  legend('bottomright',c("LH1"),bty='n',cex=1.5)
  if(plot.power.point=="YES") fun.back()
  
  #LH2
  Drop="NONE"
  handle()   
  plot.what(LH2.Sel2.brp1,ln.type[1],CLS.prob[1])
  plot.what(LH2.Sel2.brp2,ln.type[2],CLS.prob[2])
  axis(2,at=seq(0,1,by=.2),labels=seq(0,1,by=.2),cex.axis=1.3,las=1,tck=-0.03)
  axis(1,at=1:N.plaus,labels=LEG.LAB2,cex.axis=1.3,las=1,tck=-0.03)
  if(plot.power.point=="YES") fun.back()
  legend('bottomright',c("LH2"),bty='n',cex=1.5)
  mtext(expression(paste(N[0]," (in thousands of females)")),side=1,outer=T,line=2,font=1,las=0,cex=1.75)
  mtext("Probability of population recovery",side=2,outer=T,line=2,font=1,las=0,cex=1.75)
  
  dev.off()
  
  
  names(LH1.Sel1.brp1)=names(LH1.Sel1.brp2)=names(LH2.Sel2.brp1)=names(LH2.Sel2.brp2)=N.0
  
  write.table(LH1.Sel1.brp1,"Figure6.LH1.brp1.probs.csv",sep=",",row.names=T)
  write.table(LH1.Sel1.brp2,"Figure6.LH1.brp2.probs.csv",sep=",",row.names=T)
  write.table(LH2.Sel2.brp1,"Figure6.LH2.brp1.probs.csv",sep=",",row.names=T)
  write.table(LH2.Sel2.brp2,"Figure6.LH2.brp2.probs.csv",sep=",",row.names=T)
  
  
  
  #--------------- Forward projections -------------------------  #MISSING: N0 DRAWED FROM UNIFORM DISTRIBUTION, NOT CONSIDERING 8 CATCH SCENARIOS!!!!
  
  ########
  #-- 1. Four year projections for drumline program ---
  #######
  #note: projections into the future to include the capture of sharks in the drumline program
  #      Catch is the assumed drumline catch plus the average catch in commercial fisheries for last 5 years
  #       The way this works is by redoing the simulation with the extra yrs.sim years
  
  
  yrs.sim=4   #number of years projected
  N.0.future=list(N.0.Rick,N.0.Rick)  #N.0 requested by Rick
  names(N.0.future)=c("LH1","LH2")
  NN.future=length(N.0.future[[1]])
  drum.ktch=c(10,20)
  N.drum=length(drum.ktch)
  Store.Tot.Ktch=Tot.Ktch
  yr.dummy=n.Yr.prol
  
  Store.Future=vector('list',N.drum)
  names(Store.Future)=paste("D.L. catch=",drum.ktch)
  
  N.mdld=N.mdld+yrs.sim
  N.modl.yrs=c(N.modl.yrs,2013:2016)
  
  These.Yrs=match(c(Yr.proj,2013:2016),N.modl.yrs)
  Yr.proj=c(Yr.proj,2013:2016)
  
  system.time(for (h in 1:N.drum)
  {
    Tot.Ktch=Store.Tot.Ktch
    Tot.Ktch.future=matrix(rep(rowMeans(Tot.Ktch[,(length(Total.catch)-4):length(Total.catch)]),yrs.sim),ncol=yrs.sim)
    Tot.Ktch.future=round(Tot.Ktch.future+drum.ktch[h])
    Tot.Ktch=cbind(Tot.Ktch,Tot.Ktch.future)
    
    n.Yr.tot=ncol(Tot.Ktch)
    
    
    #rel abundance
    LH1.N1.Sel1.Us.ab=LH2.N1.Sel2.Us.ab=vector('list',length=NN.future)
    names(LH1.N1.Sel1.Us.ab)=N.0.future[[1]]
    names(LH2.N1.Sel2.Us.ab)=N.0.future[[2]]
    
    #Risk analysis
    for (x in 1:length(Life.hist.scenarios))
    {
      scenario=Life.hist.scenarios[[x]]
      for(j in 1:NN.future) 
      {
        if(scenario==1)N.init=N.0.future[[1]][j]
        if(scenario==2)N.init=N.0.future[[2]][j]
        
        #1. Create elements to fill in
        Pop.size.ratio=No.vec=rep(NA,length = Ns)
        Rel.abundance=Mature.rel.abundance=F_vec=matrix(NA,ncol=N.mdld,nrow = Ns)
        
        for (s in iterations)
        {
          #2.1 Sample catch to add catch uncertainty                                    
          TC=round(Tot.Ktch[s,]*Ktch.sex.ratio)
          
          #2.2. Vary vital rates for projection matrices 
          
          #draw max age sample
          if(scenario==1) A.sim=ASim.1[s]+1 #use same A.sim for all projections to keep same size matrix
          if(scenario==2) A.sim=ASim[s]+1
          
          #draw  age at maturity sample
          if(scenario==1) A.mat=AgeMatSim.1[s] 
          if(scenario==2) A.mat=AgeMatSim[s]
          
          
          #take matrix sample of same dimension
          
          #select matrices of A.sim dimensions
          Proj.Mat=r.numb=NULL
          if(scenario==1) 
          {
            condition=lapply(Proyec.matrix.1, function(x) nrow(x) == A.sim) 
            DATA=Proyec.matrix.1[unlist(condition)]
          }
          
          if(scenario==2) 
          {
            condition=lapply(Proyec.matrix, function(x) nrow(x) == A.sim) 
            DATA=Proyec.matrix[unlist(condition)]
          }
          
          #resample if there are <N.mdld
          if(Proy.type=="Stochas") r.numb=sample(1:length(DATA),N.mdld,replace=T)   #stochastic projection
          if(Proy.type=="Determ") r.numb=rep(1,N.mdld)                            #deterministic projection
          
          #keep only n.Yr.prol proj matrices
          Proj.Mat=DATA[r.numb]
          
          
          #2.2. Calculate selectivity
          Selectivity.sim=NULL
          if(scenario==1) Selectivity.sim=SelSim.1[s,]
          if(scenario==2) Selectivity.sim=SelSim[s,]
          
          Selectivity.sim=subset(Selectivity.sim,!is.na(Selectivity.sim)) #remove NAs
          Selectivity.sim=c(Selectivity.sim,Selectivity.sim[length(Selectivity.sim)])
          
          #2.3. Add harvesting                   
          harvest.matrix=function(matrix) 
          {
            
            H=diag(nrow(matrix))
            diag(H)=1-(U*Selectivity.sim) #apply U and selectivity
            MH=matrix%*%H
            return(MH)
          }
          Harvest.Proyec.mat=vector("list",length = n.Yr.tot)
          
          
          #pre-exploitation population size
          if(Proy.type=="Determ") Matrix=Proj.Mat[[1]] #previous
          if(Proy.type=="Stochas")
          {
            ind=unlist(lapply(Proj.Mat,lambda))
            ind1=sort(ind)
            ind2=round(length(ind1)/2)
            ind3=which(ind==ind1[ind2])[1]            #use mode lambda to calculate dens factor    
            Matrix=Proj.Mat[[ind3]]
            
            ind4=which(ind==ind1[ind2+1])[1]      
            ind5=which(ind==ind1[ind2-1])[1]      
            Mat.dummy=list(Proj.Mat[[ind3]],Proj.Mat[[ind4]],Proj.Mat[[ind5]])  
            n.dummy=sample(1:length(Mat.dummy),n.Yr.tot,replace=T)   
            Proj.Mat=Mat.dummy[n.dummy]
          }
          
          No=matrix(stable.stage(Matrix)) 
          Total.No=sum(No)
          Mature.No=sum(No[A.mat:A.sim])
          
          #2.4. Add density-depence in survival
          if(dens.dep=="YES")   
          {
            #1. calculated the density-dependent factor
            if (WhatClass=="ALL")fit.dens=optimize(fun.dens,lower=1e-5,upper=1e-1)
            if (WhatClass=="ONE")fit.dens=optimize(fun.dens,lower=1e-2,upper=100)
            Dens.factor=fit.dens$minimum
          }
          
          #Burn-in for stabilising population
          n.vec=vector("list",length =N.mdld)
          n.vec[[1]]=No
          
          for (p in 2:n.burn.in)
          {
            #Density-depence in survival
            dummy=Matrix                     
            if (WhatClass=="ALL")for (nn in 2:nrow(dummy))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
            if (WhatClass=="ONE")for (nn in 2:2)dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
            
            n.vec[[p]]=dummy%*%matrix(n.vec[[p-1]])          
          }
          
          
          
          #2.5. Population trajectories with exploitation
          
          #exploited population size
          Fvec=rep(NA,n.Yr.tot)
          #Fvec[1]=0     #No fishing in 1940
          for (p in 1:n.Yr.tot)
          {
            pp=p+n.burn.in  #index for population that considers burn in
            
            Project.M=Proj.Mat[[p]]
            #if(p%in%1:16) Project.M=Matrix  #maintain equilibrium conditions if no catch (1940:1955)
            
            
            #Density-depence in survival
            if(dens.dep=="YES")
            {
              if(TC[p]<sum(c(n.vec[[pp-1]]))*N.init)
              {
                U=TC[p]/(sum(c(n.vec[[pp-1]]))*N.init)      
                U=-log(1-U)             #annual to instaneous rate conversion
              }
              if(TC[p]>=sum(c(n.vec[[pp-1]]))*N.init) U=-log(1-0.9999999999)
              
              if(TC[p]>0)Harvest.Proyec.mat[[p]]=harvest.matrix(Project.M)  
              if(TC[p]==0)Harvest.Proyec.mat[[p]]=Project.M   
              
              dummy=Harvest.Proyec.mat[[p]]          
              if (WhatClass=="ALL")for (nn in 2:nrow(dummy))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[pp-1]]))
              if (WhatClass=="ONE")for (nn in 2:2)dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[pp-1]]))
              
              n.vec[[pp]]=dummy%*%matrix(n.vec[[pp-1]])          
              if(sum(n.vec[[pp]])<0) n.vec[[pp]]=rep(0,length(n.vec[[pp-1]]))
              Fvec[p]=U
              
            }
          }
          
          
          Pop.size=Mature.Pop.size=rep(0,(N.mdld))
          for(y in 1:(N.mdld))
          {
            Pop.size[y]=sum(n.vec[[y]])*N.init
            Mature.Pop.size[y]=sum(n.vec[[y]][A.mat:A.sim])*N.init
          }
          
          
          #2.6. Store outputs 
          Rel.abundance[s,]=Pop.size
          Mature.rel.abundance[s,]=Mature.Pop.size
          No.vec[s]=Mature.No
          F_vec[s,]=c(rep(0,n.burn.in),Fvec)
        }
        
        if(scenario==1)LH1.N1.Sel1.Us.ab[[j]]=list(Rel.abundance=Rel.abundance,
                                                   Mature.rel.abundance=Mature.rel.abundance,No.vec=No.vec,F_series=F_vec)
        if(scenario==2)LH2.N1.Sel2.Us.ab[[j]]=list(Rel.abundance=Rel.abundance,
                                                   Mature.rel.abundance=Mature.rel.abundance,No.vec=No.vec,F_series=F_vec)
      }
    }
    
    Store.Future[[h]]=list(LH1=LH1.N1.Sel1.Us.ab,LH2=LH2.N1.Sel2.Us.ab)
  })
  
  
  
  plot.mean.CI=function(DATA) 
  {
    DATA=DATA[,These.Yrs]
    MED=colMedians(DATA)
    LOW=apply(DATA, 2, function(x) quantile(x, 0.025))
    UP=apply(DATA, 2, function(x) quantile(x, 0.975))
    
    plot(1:n.Yr.tot,MED,type='l',lty=1,col=1,lwd=LWD,ylim=c(0,MED[1]*1.25),xlab="",ylab="",
         cex.axis=1,xaxt='n')
    YR=1:n.Yr.tot
    Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
    Biom.Vec <- c(LOW, tail(UP, 1), rev(UP), LOW[1])
    polygon(Year.Vec, Biom.Vec, col = "grey80", border = "grey20")
    lines(1:n.Yr.tot,MED,type='l',lty=1,col=1,lwd=2)
    axis(1,at=seq(1,n.Yr.tot,by=1),labels=F,las=1,tck=-0.02)
    axis(1,at=seq(match(1940,Yr.proj),n.Yr.tot,by=10),labels=F,tck=-0.04)
    axis(1,at=seq(match(1940,Yr.proj),n.Yr.tot,by=20),labels=F,tck=-0.06)
    
  }
  
  
  
  tiff(file="For Rick/N.traj.CI.future.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  par(mfrow=c(4,4),las=1,omi=c(.7,1,.1,.1),mai=c(.15,.3,.15,.25),mgp=c(1.5,0.8,0),xpd=F)
  for (j in 1:N.drum)
  {
    for (a in 1:N.drum)
    {
      for(i in 1:NN.future) 
      {
        plot.mean.CI(Store.Future[[j]][[a]][[i]]$Rel.abundance)
        legend("bottomleft",paste("N.0=",N.0.future[[a]][[i]]),bty='n')
        legend("bottomright",names(Store.Future[[j]])[a],bty='n',cex=1.25)
        legend("top",names(Store.Future)[j],bty='n')
        if(a==2&j==2)axis(1,at=seq(match(1940,Yr.proj),n.Yr.tot,by=20),
                          labels=seq(Yr.proj[match(1940,Yr.proj)],Yr.proj[length(Yr.proj)],by=20),
                          cex.axis=0.95,las=1,tck=-0.06)
        
      }
      
    }
    
    
  }
  mtext("Year",side=1,outer=T,line=1.5,font=1,las=0,cex=1.6)
  mtext("Total number of females",side=2,outer=T,line=1.8,font=1,las=0,cex=1.6)
  dev.off()
  
  
  get.what=function(DATA,scaler) 
  {
    MED=colMedians(DATA)
    Depletion=MED[length(MED)]/scaler
    return(list(Depletion=Depletion,N2016=MED[length(MED)]))
  }
  
  Dep.level.LH1.Future=Dep.level.LH2.Future=list(D.L.catch10=NA,D.L.catch20=NA)
  Dep.level.LH1.Future[[1]]=list(No.2000=NA, No.3000=NA, No.5000=NA,No.10000=NA)
  Dep.level.LH1.Future[[2]]=Dep.level.LH1.Future[[1]]
  N2016.LH1.Future=Dep.level.LH1.Future
  
  Dep.level.LH2.Future[[1]]=list(No.2000=NA, No.3000=NA, No.5000=NA,No.10000=NA)
  Dep.level.LH2.Future[[2]]=Dep.level.LH2.Future[[1]]
  N2016.LH2.Future=Dep.level.LH2.Future
  
  for (j in 1:N.drum)
  {
    for(i in 1:NN.future) 
    {
      
      Dep.level.LH1.Future[[j]][[i]]=get.what(Store.Future[[j]][[1]][[i]]$Rel.abundance,N.0.future[[1]][[i]])$Depletion   
      Dep.level.LH2.Future[[j]][[i]]=get.what(Store.Future[[j]][[2]][[i]]$Rel.abundance,N.0.future[[2]][[i]])$Depletion  
      
      
      N2016.LH1.Future[[j]][[i]]=get.what(Store.Future[[j]][[1]][[i]]$Rel.abundance,N.0.future[[1]][[i]])$N2016 
      N2016.LH2.Future[[j]][[i]]=get.what(Store.Future[[j]][[2]][[i]]$Rel.abundance,N.0.future[[2]][[i]])$N2016  
      
    }
    
  }
  
  #extract values
  Dep.LH1.future=unlist(Dep.level.LH1.Future)
  Dep.LH2.future=unlist(Dep.level.LH2.Future)
  
  N2016.LH1.future=unlist(N2016.LH1.Future)
  N2016.LH2.future=unlist(N2016.LH2.Future)
  
  write.table(Dep.LH1.future,"For Rick/Dep.LH1.future.csv",sep=",")
  write.table(Dep.LH2.future,"For Rick/Dep.LH2.future.csv",sep=",")
  write.table(N2016.LH1.future,"For Rick/N2016.LH1.future.csv",sep=",")
  write.table(N2016.LH2.future,"For Rick/N2016.LH2.future.csv",sep=",")
  
  
  
  
  ######
  #-- 2. 100 years projections for testing recovery ---
  ######
  Proy.type="Determ"
  Ns=10
  iterations=1:Ns
  
  yrs.sim=100   #number of years projected
  N.0.future=list(c(3000,6000,10000),c(3000,6000,10000))  #N.0 requested by Rick
  names(N.0.future)=c("LH1","LH2")
  NN.future=length(N.0.future[[1]])
  
  n.Yr.prol=yr.dummy
  
  N.drum=1
  Store.Future=vector('list',N.drum)
  
  
  N.modl.yrs=c(1838:Yr.end,(Yr.end+1):(Yr.end+yrs.sim))
  N.mdld=length(N.modl.yrs)
  
  These.Yrs=match(c(Yr.proj,2013:(Yr.end+yrs.sim)),N.modl.yrs)
  Yr.proj=c(1938:Yr.end,(Yr.end+1):(Yr.end+yrs.sim))
  
  CATCH="NO"  #no future catches
  #CATCH="YES"  #future catches set at average of last 4 years of catch data
  #CATCH="Low"  #future catch of only 10 individuals
  
  system.time(for (h in 1:N.drum)
  {
    Tot.Ktch=Store.Tot.Ktch
    if(CATCH=="YES")Tot.Ktch.future=matrix(rep(rowMeans(Tot.Ktch[,(length(Total.catch)-4):length(Total.catch)]),yrs.sim),ncol=yrs.sim)
    if(CATCH=="Low")Tot.Ktch.future=matrix(rep(10,yrs.sim*nrow(Tot.Ktch)),ncol=yrs.sim)
    if(CATCH=="NO")Tot.Ktch.future=matrix(rep(0,yrs.sim*nrow(Tot.Ktch)),ncol=yrs.sim)
    
    Tot.Ktch=cbind(Tot.Ktch,Tot.Ktch.future)
    
    n.Yr.tot=ncol(Tot.Ktch)
    
    
    #rel abundance
    LH1.N1.Sel1.Us.ab=LH2.N1.Sel2.Us.ab=vector('list',length=NN.future)
    names(LH1.N1.Sel1.Us.ab)=N.0.future[[1]]
    names(LH2.N1.Sel2.Us.ab)=N.0.future[[2]]
    
    #Risk analysis
    for (x in 1:length(Life.hist.scenarios))
    {
      scenario=Life.hist.scenarios[[x]]
      for(j in 1:NN.future) 
      {
        if(scenario==1)N.init=N.0.future[[1]][j]
        if(scenario==2)N.init=N.0.future[[2]][j]
        
        #1. Create elements to fill in
        Pop.size.ratio=No.vec=rep(NA,length = Ns)
        Rel.abundance=Mature.rel.abundance=F_vec=matrix(NA,ncol=N.mdld,nrow = Ns)
        
        for (s in iterations)
        {
          #2.1 Sample catch to add catch uncertainty                                    
          TC=round(Tot.Ktch[s,]*Ktch.sex.ratio)
          
          #2.2. Vary vital rates for projection matrices 
          
          #draw max age sample
          if(scenario==1) A.sim=ASim.1[s]+1 #use same A.sim for all projections to keep same size matrix
          if(scenario==2) A.sim=ASim[s]+1
          
          #draw  age at maturity sample
          if(scenario==1) A.mat=AgeMatSim.1[s] 
          if(scenario==2) A.mat=AgeMatSim[s]
          
          
          #take matrix sample of same dimension
          
          #select matrices of A.sim dimensions
          Proj.Mat=r.numb=NULL
          if(scenario==1) 
          {
            condition=lapply(Proyec.matrix.1, function(x) nrow(x) == A.sim) 
            DATA=Proyec.matrix.1[unlist(condition)]
          }
          
          if(scenario==2) 
          {
            condition=lapply(Proyec.matrix, function(x) nrow(x) == A.sim) 
            DATA=Proyec.matrix[unlist(condition)]
          }
          
          #resample if there are <N.mdld
          if(Proy.type=="Stochas") r.numb=sample(1:length(DATA),N.mdld,replace=T)   #stochastic projection
          if(Proy.type=="Determ") r.numb=rep(1,N.mdld)                            #deterministic projection
          
          #keep only n.Yr.prol proj matrices
          Proj.Mat=DATA[r.numb]
          
          
          #2.2. Calculate selectivity
          Selectivity.sim=NULL
          if(scenario==1) Selectivity.sim=SelSim.1[s,]
          if(scenario==2) Selectivity.sim=SelSim[s,]
          
          Selectivity.sim=subset(Selectivity.sim,!is.na(Selectivity.sim)) #remove NAs
          Selectivity.sim=c(Selectivity.sim,Selectivity.sim[length(Selectivity.sim)])
          
          #2.3. Add harvesting                   
          harvest.matrix=function(matrix) 
          {
            
            H=diag(nrow(matrix))
            diag(H)=1-(U*Selectivity.sim) #apply U and selectivity
            MH=matrix%*%H
            return(MH)
          }
          Harvest.Proyec.mat=vector("list",length = n.Yr.tot)
          
          
          #pre-exploitation population size
          if(Proy.type=="Determ") Matrix=Proj.Mat[[1]] #previous
          if(Proy.type=="Stochas")
          {
            ind=unlist(lapply(Proj.Mat,lambda))
            ind1=sort(ind)
            ind2=round(length(ind1)/2)
            ind3=which(ind==ind1[ind2])[1]            #use mode lambda to calculate dens factor    
            Matrix=Proj.Mat[[ind3]]
            
            ind4=which(ind==ind1[ind2+1])[1]      
            ind5=which(ind==ind1[ind2-1])[1]      
            Mat.dummy=list(Proj.Mat[[ind3]],Proj.Mat[[ind4]],Proj.Mat[[ind5]])  
            n.dummy=sample(1:length(Mat.dummy),n.Yr.tot,replace=T)   
            Proj.Mat=Mat.dummy[n.dummy]
          }
          
          No=matrix(stable.stage(Matrix)) 
          Total.No=sum(No)
          Mature.No=sum(No[A.mat:A.sim])
          
          #2.4. Add density-depence in survival
          if(dens.dep=="YES")   
          {
            #1. calculated the density-dependent factor
            if (WhatClass=="ALL")fit.dens=optimize(fun.dens,lower=1e-5,upper=1e-1)
            if (WhatClass=="ONE")fit.dens=optimize(fun.dens,lower=1e-2,upper=100)
            Dens.factor=fit.dens$minimum
          }
          
          #Burn-in for stabilising population
          n.vec=vector("list",length =N.mdld)
          n.vec[[1]]=No
          
          for (p in 2:n.burn.in)
          {
            #Density-depence in survival
            dummy=Matrix                     
            if (WhatClass=="ALL")for (nn in 2:nrow(dummy))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
            if (WhatClass=="ONE")for (nn in 2:2)dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
            
            n.vec[[p]]=dummy%*%matrix(n.vec[[p-1]])          
          }
          
          
          
          #2.5. Population trajectories with exploitation
          
          #exploited population size
          Fvec=rep(NA,n.Yr.tot)
          #Fvec[1]=0     #No fishing in 1940
          for (p in 1:n.Yr.tot)
          {
            pp=p+n.burn.in  #index for population that considers burn in
            
            Project.M=Proj.Mat[[p]]
            #if(p%in%1:16) Project.M=Matrix  #maintain equilibrium conditions if no catch (1940:1955)
            
            
            #Density-depence in survival
            if(dens.dep=="YES")
            {
              if(TC[p]<sum(c(n.vec[[pp-1]]))*N.init)
              {
                U=TC[p]/(sum(c(n.vec[[pp-1]]))*N.init)      
                U=-log(1-U)             #annual to instaneous rate conversion
              }
              if(TC[p]>=sum(c(n.vec[[pp-1]]))*N.init) U=-log(1-0.9999999999)
              
              if(TC[p]>0)Harvest.Proyec.mat[[p]]=harvest.matrix(Project.M)  
              if(TC[p]==0)Harvest.Proyec.mat[[p]]=Project.M   
              
              dummy=Harvest.Proyec.mat[[p]]          
              if (WhatClass=="ALL")for (nn in 2:nrow(dummy))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[pp-1]]))
              if (WhatClass=="ONE")for (nn in 2:2)dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[pp-1]]))
              
              n.vec[[pp]]=dummy%*%matrix(n.vec[[pp-1]])          
              if(sum(n.vec[[pp]])<0) n.vec[[pp]]=rep(0,length(n.vec[[pp-1]]))
              Fvec[p]=U
              
            }
          }
          
          
          Pop.size=Mature.Pop.size=rep(0,(N.mdld))
          for(y in 1:(N.mdld))
          {
            Pop.size[y]=sum(n.vec[[y]])*N.init
            Mature.Pop.size[y]=sum(n.vec[[y]][A.mat:A.sim])*N.init
          }
          
          
          #2.6. Store outputs 
          Rel.abundance[s,]=Pop.size
          Mature.rel.abundance[s,]=Mature.Pop.size
          No.vec[s]=Mature.No
          F_vec[s,]=c(rep(0,n.burn.in),Fvec)
        }
        
        if(scenario==1)LH1.N1.Sel1.Us.ab[[j]]=list(Rel.abundance=Rel.abundance,
                                                   Mature.rel.abundance=Mature.rel.abundance,No.vec=No.vec,F_series=F_vec)
        if(scenario==2)LH2.N1.Sel2.Us.ab[[j]]=list(Rel.abundance=Rel.abundance,
                                                   Mature.rel.abundance=Mature.rel.abundance,No.vec=No.vec,F_series=F_vec)
      }
    }
    
    Store.Future[[h]]=list(LH1=LH1.N1.Sel1.Us.ab,LH2=LH2.N1.Sel2.Us.ab)
    
  })
  
  
  
  plot.mean.CI=function(DATA) 
  {
    n.Yr.tot=length(These.Yrs)
    DATA=DATA[,These.Yrs]
    MED=colMedians(DATA)
    LOW=apply(DATA, 2, function(x) quantile(x, 0.025))
    UP=apply(DATA, 2, function(x) quantile(x, 0.975))
    
    plot(1:length(These.Yrs),MED,type='l',lty=1,col=1,lwd=LWD,ylim=c(0,MED[1]*1.25),xlab="",ylab="",
         cex.axis=1,xaxt='n')
    YR=1:length(These.Yrs)
    Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
    Biom.Vec <- c(LOW, tail(UP, 1), rev(UP), LOW[1])
    polygon(Year.Vec, Biom.Vec, col = "grey80", border = "grey20")
    lines(1:n.Yr.tot,MED,type='l',lty=1,col=1,lwd=2)
    axis(1,at=seq(1,n.Yr.tot,by=1),labels=F,las=1,tck=-0.02)
    axis(1,at=seq(match(1940,Yr.proj),n.Yr.tot,by=10),labels=F,tck=-0.04)
    axis(1,at=seq(match(1940,Yr.proj),n.Yr.tot,by=20),labels=F,tck=-0.06)
    
  }
  
  
  
  tiff(file="Future.Proj/NoCatch.future.determ1000.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
  par(mfrow=c(2,3),las=1,omi=c(.7,1,.1,.1),mai=c(.15,.3,.15,.25),mgp=c(1.5,0.8,0),xpd=F)
  for (j in 1:N.drum)
  {
    for (a in 1:2)
    {
      for(i in 1:NN.future) 
      {
        plot.mean.CI(Store.Future[[j]][[a]][[i]]$Rel.abundance)
        legend("bottomleft",paste("N.0=",N.0.future[[a]][[i]]),bty='n')
        legend("bottomright",names(Store.Future[[j]])[a],bty='n',cex=1.25)
        if(a==2&j==1)axis(1,at=seq(match(1940,Yr.proj),n.Yr.tot,by=20),
                          labels=seq(Yr.proj[match(1940,Yr.proj)],Yr.proj[length(Yr.proj)],by=20),
                          cex.axis=0.95,las=1,tck=-0.06)
        
      }
      
    }
    
    
  }
  mtext("Year",side=1,outer=T,line=1.5,font=1,las=0,cex=1.6)
  mtext("Total number of females",side=2,outer=T,line=1.8,font=1,las=0,cex=1.6)
  dev.off()
  
  
  
  
  
  
}

#######################################################################################