# ------ Find No given CNe  --------
#note: run Demography script up to #Run Risk Analysis for each scenario"

#Assumptions:
# Blower et al 2011
CNe.mean=693


# Must assume a female:male ratio of CNe (i.e. do all female and all males reproduce? if so, then 1:1)
CNe.sex.ratio=list(0.5,0.25)  #female:male sex ratios of 1:1 and 1:3


# Must assume a ratio of CNe:CNc (i.e. how much bigger is census population?)
#CNe.CNc.ratio=2 #i.e. CNc= 2 x CNe
CNe.CNc.ratio=list(1,5)

#Tot.KTCH=TOTAL.CATCH[[1]]    #define a catch series
Tot.KTCH=list(TOTAL.CATCH[[3]],TOTAL.CATCH[[2]]) #lowest and largest catch)


LH.scenario=1:2   #life history scenarios



#PROCEDURE

#Estimate No under different assumptions of genetics stuff, catch and life history

fun.findNo=function(NO)
{

    #initial population size 
    N.init=NO
    
    #2.1 Sample catch to add catch uncertainty                                    
    TC=round(Tot.Ktch*Ktch.sex.ratio)
    
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
    r.numb=rep(1,N.mdld)          #deterministic projection (doesn't work for stochastic, cannot estimate)
    
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
    Matrix=Proj.Mat[[1]] #previous
  
    No=matrix(stable.stage(Matrix)) 
    Total.No=sum(No)
    Mature.No=sum(No[A.mat:A.sim])
    
    #2.4. Add density-depence in survival
    if(dens.dep=="YES")   
    {
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
      Mature.Pop.size[y]=sum(n.vec[[y]][(A.mat+1):A.sim])*N.init
    }

  
    #Objective function
    epsilon=(Ncurrent-Mature.Pop.size[length(Mature.Pop.size)])^2 

  
  return(list(epsilon=epsilon,N.mat=Mature.Pop.size[length(Mature.Pop.size)]))
}


fn.minimize=function(NO) fun.findNo(NO=NO)$epsilon


fn.plot.like=function(No.vec,Dat)
{
  plot(No.vec,Dat$Likelihood,type='l',xlab="",cex.lab=1.25,ylab="")
  points(Dat$NO.pred,min(Dat$Likelihood),cex=2,col=2,pch=19)
  text(Dat$NO.pred,mean(Dat$Likelihood),round(Dat$NO.pred),col=2,cex=1.25)
}
  
 
fn.plot.No=function(No.vec,N.MAT)
{
  plot(No.vec,N.MAT,type='l',main="Profile for No",xlab="Assumed No",ylab="Number of mature females",cex.lab=1.25)
  abline(h=Ncurrent,lwd=2,col=2)
  text(No.vec[5],Ncurrent*1.25,"CNe",col=2)
}


#Run over all scenarios

No.vec=seq(0,20000,1000)

NOstart=10000

store4=vector('list',length(LH.scenario))
names(store4)=paste("LH",LH.scenario)
system.time(for(i in 1:length(LH.scenario))
{
  scenario=LH.scenario[i]
  store3=vector('list',length(Tot.KTCH))
  names(store3)=c("C3","C2")
  for (j in 1:length(Tot.KTCH))
  {
    TKTCH=Tot.KTCH[[j]]
    Tot.Ktch=colMeans(TKTCH)

    s=round(runif(1,1,100))        
    
    store2=vector('list',length(CNe.sex.ratio))
    names(store2)=paste("CNe.sex.ratio",CNe.sex.ratio)
    for(p in 1:length(CNe.sex.ratio))
    {
      CNEsex=CNe.sex.ratio[[p]]
      
      store1=vector('list',length(CNe.CNc.ratio))
      names(store1)=paste("CNe.CNc.ratio",CNe.CNc.ratio)
      for(x in 1:length(CNe.CNc.ratio))
      {
        #Current census number of breeding (mature) females
        CNe.sample=CNe.mean
        Ncurrent=round(CNe.sample*CNEsex*CNe.CNc.ratio[[x]])
        
        N.vec=No.vec
        if(x%in%2:3)N.vec=No.vec*3
        
        #estimate No
        N.fit <- nlminb( start=NOstart, objective=fn.minimize,
                         control=list(eval.max=500, iter.max=500  )   )
        NO.pred=N.fit$par
         
        #Profile likelihood
        PROFILE=sapply(N.vec,fun.findNo)
        Likelihood=unlist(PROFILE[1,])
        N.MAT=unlist(PROFILE[2,])
        
        
        store1[[x]]=list(NO.pred=NO.pred,Likelihood=Likelihood,N.MAT=N.MAT)
      }
      store2[[p]]=store1
    }
    store3[[j]]=store2
  }
  store4[[i]]=store3
})


KTCh.lab=c("C3","C2")

#tiff(file="Find.No.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(4,4),mai=c(.3,.3,.01,.01),oma=c(3,3,.1,.1),mgp=c(2, 1, 0)) 
for(i in 1:length(LH.scenario))
{
  for (j in 1:length(Tot.KTCH))
  {
    for(p in 1:length(CNe.sex.ratio))
    {
      for(x in 1:length(CNe.CNc.ratio))
      {
        N.vec=No.vec
        if(x%in%2:3)N.vec=No.vec*3  
        fn.plot.like(N.vec,store4[[i]][[j]][[p]][[x]])
          legend("topleft",c(paste("LH",LH.scenario[i]),paste("TC",KTCh.lab[j]),paste("Sex",CNe.sex.ratio[[p]]),
              paste("Nc.Ne",CNe.CNc.ratio[[x]])),bty='n',col="grey60",cex=0.9)
      }
    }
  }
}
mtext("Likelihood",2,outer=T,line=1.25,las=3,cex=1.5)
mtext("Female No",1,outer=T,line=1,cex=1.5)
#dev.off()
