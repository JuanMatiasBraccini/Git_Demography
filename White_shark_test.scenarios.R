Ns=1000
iterations=1:Ns


Prob.Pop.doubling=rep(NA,length = 5)
Pop.project=vector("list",length = 5)
for(aa in 1:5)
{

#2.1. Set selectivity scenarios

Selectivity.SIM=Selectivity.SIM.1=vector("list",length = Ns)

scenario.sel=1
for (s in iterations) Selectivity.SIM.1[[s]]=Sel.fn(ASim[s],LinfSim[s],kSim[s],toSim[s])

scenario.sel=2
for (s in iterations) Selectivity.SIM[[s]]=Sel.fn(ASim[s],LinfSim[s],kSim[s],toSim[s])

SelSim=add.missing.age(Selectivity.SIM)
SelSim.1=add.missing.age(Selectivity.SIM.1)

#1. Set scenarios
#biological scenario
scenario=Life.hist.scenarios[[2]]

#Selectivity scenario
scenario.sel=Sel.scenarios[[2]]

#Harvest rate scenario
scenario.U=U.scenarios[[aa]]

#N1998 scenario
scenario.N1998=N1998.scenarios[[1]]


#2. Create elements to fill in
Pop.size.ratio=rep(NA,length = Ns)
store.pop.proy=NULL
#3. Monte Carlo loop
 for (s in iterations)
 {

  #2. Vary vital rates for projection matrices 
  
  #draw max age sample
  A.sim=ASim[s] #use same A.sim for all projections to keep same size matrix
  
  #take matrix sample of same dimension
  Proj.Mat=r.numb=NULL

# 
#   if(scenario==1) 
#     {
#     condition=lapply(Proyec.matrix.1, function(x) nrow(x) == A.sim)
#     DATA=Proyec.matrix.1[unlist(condition)]
#     Proj.Mat=DATA[[1]]
#     }
#   if(scenario==2) 
#   {
#     condition=lapply(Proyec.matrix, function(x) nrow(x) == A.sim)
#     DATA=Proyec.matrix[unlist(condition)]
#     Proj.Mat=DATA[[1]]
#   }
  if(scenario==1) 
  {
    #select matrices of A.sim dimensions
    condition=lapply(Proyec.matrix.1, function(x) nrow(x) == A.sim) 
    DATA=Proyec.matrix.1[unlist(condition)]
    if(A.sim==A)r.numb=sample(1:length(DATA),n.Yr.prol,replace=T) #resample for A.sim 60 (there are <15) 
    if(A.sim<A)r.numb=sample(1:length(DATA),n.Yr.prol,replace=F)
    
    #keep only 15
    Proj.Mat=DATA[r.numb]
  }
  
  if(scenario==2) 
  {
    #select matrices of A.sim dimensions
    condition=lapply(Proyec.matrix, function(x) nrow(x) == A.sim) 
    DATA=Proyec.matrix[unlist(condition)]
    if(A.sim==A)r.numb=sample(1:length(DATA),n.Yr.prol,replace=T)  
    if(A.sim<A)r.numb=sample(1:length(DATA),n.Yr.prol,replace=F)
    
    #keep only 15
    Proj.Mat=DATA[r.numb]
  }
  
  
  #3. Calculate selectivity
  Selectivity.sim=NULL
  if(scenario.sel==1) Selectivity.sim=SelSim.1[s,]
  if(scenario.sel==2) Selectivity.sim=SelSim[s,]
  
  Selectivity.sim=subset(Selectivity.sim,!is.na(Selectivity.sim)) #remove NAs
  
  #4. Add harvesting
  harvest.matrix=function(matrix,U) 
  {
    H=diag(nrow(matrix))
    diag(H)=1-(U*Selectivity.sim) #apply U and selectivity
    MH=matrix%*%H
    return(MH)
  }
  Harvest.Proyec.mat=vector("list",length = n.Yr.prol)
  for (h in 1:n.Yr.prol) Harvest.Proyec.mat[[h]]=harvest.matrix(Proj.Mat[[h]],scenario.U[h])
#  Harvest.Proyec.mat=harvest.matrix(Proj.Mat,scenario.U[[1]][1])  
  
  #5. Project population into future
#   nn=matrix(stable.stage(Harvest.Proyec.mat)*scenario.N1998[[1]])
#   p<-pop.projection(Harvest.Proyec.mat,nn, 15)   #project population
# #  plot(p$pop.sizes) #plot pop size
#  lines(p$pop.sizes,col=4) 
  n.vec=vector("list",length = n.Yr.prol)
  n.vec[[1]]=matrix(stable.stage(Harvest.Proyec.mat[[1]])*scenario.N1998)
  for (p in 2:n.Yr.prol)
  {
    n.vec[[p]]=Harvest.Proyec.mat[[p]]%*%matrix(n.vec[[p-1]])
  }
  Pop.size=rep(0,n.Yr.prol)
  for(y in 1:n.Yr.prol) Pop.size[y]=sum(n.vec[[y]])

#   if(aa==1)plot(Pop.size,col=aa,ylim=c(0.7,max(Pop.size)))
# if(aa>1)points(Pop.size,col=aa)
  #6. Calculate population size ratio
  Pop.size.ratio[s]=Pop.size[length(Pop.size)]/Pop.size[1]
  
  store.pop.proy=rbind(store.pop.proy,Pop.size)
 }
Pop.project[[aa]]=store.pop.proy
# 
# #Calculate reference points
 Prop.Pop.double=subset(Pop.size.ratio,Pop.size.ratio>=Biom.ref.point)
 Pop.size.ratio=subset(Pop.size.ratio,!is.na(Pop.size.ratio))
 Prob.Pop.doubling[aa]=length(Prop.Pop.double)/length(Pop.size.ratio)
}

plot(Prob.Pop.doubling)

# par(mfcol=c(3,2),omi=c(.6,.9,.4,.1),mai=c(.15,.15,.15,.15))
# for (i in 1:5){
#   plot(Pop.project[[i]][1,],type='l',ylim=c(0,max(Pop.project[[i]])))
#   for(j in 2:10) lines(Pop.project[[i]][j,],type='l')
#   legend('topleft',paste("u", i))




test=NULL
for(aaa in 1:5){
  test=rbind(test,Risk.fn(Life.hist.scenarios[[2]],N1998.scenarios[[1]],
                          U.scenarios[[aaa]],Sel.scenarios[[2]]))}