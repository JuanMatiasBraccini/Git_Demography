
#script for adding density dependence to matrix projections


library(popbio)     #for solving matrices

#1. Control scenarios
#Add density dependence or not
dens.dep="YES"
#dens.dep="NO"

#where is density dependence added
WhatClass="ALL"
#WhatClass="ONE"


#2. Projection matrix
Matrix=matrix(c(0,0,2,2.1,.5,0,0,0,0,.8,0,0,0,0,.9,0),ncol=4,byrow=T)

#3. Stable age distributions
No=matrix(stable.stage(Matrix))
Total.No=sum(No)

#4. Calculate density-dependence factor
fun.dens=function(dens)
{
  #projection matrix
  MH=Matrix #add density dependence to unfished population
  
  #add survival density-dependence to projection matrix
  if (WhatClass=="ALL")for (nn in 2:nrow(MH))MH[nn,(nn-1)]=MH[nn,(nn-1)]*exp(-dens*Total.No)
  if (WhatClass=="ONE")for (nn in 2:2)MH[nn,(nn-1)]=MH[nn,(nn-1)]*exp(-dens*Total.No)
  
  #calculate lambda
  LAMBDA=lambda(MH)
  
  #obj function
  epsilon=(1-LAMBDA)^2
  
  return(list(epsilon=epsilon,LAMBDA=LAMBDA))
}

fn.F=function(dens)fun.dens(dens)$epsilon

fit.dens=optimize(fn.F,lower=0.0001,upper=10)
Dens.factor=fit.dens$minimum

#profile likelihood
a=seq(0.0001,10,.01)
b=sapply(a,fn.F)
plot(a,b,col=2,ylab="negloglike",xlab="dens-dep factor")
text(Dens.factor,min(b),round(Dens.factor,3),cex=2)

#lambda
a=seq(0,10,.1)
store.lamb=rep(NA,length(a))
for (i in 1:length(a))store.lamb[i]=fun.dens(a[i])$LAMBDA
plot(a,store.lamb,col=2,ylab="lambda",xlab="dens-dep factor")
text(Dens.factor,fun.dens(Dens.factor)$LAMBDA,round(Dens.factor,3),cex=2)


#5. Project population with or without density dependence

n.Yr.prol=50  #number of years projected


#exploited population size
n.vec=vector("list",length = n.Yr.prol)
prop.deplet=0.2
n.vec[[1]]=(No)*prop.deplet #initial population size
for (p in 2:n.Yr.prol)
{
  #no density dependence
  if(dens.dep=="NO") n.vec[[p]]=Matrix%*%matrix(n.vec[[p-1]])

  #Density-depence in year-class 1 survival
  if(dens.dep=="YES")
  {
    dummy=Matrix
    
    
    if (WhatClass=="ALL")for (nn in 2:nrow(dummy))dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
    if (WhatClass=="ONE")for (nn in 2:2)dummy[nn,(nn-1)]=dummy[nn,(nn-1)]*exp(-Dens.factor*sum(n.vec[[p-1]]))
    
    
    n.vec[[p]]=dummy%*%matrix(n.vec[[p-1]])
  }
  
}

Pop.size=rep(0,n.Yr.prol)
for(y in 1:n.Yr.prol) Pop.size[y]=sum(n.vec[[y]])
plot(Pop.size,col=2)