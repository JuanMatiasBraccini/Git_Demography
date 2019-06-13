#Create variance covariance matrix using Wintner and Cliff 1999 data

#NOTE: data are not good, cannot get model to converge!!!!

graphics.off ()
rm(list=ls(all=TRUE))

#LIBRARIES
#no libraries used in this script

#DATA SECTION
setwd("H:/Matias WA Fisheries/Analyses/Demography/White shark")

Data=read.csv("Henry.csv")


#PARAMETER SECTION
Linf=544; k=0.07; to=-3; cv=0.5;  	##CV is coef of var for length at age; lh is 50% selectivity, sh is the SD of lh
theta=c(Linf=Linf, k=k, to=to, cv=cv)		#vector of parameters of interest


#PROCEDURE SECTION
age=Data$AGE
#l=Data$TL  #total length
l=Data$PCL_cm  # PCL

plot(age,l)

vonb=function(theta)
{  					#a and l are the observed age and length data, assumed to be global in scope
	with(as.list(theta),		#theta defined as a list to be able to use the objects without naming them. A list is a bag (e.g. theta) and you can get stuff out of the bag. Can have a list inside a list
	{
	
    l.hat<<-Linf*(1-exp(-k*(age-to)))	#the von berta function
		epsilon=l-l.hat				#goodness of fit
		nloglike=-1.0 *sum(dnorm(epsilon,0,cv*l.hat,log=T))		#the likelihood function
		return(list(nloglike=nloglike,l.hat=l.hat))
	})

}

#3.b Estimating growth parameters
theta1=c(Linf=Linf, k=k, to=to)  	#new vector of only parameters of interest
cv=0.08
"objfun"=function(theta1){vonb(theta1)$nloglike}	
fit=optim(theta1,objfun, method="BFGS",hessian=T)

V=solve(fit$hessian)

