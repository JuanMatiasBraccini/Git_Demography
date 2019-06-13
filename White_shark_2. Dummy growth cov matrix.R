#Create variance covariance matrix simulating age and growth study of white shark for
#obtaining similar order or variance in parameter estimates

graphics.off ()
rm(list=ls(all=TRUE))

#LIBRARIES
#no libraries used in this script

#DATA SECTION
A=15				#maximum age
# X=7				#maximum size
# N=500				#sample size (sampled in lake)
# bw= 0.2				#bin width for selectivity, which is equal to fish measurement error. Change to see how it affects lengh-age (and hence growth parameters)
# x=seq(1,X, by=bw)		#size range for selectivity curve
# nperbin=100			#max number of fish per size bin (for subsampling the sample for ageing)
# 
# 1. Simulation phase
# the prob of sampling a fish of length (x) is a function of the selectivty, the survivorship and the age-length key

#PARAMETER SECTION
Linf=8; k=0.07; to=-3; cv=0.2; m=0.3; lh=45; sh=.5; gama=0.35		##CV is coef of var for length at age; lh is 50% selectivity, sh is the SD of lh
theta=c(Linf=Linf, k=k, to=to, cv=cv, m=m, lh=lh, sh=sh)		#vector of parameters of interest

#PROCEDURE SECTION


age=1:A
	# 1.a) VBGF
	La=Linf*(1-exp(-k*(age-to)))				

#l=jitter(rep(La,10),100)
l=rep(La,10)
l=ifelse(l>5.5,l*runif(l,exp(-0.4),exp(0.2)),ifelse(l<3,l*runif(l,exp(-0.2),exp(0.2)),
        l*runif(l,exp(-0.09),exp(0.09))))
l=ifelse(l>6.5,l*runif(l,exp(-0.2),exp(0.3)),l)
a=rep(age,10)
dummy=data.frame(a,l)
dummy=dummy[order(dummy$l),]
dummy=subset(dummy,a<10)
dummy[(nrow(dummy)+1):(nrow(dummy)+3),]=data.frame(c(11,12,14),c(5.2,5.5,5.2))
l=dummy$l
a=dummy$a
# 3. Fitting model to data
# 3.a Drawing true growth curve using given parameters
vonb=function(theta)
{						#a and l are the observed age and length data, assumed to be global in scope
	with(as.list(theta),		#theta defined as a list to be able to use the objects without naming them. A list is a bag (e.g. theta) and you can get stuff out of the bag. Can have a list inside a list
	{
		l.hat<<-Linf*(1-exp(-k*(a-to)))	#the von berta function
		epsilon=l-l.hat				#goodness of fit
		nloglike=-1.0 *sum(dnorm(epsilon,0,cv*l.hat,log=T))		#the likelihood function
		return(list(nloglike=nloglike,l.hat=l.hat))
	})

}

#3.b Estimating growth parameters
#theta1=c(Linf=6, k=0.06, to=-3,cv=0.3)		#new vector of only parameters of interest
theta1=c(Linf=7, k=0.06, to=-3)  	#new vector of only parameters of interest

"objfun"=function(theta1){vonb(theta1)$nloglike}	
fit=optim(theta1,objfun, method="BFGS", hessian=T) 
fit



#REPORT SECTION

# 3. Fitting model to data
plot(a,l, xlab="Age (years)",ylab="Length (mm)", las=1)			# plot of sampled age length data
#plot(jitter(a),jitter(l))							#this shows the number of fish at age and length
#points(a,vonb(theta)$l.hat,col=2, pch=19,cex=3)		#predicted age length relationship
											#this shows how the sample data is biased due to selectivity when compared to the true growth curve, in red (the one we simulated)
											#unbiased sample will be obtained if selectivity is 1 for all age classes by setting lh=-450 for example
lines(a[order(a)],vonb(theta)$l.hat[order(vonb(theta)$l.hat)],col=3)			#True. need to order the randomly generated ages and length
lines(a[order(a)],vonb(fit$par)$l.hat[order(vonb(fit$par)$l.hat)],col=2)		# Predicted. also need to order
legend("bottomright",c("Observation", "True", "Estimated"),lty=c(2,1,1),lwd=c(1,2,1),col=c(1,3,2),bty="n")



V=solve(fit$hessian)  #getting the variance covariance matrix

meanValues=fit$par
sd.Linf=(V[1,1])^0.5
sd.k=(V[2,2])^0.5
sd.to=(V[3,3])^0.5
