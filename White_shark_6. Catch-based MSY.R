#  SCRIPT FOR CALCULATING MSY USING CATCH-BASED METHOD

# TO DO:

rm(list=ls(all=TRUE))

# ---- DATA SECTION -----

yr=1975:2012

#Catch in numbers
#Sourced from "3.WhiteShark_demography.Leslie.R"
median.Total.catch=c(66,70,83,97,91,103,124,143,167,217,229,245,314,275,254,237,224,211,186,
              183,178,168,188,185,171,143,132,134,146,146,130,112,104,104,111,116,107,96)
upper95.Total.catch=c(86,90,109,127,123,139,168,193,226,293,310,331,425,372,343,320,303,285,
              252,247,240,228,256,251,232,195,180,182,198,198,177,152,142,141,150,158,146,130)
lower95.Total.catch=c(50,52,62,71,64,72,87,100,116,151,159,170,218,191,176,164,156,146,129,127,
              124,116,130,128,118,99,91,92,101,101,90,77,72,72,76,80,74,66)


#Reported size of caught whites, TL in metres (from Steve T)
Rep.TL=c(5.2,4.88,4.88,4.6,4.5,4.42,4.3,4.27,4.27,4,3.96,3.7,3.66,3.66,3.5,3.5,3.4,3.35,3.05,3.05,3.05,3,3,3,3,3,
         3,3,3,3,3,3,2.74,2.74,2.7,2.5,2.5,2.5,2.4,2.4,2.29,2.13,2.1,2,2,2,2,1.83,1.8,1.8,1.74,1.6,1.55,1.52,1.5,1.47,
         1.41,1.37,1.3,1.25,1.22,1.2,1.2,1.2,0.92)


# ---- PARAMETERS SECTION -----

  #r priors (productivity paper)
r.scen1=0.072
r.scen1.sd=0.028

r.scen2=0.022
r.scen2.sd=0.012

#TL-TWT (Mollet & Caillet 1996)
b=7.856
a=3.1087


#Range of K.upper tested
K.UPPER= c(50, 100, 200) 



#Start and Final biomass ranges
startbio.list   <- list(c(0.7,.9))
finalbio.list   <- list(c(0.2, 0.7))

#NUMBER OF ITERATIONS (may no need as high for gummy and whiskery)
n           <- 100000  

#PROCESS ERROR; 0 if deterministic model; 0.05 reasonable value? 0.2 is too high
sigR        <- 0.05   #source Eval Lai      


CV=0.2  #SD used for r lognormal prior

yr.n=length(yr)


#Define if using specific priors(Yes) or default priors (No)
user="Yes"
#user="No"



# ---- ANALYSIS SECTION ----

#1. Calculate catch in tons from reconstructed catch in numbers
  #1.1. histogram of reported sizes
HIST.TL=hist(Rep.TL)
Prop.TL=HIST.TL$counts/sum(HIST.TL$counts)
names(Prop.TL)=HIST.TL$mids

  #1.2. size class mean weight
Rep.TWT=b*HIST.TL$mids^a 

  #1.3. convert catch in numbers to weights using size (weight) distribution
fn.katch=function(dat)
{
  Catch.Mat=matrix(rep(dat,length(Prop.TL)),nrow=length(Prop.TL),byrow=T)
  
  Catch.Mat=Prop.TL*Catch.Mat
  Catch.Mat.kg=Catch.Mat*Rep.TWT
  return(colSums(Catch.Mat.kg)/1000)
}

median.Total.catch.tons=fn.katch(median.Total.catch)
lower.Total.catch.tons=fn.katch(lower95.Total.catch)
upper.Total.catch.tons=fn.katch(upper95.Total.catch)


  # 1.4. Resampling catches within confidence intervals
# SD.Total.catch.tons=0.17*median.Total.catch.tons
# 
# #test that SD match CI
# SamplE=function(Mean,SD) rnorm(1000,Mean,SD)
# Test=mapply(SamplE,median.Total.catch.tons,SD.Total.catch.tons)
# 
# UP=apply(Test, 2, function(x)quantile(x,probs=0.975))
# LOW=apply(Test, 2, function(x)quantile(x,probs=0.025))
# MED=apply(Test, 2, function(x)mean(x))
# 
# plot(MED,ylim=c(0,150))
# lines(UP)
# lines(LOW)
# points(median.Total.catch.tons,col=2,pch=19)
# lines(lower.Total.catch.tons,col=2)
# lines(upper.Total.catch.tons,col=2)
# 


#2. Calculate MSY


#---PROCEDURE SECTION-----

#FUNCTIONS
#Sample the catch
SamplE=function(Mean,SD) rnorm(1,Mean,SD)

#2.1. Default
#note: this assumes C trend is Abundance trend proxy, it ignores targeting changes, etc

#2.1.1 biomass range at start of time series
B.i.def.fn=function(DAT)if(DAT$ct[nrow(DAT)]/max(DAT$ct)<0.5)B.i=c(.5,.9)else{B.i=c(.3,.6)}

#2.1.2 biomass range at after last catch 
B.f.def.fn=function(DAT)if(DAT$ct[nrow(DAT)]/max(DAT$ct)>0.5)B.f=c(.3,.7)else{B.f=c(.01,.4)}


#Surplus production model for calculating biomass
.schaefer  <- function(theta)
{
  with(as.list(theta), {  ## for all combinations of ri & ki
    bt=vector()
    ell = 0  ## initialize ell
    for (j in startbt)
    {
      if(ell == 0) 
      {
        bt[1]=j*k*exp(rnorm(1,0, sigR))  ## set biomass in first year
        for(i in 1:nyr) ## for all years in the time series
        {
          xt=rnorm(1,0, sigR)
          bt[i+1]=(bt[i]+r*bt[i]*(1-bt[i]/k)-ct[i])*exp(xt) ## calculate biomass as function of previous year's biomass plus net production minus catch
        }
        
        #Bernoulli likelihood, assign 0 or 1 to each combination of r and k
        ell = 0
        if(bt[nyr+1]/k>=lam1 && bt[nyr+1]/k <=lam2 && min(bt) > 0 && max(bt) <=k && bt[which(yr==interyr)]/k>=interbio[1] && bt[which(yr==interyr)]/k<=interbio[2]) 
          ell = 1
      }  
    }
    return(list(ell=ell))
    
    
  })
}

#This function conducts the stock reduction analysis for N trials
#args:
#  theta - a list object containing:
#    r (lower and upper bounds for r)
#  	k (lower and upper bounds for k)
#		lambda (limits for current depletion)
sraMSY	<-function(theta, N)
{
  with(as.list(theta), 
{
  ## get N values of r, assign to ri
  #lognormal r  
  if(user=="Yes")
  {
    logSD=sqrt(log(1+((CV*mean(r))/mean(r))^2))
    ri = rlnorm(N, log(mean(r)), logSD)
  } 
  #uniform r   
  if(user=="No") ri = exp(runif(N, log(r[1]), log(r[2])))  
  
  ## get N values between k[1] and k[2], assing to ki
  ki = exp(runif(N, log(k[1]), log(k[2])))  
  
  ## assign ri, ki, and final biomass range to itheta
  itheta=cbind(r=ri,k=ki, lam1=lambda[1],lam2=lambda[2], sigR=sigR) 
  
  ## call Schaefer function with parameters in itheta
  M = apply(itheta,1,.schaefer) 
  i=1:N
  ## prototype objective function
  get.ell=function(i) M[[i]]$ell
  ell = sapply(i, get.ell) 
  return(list(r=ri,k=ki, ell=ell))	
})
}




#2.1. Sensitivity Analysis

#Put sensitivity data in list
r.Sens=c(r.scen1,r.scen2)
r.Sens.sd=c(r.scen1.sd,r.scen2.sd)
Catch.matrix=matrix(c(lower.Total.catch.tons,median.Total.catch.tons,
                      upper.Total.catch.tons),ncol=3)

N.sens.test=length(K.UPPER)*length(r.Sens)*ncol(Catch.matrix)
Sens.mat=as.data.frame(expand.grid(K.UPPER,r.Sens,1:ncol(Catch.matrix)))
colnames(Sens.mat)=c("K.upper","r","dummy catch")
Sens.mat$r.sd=with(Sens.mat,ifelse(r==r.scen1,r.scen1.sd,r.scen2.sd))


List.sens=vector('list',length=N.sens.test)
N.sens=length(List.sens)
for (i in 1:N.sens)
{

  List.sens[[i]]$K.upper=Sens.mat[i,1]
  List.sens[[i]]$r=Sens.mat[i,2]
  List.sens[[i]]$r.sd=Sens.mat[i,4]
  List.sens[[i]]$catch=Catch.matrix[,Sens.mat[i,3]]

}


#Run sensitivity analysis
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')

setwd(handl_OneDrive("Analyses/Demography/White shark/Outputs/MSY_Catch/Sensitivities"))


STORE.sens=vector("list",length=N.sens)

system.time(for (i in 1:N.sens)
{
  K.upper=List.sens[[i]]$K.upper
  r=List.sens[[i]]$r
  R.SD=List.sens[[i]]$r.sd
  KTCH=List.sens[[i]]$catch
  
  K.MAX=max(KTCH)
  
  #Parameters
  Sp=paste("sensitivity",i)
  
  
  cdat=data.frame(yr=yr,ct=KTCH,stock=rep(Sp,yr.n),res=rep("Very low",yr.n))

    r.list <- list(quantile(rnorm(10000,r,R.SD),probs=c(0.015,0.975))) 
  k.list <- list(c(K.MAX,K.upper*K.MAX))
  
  
  
  ## Loop through stocks
  stock_id <- unique(as.character(cdat$stock))
  stock=stock_id
    p=1
    
    yr   <- cdat$yr[as.character(cdat$stock)==stock]
    ct   <- as.numeric(cdat$ct[as.character(cdat$stock)==stock])
    
    #sample the catch
    # ct=mapply(SamplE,median.Total.catch.tons,SD.Total.catch.tons)
    
    res  <- unique(as.character(cdat$res[as.character(cdat$stock)==stock])) ## resilience 
    nyr  <- length(yr)    ## number of years in the time series
    
    #  cat("\n","Stock",stock,"\n")
    #  flush.console()
    
    ## PARAM SECTION
    #user defined pars
    if(user=="Yes")
    {
      start_r=r.list[[p]]
      start_k=k.list[[p]]
      startbio=startbio.list[[p]]
      finalbio=finalbio.list[[p]]
      
    }else
      #defaults
    {
      start_r  <- if(res == "Very low")             ## for unknown resilience 
      {c(0.015, 0.1)}else if
      (res == "Low")
      {c(0.05,0.5)}else if
      (res == "High")
      {c(0.6,1.5)}else
      {c(0.2,1)} 
      
      start_k     <- c(max(ct),100*max(ct))          ## for unknown upper K 
      
      startbio    <- if(ct[1]/max(ct) < 0.5) {c(0.5,0.9)} else {c(0.3,0.6)} # for unknown startbio
      finalbio    <- if(ct[nyr]/max(ct) > 0.5) {c(0.3,0.7)} else {c(0.01,0.4)} # for unknown finalbio 
      
    }
    
    interyr   <- yr[2]   ## interim year within time series for which biomass estimate is available;
    ## set to yr[2] if no estimates are available
    interbio   <- c(0, 1) ## biomass range for interim year, as fraction of k;
    ## set to 0 and 1 if not available
    
    startbt     <- seq(startbio[1], startbio[2], by = 0.05) ## apply range of start biomass in steps of 0.05	
    parbound <- list(r = start_r, k = start_k, lambda = finalbio, sigR)
    
    
    #print out some data and assumption (checks)
    cat("Last year =",max(yr),", last catch (tons)=",ct[nyr],"\n")
    cat("Resilience =",res,"\n")
    cat("Process error =", sigR,"\n")
    cat("Assumed initial biomass (B/k) =", startbio[1],"-", startbio[2], " k","\n")
    cat("Assumed intermediate biomass (B/k) in", interyr, " =", interbio[1],"-",interbio[2]," k","\n")
    cat("Assumed final biomass (B/k) =", parbound$lambda[1],"-",parbound$lambda[2]," k","\n")
    cat("Initial bounds for r =", parbound$r[1], "-", parbound$r[2],"\n")
    cat("Initial bounds for k (tons)=", format(parbound$k[1], digits=3), "-", format(parbound$k[2],digits=3),"\n")
    
    #flush.console()
    
    
    
    
    ## MAIN
    R1 = sraMSY(parbound, n)  
    
    
    ## Get statistics on r, k, MSY and determine new bounds for r and k
    r1 	<- R1$r[R1$ell==1]
    k1 	<- R1$k[R1$ell==1]
    msy1  <- r1*k1/4
    mean_msy1 <- exp(mean(log(msy1))) 
    
    if(user=="Yes") max_k1 <- max(k1[r1*k1/4<mean_msy1])                            #REVIEW!!!!!
    
    if(!user=="Yes")
    {
      max_k1a  <- min(k1[r1<1.1*parbound$r[1]]) ## smallest k1 near initial lower bound of r
      max_k1b  <- max(k1[r1*k1/4<mean_msy1]) ## largest k1 that gives mean MSY
      max_k1 <- if(max_k1a < max_k1b) {max_k1a} else {max_k1b}
    }
    
    if(length(r1)<10)
    {
      cat("Too few (", length(r1), ") possible r-k combinations, check input parameters","\n")
      #   flush.console()
    }
    
    if(length(r1)>=10)
    {
      ## set new upper bound of r to 1.2 max r1
      parbound$r[2] <- 1.2*max(r1)
      ## set new lower bound for k to 0.9 min k1 and upper bound to max_k1 
      parbound$k 	  <- c(0.9 * min(k1), max_k1)
      
      
      cat("First MSY =", format(mean_msy1, digits=3),"\n")
      cat("First r =", format(exp(mean(log(r1))), digits=3),"\n")
      cat("New upper bound for r =", format(parbound$r[2],digits=2),"\n")	
      cat("New range for k (tons)=", format(parbound$k[1], digits=3), "-", format(parbound$k[2],digits=3),"\n")
      
      
      ## Repeat analysis with new r-k bounds
      R1 = sraMSY(parbound, n)
      
      ## Get statistics on r, k and msy
      r = R1$r[R1$ell==1]
      k = R1$k[R1$ell==1]
      k = k #convert to tonnes
      msy = r * k / 4
      mean_ln_msy = mean(log(msy))
      
      ct1=ct #convert back to tonnes
      
      ## plotting
      Mean.MSY=exp(mean(log(msy)))
  
      tiff(file=paste("model.output/CatchMSY_Multi.",Sp[p],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
      par(mfcol=c(2,3))    
      plot(yr, ct1, type="l", ylim = c(0, max(ct1)), xlab = "Year", ylab = "Catch (t)", main = stock,lwd=2)
      abline(h=exp(mean(log(msy))),col="black", lwd=2)
      abline(h=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
      abline(h=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
      
      
      hist(r, freq=F, xlim=c(0, 1.2 * max(r)), main = "")
      abline(v=exp(mean(log(r))),col="black",lwd=2)
      abline(v=exp(mean(log(r))-2*sd(log(r))),col="red")
      abline(v=exp(mean(log(r))+2*sd(log(r))),col="red")
      
      plot(r1, k1, xlim = start_r, ylim = start_k, xlab="r", ylab="k (t)")
      
      hist(k, freq=F, xlim=c(0, 1.2 * max(k)), xlab="k (t)", main = "")
      abline(v=exp(mean(log(k))),col="black", lwd=2)	
      abline(v=exp(mean(log(k))-2*sd(log(k))),col="red")
      abline(v=exp(mean(log(k))+2*sd(log(k))),col="red")
      
      
      plot(log(r), log(k),xlab="ln(r)",ylab="ln(k)")
      abline(v=mean(log(r)))
      abline(h=mean(log(k)))
      abline(mean(log(msy))+log(4),-1, col="black",lwd=2)
      abline(mean(log(msy))-2*sd(log(msy))+log(4),-1, col="red")
      abline(mean(log(msy))+2*sd(log(msy))+log(4),-1, col="red")
      
      hist(msy, freq=F, xlim=c(0, 1.2 * max(msy)), xlab="MSY (t)",main = "")
      abline(v=exp(mean(log(msy))),col="black", lwd=2)
      abline(v=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
      abline(v=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
      dev.off()
      
      cat("Possible combinations = ", length(r),"\n")
      cat("geom. mean r =", format(exp(mean(log(r))),digits=3), "\n")
      cat("r +/- 2 SD =", format(exp(mean(log(r))-2*sd(log(r))),digits=3),"-",format(exp(mean(log(r))+2*sd(log(r))),digits=3), "\n")
      cat("geom. mean k (tons)=", format(exp(mean(log(k))),digits=3), "\n")
      cat("k +/- 2 SD (tons)=", format(exp(mean(log(k))-2*sd(log(k))),digits=3),"-",format(exp(mean(log(k))+2*sd(log(k))),digits=3), "\n")
      cat("geom. mean MSY (tons)=", format(exp(mean(log(msy))),digits=3),"\n")
      cat("MSY +/- 2 SD (tons)=", format(exp(mean_ln_msy - 2 * sd(log(msy))),digits=3), "-", format(exp(mean_ln_msy + 2 * sd(log(msy))),digits=3), "\n")
      
      
      ## Export results
      outfile  <- paste("model.output/CatchMSY_Output.",Sp[p],".csv",sep="")
      output = data.frame(Species=stock, Proc.Err=sigR, B.init.low=startbio[1], B.init.high=startbio[2], 
                          B.2.low=interbio[1], B.2.high=interbio[2], 
                          B.fin.low=finalbio[1], B.fin.high=finalbio[2], Min.yr=min(yr), Max.yr=max(yr), 
                          Resilience=res, Catch.max=max(ct), Catch.init=ct[1], Catch.fin=ct[nyr], 
                          Lenth.r=length(r), Mean.r=exp(mean(log(r))), SD.r=sd(log(r)), Min.r=min(r), CI.95.low.r=quantile(r,0.025), 
                          Median.r=median(r), CI.95.up.r=quantile(r,0.975), Max.r=max(r), 
                          Mean.k=exp(mean(log(k))), SD.k=sd(log(k)), Min.k=min(k), CI.95.low.k=quantile(k, 0.025), 
                          Median.k=median(k), CI.95.up.k=quantile(k, 0.975),Max.k=max(k),
                          Median.MSY=median(msy),Mean.MSY=Mean.MSY,SD.MSY=sd(log(msy)),Min.MSY=min(msy),Max.MSY=max(msy),
                          CI.95.low.MSY=quantile(msy,0.025),CI.95.up.MSY=quantile(msy,0.975)) 
      
      
       
      write.table(output, file = outfile, sep = ",",row.names = FALSE, col.names =T)
      
    }
  
  STORE.sens[[i]]=list(msy=msy,r=r,k=k)
  
  })  
  
 


#REPORT
LWD=2.5
Year.Vec <- c(yr, tail(yr, 1), rev(yr), yr[1]) #year vector for polygons
Biom.Vec <- c(lower.Total.catch.tons, tail(upper.Total.catch.tons, 1), 
              rev(upper.Total.catch.tons), lower.Total.catch.tons[1])


#Reconstructed catches
tiff(file="figure1.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mai=c(.9,.9,.1,.1),las=1)
plot(yr,median.Total.catch.tons,type='l',lwd=LWD,xlab="Year",ylab="Catch (tons)",
     cex.axis=1.35,cex.lab=1.75,ylim=c(0,max(upper.Total.catch.tons)*1.05))
polygon(Year.Vec, Biom.Vec, col = "grey60", border = NA)
lines(yr,median.Total.catch.tons,lwd=LWD)
dev.off()

#MSY

STORE.msy=as.data.frame(matrix(NA,nrow=N.sens,ncol=4))
colnames(STORE.msy)=c('lower','mean','upper','Sensitivity')
for (i in 1:N.sens)
{
  emeswhy=STORE.sens[[i]]$msy
  mean_ln_msy = mean(log(emeswhy))
  STORE.msy$mean[i]=exp(mean(log(emeswhy)))
  STORE.msy$lower[i]=exp(mean_ln_msy - 2 * sd(log(emeswhy)))
  STORE.msy$upper[i]=exp(mean_ln_msy + 2 * sd(log(emeswhy)))
  STORE.msy$Sensitivity[i]=paste("k.upper=",List.sens[[i]]$K.upper,"r=",List.sens[[i]]$r)
}
write.csv(STORE.msy,"STORE.msy.csv")

UPPER.y=max(upper.Total.catch.tons)*1.2
tiff(file="figure2.Lower.K.upper.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(2,3),mai=c(.1,.1,.1,.1),oma=c(4,6,2,.01),las=1,mgp=c(1,1.05,0))
for (i in c(1,4,7,10,13,16))
{
  emeswhy=STORE.sens[[i]]$msy
  mean_ln_msy = mean(log(emeswhy))
  ktch=List.sens[[i]]$catch
  
  MED.msy=exp(mean(log(emeswhy)))
  LOW.msy=exp(mean_ln_msy - 2 * sd(log(emeswhy)))
  UP.msy=exp(mean_ln_msy + 2 * sd(log(emeswhy)))
  
  plot(yr, ktch, type="l", ylim = c(0, UPPER.y), 
       ylab="",xlab="",xaxt="n",yaxt="n",cex.main=1.5,
       main = "",lwd=3)
  
  axis(1,at=yr,labels=F,tck=-0.015)
  axis(2,at=seq(0, UPPER.y,by=20),labels=F,tck=-0.04)
  axis(1,at=seq(yr[1],yr[length(yr)],5),labels=F,tck=-0.03,cex.axis=1.25)
  axis(1,at=seq(yr[6],yr[length(yr)],10),labels=F,tck=-0.06,cex.axis=1.6)
  
  if(i%in%c(4,10,16))axis(1,at=seq(yr[6],yr[length(yr)],10),labels=seq(yr[6],
                    yr[length(yr)],10),tck=-0.06,cex.axis=1.5)
  if(i%in%c(1,4))axis(2,at=seq(0, UPPER.y,by=20),
                labels=seq(0, UPPER.y,by=20),cex.axis=1.5,tck=-0.04)
  
  
  abline(h=MED.msy,col="grey60", lwd=3)
  abline(h=LOW.msy,col="grey60", lwd=3,lty=2)
  abline(h=UP.msy,col="grey60", lwd=3,lty=2)
  
  if(i==1)mtext("Low CI catch",side=3,line=.5,font=1,las=1,cex=1.45,col='black')
  if(i==7)mtext("Mean catch",side=3,line=.5,font=1,las=1,cex=1.45,col='black')
  if(i==13)mtext("Up CI catch",side=3,line=.5,font=1,las=1,cex=1.45,col='black')
  if(i==1)mtext("LH1",side=2,line=2.8,font=1,las=1,cex=1.45,col='black')
  if(i==4)mtext("LH2",side=2,line=2.8,font=1,las=1,cex=1.45,col='black')
  
  
}
mtext("Year",1,3,outer=T,cex=1.75)
mtext("Catch (t)",2,4,outer=T,las=3,cex=1.75)
dev.off()


tiff(file="figure3.Mid.K.upper.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(2,3),mai=c(.1,.1,.1,.1),oma=c(4,6,2,.01),las=1,mgp=c(1,1.05,0))
for (i in c(2,5,8,11,14,17))
{
  emeswhy=STORE.sens[[i]]$msy
  mean_ln_msy = mean(log(emeswhy))
  ktch=List.sens[[i]]$catch
  
  MED.msy=exp(mean(log(emeswhy)))
  LOW.msy=exp(mean_ln_msy - 2 * sd(log(emeswhy)))
  UP.msy=exp(mean_ln_msy + 2 * sd(log(emeswhy)))
  
  
  plot(yr, ktch, type="l", ylim = c(0, UPPER.y), 
       ylab="",xlab="",xaxt="n",yaxt="n",cex.main=1.5,
       main = "",lwd=3)
  
  axis(1,at=yr,labels=F,tck=-0.015)
  axis(2,at=seq(0, UPPER.y,by=20),labels=F,tck=-0.04)
  axis(1,at=seq(yr[1],yr[length(yr)],5),labels=F,tck=-0.03,cex.axis=1.25)
  axis(1,at=seq(yr[6],yr[length(yr)],10),labels=F,tck=-0.06,cex.axis=1.6)
  
  if(i%in%c(5,11,17))axis(1,at=seq(yr[6],yr[length(yr)],10),labels=seq(yr[6],
                        yr[length(yr)],10),tck=-0.06,cex.axis=1.5)
  if(i%in%c(2,5))axis(2,at=seq(0, UPPER.y,by=20),
                      labels=seq(0, UPPER.y,by=20),cex.axis=1.5,tck=-0.04)
  
  
  abline(h=MED.msy,col="grey60", lwd=3)
  abline(h=LOW.msy,col="grey60", lwd=3,lty=2)
  abline(h=UP.msy,col="grey60", lwd=3,lty=2)
  
  if(i==2)mtext("Low CI catch",side=3,line=.5,font=1,las=1,cex=1.45,col='black')
  if(i==8)mtext("Mean catch",side=3,line=.5,font=1,las=1,cex=1.45,col='black')
  if(i==14)mtext("Up CI catch",side=3,line=.5,font=1,las=1,cex=1.45,col='black')
  if(i==2)mtext("LH1",side=2,line=2.8,font=1,las=1,cex=1.45,col='black')
  if(i==5)mtext("LH2",side=2,line=2.8,font=1,las=1,cex=1.45,col='black')
  
  
}
mtext("Year",1,3,outer=T,cex=1.75)
mtext("Catch (t)",2,4,outer=T,las=3,cex=1.75)
dev.off()


tiff(file="figure4.High.K.upper.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(2,3),mai=c(.1,.1,.1,.1),oma=c(4,6,2,.01),las=1,mgp=c(1,1.05,0))
for (i in c(3,6,9,12,15,18))
{
  emeswhy=STORE.sens[[i]]$msy
  mean_ln_msy = mean(log(emeswhy))
  ktch=List.sens[[i]]$catch
  
  MED.msy=exp(mean(log(emeswhy)))
  LOW.msy=exp(mean_ln_msy - 2 * sd(log(emeswhy)))
  UP.msy=exp(mean_ln_msy + 2 * sd(log(emeswhy)))
  
  plot(yr, ktch, type="l", ylim = c(0, UPPER.y), 
       ylab="",xlab="",xaxt="n",yaxt="n",cex.main=1.5,
       main = "",lwd=3)
  
  axis(1,at=yr,labels=F,tck=-0.015)
  axis(2,at=seq(0, UPPER.y,by=20),labels=F,tck=-0.04)
  axis(1,at=seq(yr[1],yr[length(yr)],5),labels=F,tck=-0.03,cex.axis=1.25)
  axis(1,at=seq(yr[6],yr[length(yr)],10),labels=F,tck=-0.06,cex.axis=1.6)
  
  if(i%in%c(6,12,18))axis(1,at=seq(yr[6],yr[length(yr)],10),labels=seq(yr[6],
                        yr[length(yr)],10),tck=-0.06,cex.axis=1.5)
  if(i%in%c(3,6))axis(2,at=seq(0, UPPER.y,by=20),
                      labels=seq(0, UPPER.y,by=20),cex.axis=1.5,tck=-0.04)
  
  
  abline(h=MED.msy,col="grey60", lwd=3)
  abline(h=LOW.msy,col="grey60", lwd=3,lty=2)
  abline(h=UP.msy,col="grey60", lwd=3,lty=2)
  
  if(i==3)mtext("Low CI catch",side=3,line=.5,font=1,las=1,cex=1.45,col='black')
  if(i==9)mtext("Mean catch",side=3,line=.5,font=1,las=1,cex=1.45,col='black')
  if(i==15)mtext("Up CI catch",side=3,line=.5,font=1,las=1,cex=1.45,col='black')
  if(i==3)mtext("LH1",side=2,line=2.8,font=1,las=1,cex=1.45,col='black')
  if(i==6)mtext("LH2",side=2,line=2.8,font=1,las=1,cex=1.45,col='black')
  
  
}
mtext("Year",1,3,outer=T,cex=1.75)
mtext("Catch (t)",2,4,outer=T,las=3,cex=1.75)
dev.off()
