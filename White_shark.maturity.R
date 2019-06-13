
#Maturity data (Francis 1996)
Fem.maturity=data.frame(Size.class=c("300-325","325-350","350-375","375-400","400-425","425-450"
                               ,"450-475","475-500","500-525","525-550","550-575"),
                        mid.point=c(313,338,363,388,413,438,463,488,513,538,563),
                        n.imm=c(1,0,3,5,3,0,1,1,0,0,0),n.mat=c(1,0,0,0,3,2,3,1,4,1,2))
Fem.maturity$prop.mat=with(Fem.maturity,n.mat/(n.imm +n.mat))
Fem.maturity$prop.mat=ifelse(is.na(Fem.maturity$prop.mat),0,Fem.maturity$prop.mat)


theta=c(L50=450,L95=500)
fun.mat=function(theta)
{
  Pred.mat=1*(1+exp((-log(19)*((Fem.maturity$mid.point-theta[1])/(theta[2]-theta[1])))))^-1
  #obj function
  epsilon=sum((Fem.maturity$prop.mat-Pred.mat)^2)
  
  return(epsilon)
}

#ad hoc (Terry Walker's approach)
fit=nlm(fun.mat, theta)
fit.preds=1*(1+exp((-log(19)*((Fem.maturity$mid.point-fit$estimate[1])/(fit$estimate[2]-fit$estimate[1])))))^-1 

#logistic glm approach
glm1 <-glm(prop.mat~mid.point,data=Fem.maturity,family=binomial)
pred.glm=predict(glm1,type ="response")


#get L50 and L95
lrPerc <-function(cf,p) (log(p/(1-p))-cf[1])/cf[2]
L50 <-lrPerc(coef(glm1),0.5)
L95 <-lrPerc(coef(glm1),0.95)


plot(Fem.maturity$mid.point,Fem.maturity$prop.mat)
lines(Fem.maturity$mid.point,fit.preds,col=2)
lines(Fem.maturity$mid.point,pred.glm,col=3)

lines(Fem.maturity$mid.point,pred.manual,col=4)

Pred.mat=exp(coef(glm1)[1]+coef(glm1)[2]*Fem.maturity$mid.point)/(1+exp(coef(glm1)[1]+coef(glm1)[2]*Fem.maturity$mid.point))



