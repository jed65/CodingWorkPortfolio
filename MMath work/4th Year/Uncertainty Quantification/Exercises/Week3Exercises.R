#### Ex.1 Section 1: Inverse UQ ####
ODEsol=function(t,params) #function that finds the values fx solving the ODE using analytical solution
{
  c=params[1]; r=params[2]; a=params[3];
  fx=c/(1+a*exp(-r*t))
  return(list(fx))
}
tt<-1:10 #vector of days we're interested in
init.params<-c(c=100,r=2,a=49) #third parameter a equivalent to setting f(0)=2
sim.data=as.data.frame(ODEsol(tt,init.params)) #make dataframe from output 
set.seed(12345) #set seed for reproducibility
sim.data$noisy=sim.data[,1]+rnorm(nrow(sim.data),0,20) #add noise column
colnames(sim.data)=c("growth","noisy.growth") #name the columns for future use
param.start=list(c=max(sim.data$noisy.growth),r=2) #new start parameters
nls.fit=nls(noisy.growth~I(c/(1+((c-2)/2)*exp(-r*tt))),data=sim.data,start=param.start) #this finds parameter estimates for model given (RHS equation)
coef(nls.fit) #next two lines give value of coefficients c and r and variances
vcov(nls.fit)
nls.estimate=fitted(nls.fit) #use the estimated parameters to fit the data
plot(noisy.growth~tt,data=sim.data,xlab="Day",ylab="Growth") #rest of code produces plot
lines(noisy.growth~tt,data=sim.data)
lines(growth~tt,data=sim.data,col="blue")
lines(tt,nls.estimate,col="red")
legend('bottomright',legend=c("Data","ODE Solution","NLS Fitted"),col=c("black","blue","red"),pt.lwd=1,lwd=1,bty="n",text.font=.5)
#### end ####
#### Ex.1 Section 2: Forward UQ ####
c.hat=coef(nls.fit)[1] #extract means/variances of fitted c and r
r.hat=coef(nls.fit)[2]
sd.chat=sqrt(vcov(nls.fit)[1,1])
sd.rhat=sqrt(vcov(nls.fit)[2,2])
set.seed(54321)
c.sample=rnorm(20,mean=c.hat,sd=sd.chat) #get sample from asymptotic distributions of c/r
r.sample=rnorm(20,mean=r.hat,sd=sd.rhat)
sols=matrix(0,nrow=10,ncol=20) #create matrix of zeros with 10 rows (no. of days) and 20 columns (no. of samples)
for (i in 1:20)
{
  sols[,i]=unlist(ODEsol(tt,c(c.sample[i],r.sample[i],(c.sample[i]-2)/2))) #unlist function creates vector or matrix from list
}
matplot(sols,col="orange",type="l",ylim=c(1,115),xlab="Day",ylab="Growth")
points(noisy.growth~tt,data=sim.data)
lines(growth~tt,data=sim.data,col="blue")
lines(tt,nls.estimate,col="red")
#### end ####
#### Ex.2 ####
conc<-c(2.857,5.005,7.519,22.102,27.770,39.198,45.483,203.784)
rate<-c(14.583,24.741,31.346,72.970,77.501,96.088,96.966,108.884)
solution=matrix(0,nrow=length(rate),ncol=2)
solution[,1]=conc 
solution[,2]=rate #matrix of values formed 
observed.data<-as.data.frame(solution)
colnames(observed.data)=c("Concentration","Rate")
param.begin=list(theta1=119.91,theta2=20.64) #parameters found by solving linear system from observed
nls.ratefit=nls(rate~I((theta1*conc)/(theta2+conc)),data=observed.data,start=param.begin) #fit model
coef(nls.ratefit)
vcov(nls.ratefit)
nls.rateestimate=fitted(nls.ratefit)
plot(rate~conc,data=observed.data,xlab="Concentration",ylab="Rate") #rest of code produces plot
lines(rate~conc,data=observed.data)
lines(conc,nls.rateestimate,col="red")
param.begin2=list(theta1=151.74,theta2=25.691) #changing values of the initial parameters doesn't change the outcome of the Gauss-Newton algorithm
nls.ratefit2=nls(rate~I((theta1*conc)/(theta2+conc)),data=observed.data,start=param.begin2) #fit model
coef(nls.ratefit2)
vcov(nls.ratefit2)
nls.rateestimate2=fitted(nls.ratefit2)
plot(rate~conc,data=observed.data,xlab="Concentration",ylab="Rate") #rest of code produces plot
lines(rate~conc,data=observed.data)
lines(conc,nls.rateestimate2,col="red")
observed.data$new.Rate=1/observed.data[,2]
observed.data$new.Concentration=1/observed.data[,1]
colnames(observed.data)=c("Concentration","Rate","new.Rate","new.Concentration")
linearFit<-lm(new.Rate~new.Concentration,data=observed.data)
summary(linearFit)
plot(new.Rate~new.Concentration,data=observed.data,xlab="1/Concentration",ylab="1/Rate")
lines(new.Rate~new.Concentration,data=observed.data)
lineEstimate=fitted(linearFit)
lines(observed.data$new.Concentration,lineEstimate,col="red")
RateSol=function(x,params)
{
  theta1=params[1];theta2=params[2];
  fxRate=(theta1*x)/(x+theta2)
  return(list(fxRate))
}
RateModel=as.data.frame(RateSol(observed.data$Concentration,c(148.24,26.04)))
RateModel[,2]=observed.data$Concentration
RateModel[,3]=observed.data$Rate
colnames(RateModel)=c("ModelRate","Concentration","ObservedRate")
plot(ObservedRate~Concentration,data=RateModel,xlab="Concentration",ylab="Rate")
points(ModelRate~Concentration,data=RateModel)
lines(ObservedRate~Concentration,data=RateModel)
lines(ModelRate~Concentration,data=RateModel,col="red")
pred1<-predict(nls.ratefit,list(conc=c(0.8))) #use if you just need predictions
library(propagate)
conc1<-c(2.857,5.005,7.519,22.102,27.770,39.198,45.483,203.784)
rate.pred<-predictNLS(nls.ratefit,newdata=data.frame(conc=conc1),interval="prediction",alpha=0.05,nsim=10000)$summary
matlines(conc1,rate.pred[,c("Sim.2.5%","Sim.97.5%")],col="blue",lty="solid") #this plots the interval
rate.pred[,c("Sim.2.5%","Sim.97.5%")] #gives the lower and upper bound for the point

library(numDeriv)
sum<-0
for (i in 1:length(nls.rateestimate))
{
  sum<-sum+(rate[i]-nls.rateestimate[i])^2
}
residsigma<-(1/(length(rate)-2))*sum
f1<-function(x,theta2){
  return((x*0.8)/(theta2+0.8))
  }
j0<-grad(func=f1,x=126.03243,theta2=17.07868)
f2<-function(x,theta2,y){
  return((x*y)/(theta2+y))
}
J0<-rep(0,8) #????
for (i in 1:8)
{
  J0[i]=grad(func=f2(x,theta2,rate[i]),x=126.03243,theta2=17.07868)
}


library(matlib)
halfwidth<-qt(0.975,6)*sqrt(residsigma*j0*inv(t(J0)%*%J0)*j0)
interval<-matrix(0,nrow=1,ncol=3)
interval[1]<-pred1-halfwidth
interval[2]<-pred1
interval[3]<-pred1+halfwidth
interval