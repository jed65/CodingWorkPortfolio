#### Ex.1 ####
Lklh<-function(x,n,mu){
  likelihood<-1
  for (i in 1:n){
    likelihood<-likelihood*(2*pi)^(-0.5)*exp(-0.5*(x[i]-mu)^2)
  }
  return(-likelihood)
}
n<-100
sampleHolder<-matrix(0,nrow=20,ncol=n)
for (j in 1:20){
  x<-rnorm(n,0.5,1)
  sampleHolder[j,]<-x
}
mu<-seq(-1,1,length.out=n)
Likelihood<-rep(0,length(mu))
LikelihoodHolder<-matrix(0,nrow=20,ncol=length(mu))
for (j in 1:20){
  for (i in 1:100){
    Likelihood[i]<-Lklh(sampleHolder[j,],n,mu[i])*(-1)
  }
  LikelihoodHolder[j,]<-Likelihood
  Likelihood<-rep(0,length(mu))
}
plot(mu,LikelihoodHolder[5,],col="blue",xlab="mu",ylab="L(mu,1)",cex=0.1)
lines(mu,LikelihoodHolder[1,])
lines(mu,LikelihoodHolder[2,])
lines(mu,LikelihoodHolder[3,])
lines(mu,LikelihoodHolder[4,])
lines(mu,LikelihoodHolder[5,])

mleHolder<-rep(0,5)
arithmeticMeanHolder<-rep(0,5)
for (j in 1:5){
  mle<-optim(par=c(0.5),Lklh,x=sampleHolder[j,],n=100)
  mleHolder[j]<-mle$par
  arithmeticMean<-sum(sampleHolder[j,])/n
  arithmeticMeanHolder[j]<-arithmeticMean
}
matplot(mu,t(LikelihoodHolder),type="l")

#### Ex.2 ####
logL<-function(x,n,mu,sigma){
  logLik<- (-1)*n*log(2*sigma)
  for (i in 1:n){
    logLik<-logLik-(1/(2*sigma))*abs(x[i]-mu)
  }
  return(logLik)
}
n<-50
x<-rLaplace(n,1.5,2)
mu<- -1
sigma<-seq(0.1,5,by=0.1)
Likelihood<-rep(0,length(sigma))
for (i in 1:length(sigma)){
  Likelihood[i]<-logL(x,n,mu,sigma[i])
}
plot(sigma,Likelihood,xlab="Sigma",ylab="log(Likelihood)",ylim=c(-250,0))
mu2<- seq(-1,1,by=0.2)
for (j in 1:length(mu2)){
  for (i in 1:length(sigma)){
    Likelihood[i]<-logL(x,n,mu2[j],sigma[i])
  }
  lines(sigma,Likelihood)
}
#seems that increasing mu tends to increase the mle value of sigma
#### Ex.3 ####
n<-seq(1,1000,by=20) #create vector of n values to test
Error<-rep(0,length(n)) #to hold the errors for each n
for (i in 1:length(n)){
  x<-rbinom(i,1,0.6) #generate random sample
  theta_m<-sum(x)/i #find arithmetic mean
  Error[i]<-abs(theta_m-0.6) #find absolute error
}
plot(n,Error,xlab="n",ylab="E(n)") #produce plot
nstar<-250
counter<-0
for (i in 1:10000){
  x<-rbinom(nstar,1,0.6)
  theta_m<-sum(x)/nstar
  Jtheta<-sum(x)/(theta_m^2)+(nstar-sum(x))/((1-theta_m)^2)
  lowCI<-theta_m-1.96/sqrt(Jtheta) #calculate lower bound
  upCI<-theta_m+1.96/sqrt(Jtheta) #calculate upper bound
  if (lowCI<0.6 && upCI>0.6){ #add 1 if 0.6 between bounds
    counter<-counter+1
  }
}
proportion<-counter/10000 #find proportion
proportion
#### Ex.4 ####
theta<-runif(1,0,1) #get random theta from uniform [0,1]
theta
n<-50
D<-rbinom(n,1,theta)
x<-seq(0,1,by=0.01)
BetaForPlot<-dbeta(x,1+sum(D),n-sum(D))
plot(BetaForPlot~x,col="blue",lty=1)
x[which(BetaForPlot==max(BetaForPlot))]
#Get a second interval (credible)
b<-qbeta(0.98,1+sum(D),n-sum(D)) #use in-built function at certain quantiles with difference 0.95
a<-qbeta(0.03,1+sum(D),n-sum(D))
FindInterval<-rbeta(10000,1+sum(D),n-sum(D)) #find interval from posterior distribution sample
FindInterval<-sort(FindInterval)
a<-FindInterval[10000*0.025]
b<-FindInterval[10000*0.975]
#### Ex.5 ####
n<-seq(20,60,by=20)
mu<-0
alpha<-4
beta<-5
supportx<-seq(0.1,2.1,by=0.01)
betas<-matrix(0,nrow=3,ncol=length(supportx))
for (i in n){
  x<-rnorm(i,mu,1) #tau=1
  sum<-0
  for (j in 1:i){
    sum<-sum+(x[j]-mu)^2
  }
  betanew<-beta+0.5*sum
  for (k in 1:length(supportx)){
    betas[i/20,k]<-dbeta(supportx[k],alpha,betanew)
  }
}
plot(supportx,betas[1,],col="blue",xlab="tau",ylab="Probability")
lines(supportx,betas[2,],col="red")
lines(supportx,betas[3,],col="orange")