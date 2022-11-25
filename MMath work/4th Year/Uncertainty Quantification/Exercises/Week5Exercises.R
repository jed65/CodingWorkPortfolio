#### Snippet from notes ####
library(mvtnorm) #R package for simulating multivariate Gaussians
calcSigma <- function(X1,X2,theta) { # Function takes in a discretised grid
  #input X1 x X2 and computes covariance matrix using covariance function k
  v1=theta[1]
  v2=theta[2]
  alpha=theta[3]
  r=theta[4]
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1)) #fill with zeroes, give rows equal to number of elements of x1
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- v1*exp(-(abs(X1[i]-X2[j])/r)^alpha)+v2*(i==j) #use formula given, remember boolean are 1 or 0 in R
    }
  }
  return(Sigma) 
}
plotGP=function(theta=c(1,.01,2,2))
{
  f=data.frame(x=c(0.9,3.8,5.2,6.1,7.5,9.6), y=c(0.1,1.2,2.1,1.1,1.5,1.2)) #this is the data given
  par(mfrow=c(1,4)) #this command creates a subplot grid
  plot(f$x,f$y,col=2,pch=10,ylim=c(-2.5,2.5),xlab="",ylab="") #plot data points on their own
  x.star=seq(0,10,length=50) #this is the vector of points that we're interested in
  Sig=calcSigma(x.star,x.star,theta)
  prior.sample=matrix(0,nrow=length(x.star),ncol=20)
  for(i in 1:ncol(prior.sample))
  {
    prior.sample[,i]=rmvnorm(1,mean=rep(0,length(x.star)), sigma=Sig) #this is prior distribution based off this covariance function
  }
  matplot(x.star,prior.sample,type="l",xlab="",ylab="",
          main="prior sample",col="orange") #plot 2/4
  x=f$x 
  k.xx <- calcSigma(x,x,theta) #get covariance matrix ingredients to form posterior
  k.xxs <- calcSigma(x,x.star,theta)
  k.xsx <- calcSigma(x.star,x,theta)
  k.xsxs <- calcSigma(x.star,x.star,theta)
  n.samples=20 #number of samples 
  # The standard deviation of the noise
  sigma.n=theta[2]
  # calculate the mean and covariance functions
  f.star.bar <- k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%f$y
  cov.f.star <- k.xsxs - k.xsx%*%solve(k.xx +
                                         sigma.n^2*diag(1, ncol(k.xx)))%*%k.xxs #calculate updated mean/covariance
  # calculate the sample functions
  values <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples) 
  for (i in 1:n.samples) {
    values[,i] <- rmvnorm(1, f.star.bar, cov.f.star) #generate simulations from posterior
  }

  matplot(x.star,values,type="l",xlab="",ylab="",
          main="posterior sample",xaxs="i",yaxs="i",col="orange") #plot 3/4
  plot(x.star,f.star.bar,type="l",xlab="",ylab="",
       main="posterior mean+confidence band",col=4,ylim=c(-2.5,2.5))
  y.up=f.star.bar+1.96*diag(cov.f.star)^0.5
  y.low=f.star.bar-1.96*diag(cov.f.star)^0.5 #get confidence bands and plot them (4/4)
  points(f$x,f$y,col=2)
  lines(x.star,y.low,col="dark grey",lty=2)
  lines(x.star,y.up,col="dark grey",lty=2)
}
plotGP()
#### Ex.1 ####
xstar<-seq(0,10,length=50) #50 points in [0,10] equally spread apart
calcCov<-function(x1,x2){ #this function gives us covariance matrices
  CovarianceMatrix<-matrix(rep(0,length(x1)*length(x2)),nrow=length(x1))
  for (i in 1:nrow(CovarianceMatrix)){
    for (j in 1:ncol(CovarianceMatrix)){
      CovarianceMatrix[i,j]<-exp(-0.5*(x1[i]-x2[j])^2)
    }
  }
  return(CovarianceMatrix)
}
D=data.frame(x=c(0.9,3.8,5.2,6.1,7.5,9.6),y=c(0.1,1.2,2.1,1.1,1.5,1.2)) #data given
x<-D$x
y<-D$y
k.xstarx<-calcCov(xstar,x)
k.xx<-calcCov(x,x)
noiseSigma<-0.15
mstar<-k.xstarx %*% solve(k.xx+noiseSigma^2*diag(1,ncol(k.xx))) %*% y #get posterior mean
par(mfrow=c(1,1))
plot(D$x,D$y,col=2,pch=10,ylim=c(-2.5,2.5),xlab="",ylab="") #plot data points and mean curve
lines(xstar,mstar,col="blue")
k.xxstar<-calcCov(x,xstar)
k.xstarxstar<-calcCov(xstar,xstar)
covstar<-k.xstarxstar-k.xstarx  %*% solve(k.xx+noiseSigma^2*diag(1,ncol(k.xx))) %*% k.xxstar #get posterior covariance
fupper<-mstar+1.96*diag(covstar)^(0.5)
flower<-mstar-1.96*diag(covstar)^(0.5)
lines(xstar,fupper,col="dark grey",lty=2)
lines(xstar,flower,col="dark grey",lty=2) #add the upper and lower bounds
xnew<-seq(10,11.5,length=10) #create vector for new inputs
k.xnewx<-calcCov(xnew,x)
k.xnewxnew<-calcCov(xnew,xnew)
k.xxnew<-calcCov(x,xnew)
mstarnew<-k.xnewx %*% solve(k.xx+noiseSigma^2*diag(1,ncol(k.xx))) %*% y #get new posterior mean
covstarnew<-k.xnewxnew-k.xnewx  %*% solve(k.xx+noiseSigma^2*diag(1,ncol(k.xx))) %*% k.xxnew #new posterior covariance
fuppernew<-mstarnew+1.96*diag(covstarnew)^(0.5) #upper bound
flowernew<-mstarnew-1.96*diag(covstarnew)^(0.5) #lower bound
plot(D$x,D$y,col=2,pch=10,ylim=c(-2.5,2.5),xlim=c(0,12),xlab="",ylab="") #plot, include extrapolatory region
lines(xstar,mstar,col="blue")
lines(xnew,mstarnew,col="orange")
lines(xstar,fupper,col="dark grey",lty=2)
lines(xstar,flower,col="dark grey",lty=2)
lines(xnew,flowernew,col="green",lty=2)
lines(xnew,fuppernew,col="green",lty=2)
calcCov2<-function(x1,x2,theta){ #this function gives us covariance matrices where we have general hyperparameters
  CovarianceMatrix<-matrix(rep(0,length(x1)*length(x2)),nrow=length(x1))
  r=theta[1]
  alpha=theta[2]
  for (i in 1:nrow(CovarianceMatrix)){
    for (j in 1:ncol(CovarianceMatrix)){
      CovarianceMatrix[i,j]<-exp((-(x1[i]-x2[j])/r)^alpha)
    }
  }
  return(CovarianceMatrix)
}
#### Ex.2 ####
library(deSolve)
params<-c(alpha=2/3,beta=4/3,delta=1,gamma=1) #define parameters of LV model
initial<-c(X=0.9,Y=1.5) #initial conditions
LV<-function(t,initial,params){ #function to give back list of dx/dy
  with(as.list(c(initial, params)),{
  dX<-params[1]*X-params[2]*X*Y
  dY<-params[3]*X*Y-params[4]*Y
  list(c(dX,dY))
  })
}
times<-seq(0,20,by=5)
PredPreyData<-as.data.frame(ode(y=initial,times=times,func=LV,parms=params)) #make dataframe with solutions
colnames(PredPreyData)=c("t","X","Y")
PredPreyData$noisyX<-PredPreyData$X+rnorm(length(times),0,0.15) #add noise to predator
PredPreyData$noisyY<-PredPreyData$Y+rnorm(length(times),0,0.2) #add noise to prey
plot(PredPreyData$noisyX~PredPreyData$t,data=PredPreyData,col="blue",cex=0.5,ylim=c(-3,3),xlab="Time",ylab="Population") #plot predators
points(PredPreyData$noisyY~PredPreyData$t,data=PredPreyData,col="red",cex=0.5) #plot prey
#Need to predict the parameter values theta
MarginalLklh<-function(par,y){ #form marginal log likelihood
  times<-seq(0,20,by=5)
  K<-matrix(rep(0,length(times)*length(times)),nrow=length(times))
  for (i in 1:nrow(K)){
    for (j in 1:nrow(K)){
      K[i,j]=exp(-((times[i]-times[j])/par[1])^par[2])+0.15*(i==j)
    }
  }
  marg<-0.5*y %*% solve(K) %*% y+0.5*log(det(K))+0.5*length(y)*log(2*pi)
  return(marg)
}
thetaMax<-optim(par=c(2,2),MarginalLklh,y=PredPreyData$noisyX)

calcCovariance <- function(X1,X2,theta) { # Function takes in a discretised grid
  #input X1 x X2 and computes covariance matrix using covariance function k
  r=theta[1]
  alpha=theta[2]
  noiseVar=theta[3]
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1)) #fill with zeroes, give rows equal to number of elements of x1
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-((X1[i]-X2[j])/r)^alpha)+noiseVar*(i==j) #use formula given, remember boolean are 1 or 0 in R
    }
  }
  return(Sigma) 
}
tstar<-seq(0,20,by=0.1)
theta<-c(thetaMax$par[1],thetaMax$par[2],0.15)
k.tstar_t<-calcCovariance(tstar,times,theta)
k.tt<-calcCovariance(times,times,theta)
k.tstar_tstar<-calcCovariance(tstar,tstar,theta)
k.t_tstar<-calcCovariance(times,tstar,theta)
m_starPredator<-k.tstar_t%*%solve(k.tt+0.15*diag(1,ncol(k.tt)))%*%PredPreyData$noisyX
cov_starPredator<-k.tstar_tstar-k.tstar_t %*% solve(k.tt+0.15*diag(1,ncol(k.tt))) %*% k.t_tstar
upperCIPredator<-m_starPredator+1.96*diag(cov_starPredator)^(0.5)
lowerCIPredator<-m_starPredator-1.96*diag(cov_starPredator)^(0.5)
lines(tstar,m_starPredator,col="black")
lines(tstar,lowerCIPredator,col="dark grey",lty=2)
lines(tstar,upperCIPredator,col="dark grey",lty=2)
tnew<-seq(20,25,length.out=10)
k.tnew_t<-calcCovariance(tnew,times,theta)
k.tnew_tnew<-calcCovariance(tnew,tnew,theta)
k.t_tnew<-calcCovariance(times,tnew,theta)
k.tt<-calcCovariance(times,times,theta)
m_starNewPredator<-k.tnew_t%*%solve(k.tt+0.15*diag(1,ncol(k.tt)))%*%PredPreyData$noisyX
cov_starNewPredator<-k.tnew_tnew-k.tnew_t %*% solve(k.tt+0.15*diag(1,ncol(k.tt))) %*% k.t_tnew
upperCIPredatorNew<-m_starNewPredator+1.96*diag(cov_starNewPredator)^(0.5)
lowerCIPredatorNew<-m_starNewPredator-1.96*diag(cov_starNewPredator)^(0.5)
plot(PredPreyData$noisyX~PredPreyData$t,data=PredPreyData,col="blue",cex=0.5,ylim=c(-4,4),xlim=c(0,25),xlab="Time",ylab="Population") #plot predators
lines(tnew,m_starNewPredator,col="green")
lines(tnew,lowerCIPredatorNew,col="orange",lty=2)
lines(tnew,upperCIPredatorNew,col="orange",lty=2)