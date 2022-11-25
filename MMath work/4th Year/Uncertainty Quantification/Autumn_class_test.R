#### Part 1 (up to 5(a)) ####
library(deSolve) #going to use ODE
n<-10 #fix as in the question
alphas<-seq(3,10,length.out=n) #get 10 equally spaced alpha in the given region
times<-seq(0,20,length.out=n) #get 10 equally spaced points in the domain
LV<-function(t,initial,params){ #function to give back list of dx/dy
  with(as.list(c(initial, params)),{
    dX<-params[1]*X-0.5*X*Y
    dY<-0.6*X*Y-Y
    list(c(dX,dY))
  })
}
initial<-c(X=2,Y=8) #initial values
PreyValues<-matrix(0,nrow=10,ncol=10) #define matrix to hold the values
QoI<-matrix(0,nrow=10,ncol=1)
for (i in 1:10){ 
  params<-c(alpha=alphas[i]) #this is the estimate for alpha
  sols<-as.data.frame(ode(y=initial,times=times,func=LV,parms=c(alpha=alphas[i]))) #solve system for this alpha
  PreyValues[,i]<-sols$X #assign column in this matrix to the solution
  QoI[i]<-PreyValues[10,i] #assign quantity of interest
}
d<-data.frame(alpha=alphas,eps=QoI) #create dataframe to store the data
plot(eps~alpha,data=d,ylim=c(-2,4)) #have a look at the data, this is also plot 1 and should be run with lines() given later

calcCovariance <- function(X1,X2) { # Function takes in a discretised grid
  #input X1 x X2 and computes covariance matrix using covariance function k
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1)) #fill with zeroes
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-0.5*(X1[i]-X2[j])^2)+0.15*(i==j) #use formula and parameters given
    }
  }
  return(Sigma) 
}
alphastar<-seq(0,10,by=0.1) #define region to interpolate/fit GP to
k.alphastar_alpha<-calcCovariance(alphastar,alphas) #obtain the quantities we need
k.alphaalpha<-calcCovariance(alphas,alphas)
k.alphastar_alphastar<-calcCovariance(alphastar,alphastar)
k.alpha_alphastar<-calcCovariance(alphas,alphastar)
m_star<-k.alphastar_alpha%*%solve(k.alphaalpha+0.15*diag(1,ncol(k.alphaalpha)))%*%d$eps #posterior mean and covariance
cov_star<-k.alphastar_alphastar-k.alphastar_alpha%*% solve(k.alphaalpha+0.15*diag(1,ncol(k.alphaalpha)))%*% k.alpha_alphastar
upper<-m_star+qnorm(0.995)*diag(cov_star)^(0.5) #99% confidence interval
lower<-m_star-qnorm(0.995)*diag(cov_star)^(0.5)
lines(alphastar,m_star,col="blue") #add these to the plot
lines(alphastar,lower,col="dark grey",lty=2)
lines(alphastar,upper,col="dark grey",lty=2)

#### Part 5b ####
alphanew<-c(5.5) #this is the point we want to predict at 
k.alphanew_alpha<-calcCovariance(alphanew,alphas) #obtain the quantities we need
k.alphaalpha<-calcCovariance(alphas,alphas)
k.alphanew_alphanew<-calcCovariance(alphanew,alphanew)
k.alpha_alphanew<-calcCovariance(alphas,alphanew)
m_starnew<-k.alphanew_alpha%*%solve(k.alphaalpha+0.15*diag(1,ncol(k.alphaalpha)))%*%d$eps #obtain prediction mean and covariance
cov_starnew<-k.alphanew_alphanew-k.alphanew_alpha%*% solve(k.alphaalpha+0.15*diag(1,ncol(k.alphaalpha)))%*% k.alpha_alphanew
SampleAlphaPrediction<-rnorm(10000,mean=m_starnew,sd=sqrt(cov_starnew)) #get large sample
hist(SampleAlphaPrediction,xlab="g(5.5)",main="Histogram of g(5.5) distribution",col="blue") #plot histogram (normal distribution as prediction from GP)

#### Interpretation of findings ####
# It's first worth noting some points about the data. It seems to follow the same
# pattern for all points in the domain of alpha except for the point around 5.5.
# This point has an unusually large distance (residual) from the fitted posterior 
# mean compared to the other points. The effect of this can be seen in the confidence
# band which has a small bump above the point, creating slightly greater width than elsewhere 
# in the domain. I expected this point to have a larger effect on the width of the confidence
# interval which is worrying and a reason for suspicion. This point is much closer to the upper
# limit of the confidence interval. A histogram is worth looking at for the random 
# variable (prediction) at this point. The histogram displays how large the variance 
# of the prediction is at g(5.5), with values seen in the range 0 to approximately 3.5.
# This doesn't give us a great deal more information than what can be seen in the first plot.
# Linking back to the first plot, the data value is close to 3, at the upper end of the 
# prediction interval from the histogram. More points should potentially be considered near 
# to 5.5 to get a realistic idea of how the function behaves near to this point as a consequence 
# of the factors discussed, as it doesn't follow the trend of the rest of the data and the 
# GP model has not adapted too much given this point.
