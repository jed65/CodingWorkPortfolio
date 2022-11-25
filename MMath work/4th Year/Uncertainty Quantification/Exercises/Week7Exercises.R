#### Ex.1 ####

J<-seq(10,210,by=25)
Sim_squareError<-rep(0,length(J))
Theoretical_squareError<-rep(0,length(J))
times<-seq(0,15,length.out=1000)
m_u<- (-1)*1.25*(1/times)*(exp(-0.9*times)-exp(-0.1*times))
m_u[0]<-0
uMC<-rep(0,length(times))
uLambdaJ<-rep(0,length(times))
m_u_at_5<-0.25*(exp(-0.5)-exp(-4.5))

for (k in 1:length(J)){
expectation_u<-0
uLambdas<-matrix(0,nrow=J[k],ncol=length(times))
for (v in 1:10){
  lambda<-runif(J[k],min=0.1,max=0.9)
  for (i in 1:J[k]){
    for (j in 1:length(times)){
      uLambdaJ[j]<-exp(-1*lambda[i]*times[j])
    }
    uMC<-uMC+uLambdaJ
    uLambdas[i,]=uLambdaJ
  }
  uMC<-uMC*(1/J[k])
  expectation_u<-expectation_u+(uMC[334]-m_u_at_5)^2
}
varMC<-rep(0,length(times))
for (i in 1:J[k]){
  varMC<-varMC+(uLambdas[i,]-uMC)^2
}
varMC<-varMC*(1/(J[k]-1))
Theoretical_squareError[k]<-varMC[334]/J[k]

expectation_u<-expectation_u/20
Sim_squareError[k]<-expectation_u
}
#Takes a minute for the above to run



plot(times,uMC,xlab="t",ylab="u",lty="solid",cex=0.01)
upper<-uMC+1.96*sqrt(varMC/J)
lower<-uMC-1.96*sqrt(varMC/J)
lines(times,upper,col="dark grey",lty=2)
lines(times,lower,col="dark grey",lty=2) #looks similar to the plot in the notes
lines(times,m_u,col="red") #this is the analytical solution
hist(uLambdas[,334],xlab="u(5,lambda)",main="Histogram of u(5,lambda) values",xlim=c(exp(-4.5),exp(-0.5)))
u<-seq(exp(-4.5),exp(-0.5),length.out=100)
u_5_density<-1/(2*0.4*5*u)
lines(u,u_5_density,col="blue",lty="solid") #looks like they match up quite well
par(mfrow=c(1,2))
plot(log(J),log(Sim_squareError),col="blue",main="Log-log for MC Simulation")
plot(log(J),log(Theoretical_squareError),col="red",main="Log-log for Theoretical")

#### Ex.2 ####
B<-12
lambda<-runif(B,min=0.1,max=0.9)
lambdaRange<-seq(0.1,0.9,length.out=B)
Q_5<-rep(0,B)
for (j in 1:B){
 Q_5[j]<-exp(-1*lambda[j]*5)
}
par(mfrow=c(1,1))
plot(lambda,Q_5,cex=0.5,xlab="Lambda",ylab="u(5,Lambda)")
lines(lambda[order(lambda)],Q_5[order(lambda)],col="black")

#Now fit a GP
calcSigma <- function(X1,X2) { # Function takes in a discretised grid
  #input X1 x X2 and computes covariance matrix using covariance function k
  v1=1
  v2=0.01
  alpha=2
  r=2
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1)) #fill with zeroes, give rows equal to number of elements of x1
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- v1*exp(-(abs(X1[i]-X2[j])/r)^alpha)+v2*(i==j) #use formula given, remember boolean are 1 or 0 in R
    }
  }
  return(Sigma) 
}
lambdaStar<-seq(0.1,0.9,length.out=50)
k.lstar_l<-calcSigma(lambdaStar,lambda)
k.ll<-calcSigma(lambda,lambda)
k.lstar_lstar<-calcSigma(lambdaStar,lambdaStar)
k.l_lstar<-calcSigma(lambda,lambdaStar)
m_star<-k.lstar_l%*%solve(k.ll+0.01*diag(1,ncol(k.ll)))%*%Q_5
lines(lambdaStar,m_star,col="blue")

