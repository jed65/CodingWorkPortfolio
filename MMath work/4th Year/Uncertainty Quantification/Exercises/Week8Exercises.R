#### Ex.3 ####
J<-1000
LambdaOdd<-runif(J,0.1,0.9)
LambdaEven<-rep(1,J)-LambdaOdd
t<-seq(0,15,length.out=150)
#Solve for each value
OddSolutions<-matrix(0,nrow=J,ncol=length(t)) #each row is an odd lambda value
EvenSolutions<-matrix(0,nrow=J,ncol=length(t)) #each row is an even lambda value
for (i in 1:J){
  OddSolutions[i,]<-exp(-1*LambdaOdd[i]*t) #fill with solution
  EvenSolutions[i,]<-exp(-1*LambdaEven[i]*t)
}
m_u<-(-1/(2*0.4*t))*(exp(-0.9*t)-exp(-0.1*t))
covariance<-rep(0,length(t))
for (j in 1:J){
  covariance<-covariance+(OddSolutions[j,]-m_u)*(EvenSolutions[j,]-m_u)
}
covariance<-(1/J)*covariance
covariance #can clearly see that all the elements here are negative
W_J<-rep(0,length(t))
for (k in 1:J){
  W_J<-W_J+OddSolutions[k,]
  W_J<-W_J+EvenSolutions[k,]
}
W_J<-(1/(2*J))*W_J
ReducedVariance<-rep(0,length(t))
ReducedVariance<-(W_J-m_u)^2
sigmasquare<-(-1/(4*0.4*t))*(exp(-1.8*t)-exp(-0.2*t))-(1/(4*0.16*t^2))*(exp(-0.9*t)-exp(-0.1*t))^2
sigmasquare<-(1/(2*J))*sigmasquare
sigmasquare-ReducedVariance #differences are order 10^(-6), very small amounts, need at least J=1000

#### Ex.4 ####
#Need to form the Van der Corput sequence
binary=function(n,K,d){ #Get terms in binary expansion
d[K+1]<-1 #in all cases 1 is the final term
n<-n-2^(K)
for (i in 0:K-1)
{
  if (2^(i)<=n){
    nextK<-i
  }
}
if (n==0){
  return(d)
}
else {
binary(n,nextK,d)
}
}
#Now form Van der Corput sequence
error<-rep(0,15)
for (N in 1:15){
  
eps<-rep(0,N)
for (n in 1:N)
{
  for (i in 0:3)
  {
    if (n>=2^(i)){
       K<-i
    }
  }
  d<-rep(0,K+1)
  d<-binary(n,K,d)
  term<-0
  for (j in 0:K){
    term<-term+d[j+1]*2^(-j-1)
  }
  eps[n]<-term
}
lambdaValues<-0.1*rep(1,N)+0.8*eps
times<-seq(0,15,by=0.1)
uSamples<-matrix(0,nrow=N,ncol=length(times))
uQMC<-rep(0,length(times))
for (i in 1:N){
  uSamples[i,]<-exp(-lambdaValues[i]*times)
  uQMC<-uQMC+(1/N)*uSamples[i,]
}
muTrue<- -1.25*(1/times)*(exp(-0.9*times)-exp(-0.1*times))
muTrue[1]<-1
error[N]<-max(abs((muTrue-uQMC))) #this shows that the method works effectively a bit better
}

plot(times,uQMC,cex=0.2)
lines(times,uQMC,col="black")
lines(times,muTrue,lty="solid",col="blue")
NValues<-seq(1,15,by=1)
plot(NValues,error)

