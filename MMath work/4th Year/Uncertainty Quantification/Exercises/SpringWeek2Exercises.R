#### Exercise 1 ####
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

eps<-rep(0,9)
for (n in 9)
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

#Or... Use this command in R from the vipor package
#### Exercise 2 ####
library(vipor)
J<-seq(10,10000,by=10)
uQMC<-rep(0,length(J))
for (i in 1:length(J)){
  eps2<-vanDerCorput(J[i],base=2,start=1)
  lambda<-0.1*rep(1,J[i])+0.8*eps2
  u5<-0
  for (j in 1:J[i]){
    u5<-u5+(1/J[i])*exp(-lambda[j]*5)
  }
  uQMC[i]<-u5
}
m_u5<-rep((-0.25)*(exp(-4.5)-exp(-0.5)),length(J))
errorQMC<-(uQMC-m_u5)^2
plot(log(J),log(errorQMC),cex=0.5)
plot(J,errorQMC,log="xy")
OrderlogJoverJ<-8e-5*log(J)/J
lines(J,OrderlogJoverJ,col="red",lty="dashed")

#### Exercise 4 ####
EvaluateU<-function(t,lambda){ #function to evaluate uprime(t,lambda)
  lambdaPrime<-0.5-0.4+2*0.4*lambda
  u<-exp(-lambdaPrime*t)
  return(u)
}

CTR<-function(J,t){
  delta<-1/J
  jValues<-delta*seq(0,J,by=1)
  Ictr<-0
  for (i in 1:length(jValues)){
    if (jValues[i]==0){
      Ictr<-Ictr+0.5*delta*EvaluateU(t,jValues[i])
    }
    else if (jValues[i]==1){
      Ictr<-Ictr+0.5*delta*EvaluateU(t,jValues[i])
    }
    else {
      Ictr<-Ictr+delta*EvaluateU(t,jValues[i])
    }
  }
  return(Ictr)
}

CSR<-function(J,t){
  delta<-1/J
  jValues<-delta*seq(0,J,by=1)
  Icsr<-0
  for (i in 1:length(jValues)){
    if (jValues[i]==0){
      Icsr<-Icsr+(delta/3)*EvaluateU(t,jValues[i])
    }
    else if (i%%2==0){
      Icsr<-Icsr+(delta/3)*4*EvaluateU(t,jValues[i])
    }
    else if (jValues[i]==1){
      Icsr<-Icsr+(delta/3)*EvaluateU(t,jValues[i])
    }
    else {
      Icsr<-Icsr+(delta/3)*2*EvaluateU(t,jValues[i])
    }
  }
  return(Icsr)
}

#Composite Trapezium Rule first
J<-seq(10,500,by=10)
t<-5
m_u5<-rep((-0.25)*(exp(-0.9*t)-exp(-0.1*t)),length(J))
CTRApprox<-rep(0,length(J))
for (i in 1:length(J)){
  CTRApprox[i]<-CTR(J[i],5)
}
errorCTR<-abs(m_u5-CTRApprox)
plot(J,errorCTR,log="xy")
JandError<-data.frame(J,errorCTR)
Order1overJSquare<-0.19873*(1/J)^2
lines(J,Order1overJSquare,col="red",lty="dashed")

#Composite Simpson's rule next
J<-seq(10,500,by=10)
t<-5
m_u5<-rep((-0.25)*(exp(-0.9*t)-exp(-0.1*t)),length(J))
CSRApprox<-rep(0,length(J))
for (i in 1:length(J)){
  CSRApprox[i]<-CSR(J[i],5)
}
errorCSR<-abs(m_u5-CSRApprox)
plot(J,errorCSR,log="xy")
JandError<-data.frame(J,errorCSR)
Order1overJFour<-0.2117*(1/J)^4
lines(J,Order1overJFour,col="red",lty="dashed")

#### Exercise 5 ####
library(pracma) #contains function gaussLegendre which gives weights and nodes for gauss quadrature

EvaluateUForGQ<-function(t,lambda){
  lambdaPrime<-0.5+0.4*lambda
  u<-exp(-lambdaPrime*t)
  return(u)
}

GaussQuadrature<-function(J,t){
  nodesAndWeights<-gaussLegendre(J,a=-1,b=1)
  x<-nodesAndWeights$x
  w<-nodesAndWeights$w
  Igq<-0
  for (i in 1:length(x)){
    Igq<-Igq+w[i]*EvaluateUForGQ(t,x[i])
  }
  return(Igq)
}

J<-c(2,5,10,15,20,25,50,100,500)
t<-5
m_u5<-rep((-0.25)*(exp(-0.9*t)-exp(-0.1*t)),length(J))
GQApprox<-rep(0,length(J))
for (i in 1:length(J)){
  GQApprox[i]<-0.5*GaussQuadrature(J[i],5)
}
errorGQ<-abs(m_u5-GQApprox)
plot(J,errorGQ,log="xy") #Machine precision about reached after using 10 points
.Machine$double.eps #This is machine precision, which is 2.22e-16 on my machine

#### Exercise 6 ####
EvaluateUExercise6<-function(beta,lambda){
  u<-1/(1+exp(-2*beta*(lambda+0.5)))-1/(1+exp(-2*beta*(lambda-0.5)))
  return(u)
}

CTRExercise6<-function(J,beta){ #repeated from before but replaces t
  delta<-1/J
  jValues<-delta*seq(0,J,by=1)
  Ictr<-0
  for (i in 1:length(jValues)){
    if (jValues[i]==0){
      Ictr<-Ictr+0.5*delta*EvaluateUExercise6(beta,jValues[i])
    }
    else if (jValues[i]==1){
      Ictr<-Ictr+0.5*delta*EvaluateUExercise6(beta,jValues[i])
    }
    else {
      Ictr<-Ictr+delta*EvaluateUExercise6(beta,jValues[i])
    }
  }
  return(Ictr)
}

GaussQuadratureExercise6<-function(J,beta){ #again repeated from before with t replaced
  nodesAndWeights<-gaussLegendre(J,a=-1,b=1)
  x<-nodesAndWeights$x
  w<-nodesAndWeights$w
  Igq<-0
  for (i in 1:length(x)){
    Igq<-Igq+w[i]*EvaluateUExercise6(beta,x[i])
  }
  return(Igq)
}

beta<-50 #Can set this to 5 or 50
J<-c(2,5,10,15,20,25,50,100,500)
ExactFromCTR<-rep(CTRExercise6(1000000,beta),length(J))
CTRApprox<-rep(0,length(J))
GQApprox<-rep(0,length(J))
for (i in 1:length(J)){
  CTRApprox[i]<-CTRExercise6(J[i],beta)
  GQApprox[i]<-0.5*GaussQuadratureExercise6(J[i],beta) #for some reason need to multiply by 0.5 again
}
errorCTR<-abs(ExactFromCTR-CTRApprox)
errorGQ<-abs(ExactFromCTR-GQApprox)
plot(J,errorGQ,log="xy",col="red") #points from Gauss quadrature in red
points(J,errorCTR,col="black") #points from trapezium rule in black, when beta=50 result is exact!


