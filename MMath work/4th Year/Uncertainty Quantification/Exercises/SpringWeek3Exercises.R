#### Exercise 2 ####
EvaluateIntegrand<-function(t,lambda,degree){ #function to evaluate integrand
  if (degree==1){
    integrand<-exp(-lambda*t)*(lambda-0.5)*1.25
  }
  else if (degree==2){
    integrand<-exp(-lambda*t)*(lambda^2-lambda+(59/300))*1.25
  }
  return(integrand)
}


CTR<-function(J,t,degree){
  delta<-0.8/J
  jValues<-rep(0.1,J+1)+delta*seq(0,J,by=1)
  Ictr<-0
  for (i in 1:length(jValues)){
    if (jValues[i]==0.1){
      Ictr<-Ictr+0.5*delta*EvaluateIntegrand(t,jValues[i],degree)
    }
    else if (jValues[i]==0.9){
      Ictr<-Ictr+0.5*delta*EvaluateIntegrand(t,jValues[i],degree)
    }
    else {
      Ictr<-Ictr+delta*EvaluateIntegrand(t,jValues[i],degree)
    }
  }
  return(Ictr)
}

t<-5
epsilon<-seq(-1,1,by=0.05)
uExact<-exp(-(0.5+0.4*epsilon)*t)
u0<-rep((5/(4*t))*(exp(-0.1*t)-exp(-0.9*t)),length(uExact))
plot(epsilon,uExact,ylim=c(-0.5,1))
lines(epsilon,uExact,col="dark blue")
lines(epsilon,u0,col="red")
u1<-(75/4)*CTR(100,t,1)
u1DiscreteApprox<-u0+u1*(0.4*epsilon)
lines(epsilon,u1DiscreteApprox,col="black")
u2<-(28125/64)*CTR(100,t,2)
u2DiscreteApprox<-u0+u1*(0.4*epsilon)+u2*((4/25)*epsilon^2-(4/75))
lines(epsilon,u2DiscreteApprox,col="green")

#### Exercise 3- TO DO ####


#### Exercise 4 ####
#Implements Stochastic Galerkin for N=1
u0SG<-(8/(41+5*sqrt(41)))*exp(((-10+sqrt(41))/30)*t)+((41+5*sqrt(41))/82)*exp(((-10-sqrt(41))/30)*t)
u1SG<-(-2*sqrt(41)/41)*exp(((-10+sqrt(41))/30)*t)+(2*sqrt(41)/41)*exp(((-10-sqrt(41))/30)*t)
u1StochasticGalerkin<-u0SG+u1SG*epsilon
lines(epsilon,u1StochasticGalerkin,col="purple",lty="dashed")