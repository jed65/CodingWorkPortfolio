#### Exercise 1 ####
library(MASS)
library(mvtnorm)
#everything done here in one function, in interval [0,5] and 1000 points
gaussprocess<-function(from=0,to=1,K=function(s,t){min(s,t)},m=100){
  t<-seq(from=from,to=to,length.out=m)
  Sigma<-sapply(t,function(s1){
    sapply(t,function(s2){
      K(s1,s2)
    })
  })
  path<-mvrnorm(mu=rep(0,times=m),Sigma=Sigma)
  return(data.frame("t"=t,"xt"=path))
}
Kdata<-gaussprocess(from=0,to=1,K=function(s,t){min(s,t)}) #extract dataframe
plot(xt~t,data=Kdata,xlab="t",ylab="xt") #produce plot
lines(Kdata$t,Kdata$xt,col="blue")
#To change the kernel function and implement, remember to change in both function and Kdata!

#### Exercise 2 ####
N<-1000 #this is the number of points over the interval
cjVec<-double(N) #first create the vector of independent cj's
for (j in 1:N){
  lambda<-(j*pi)^(-1)
  cjVec[j]<-rnorm(1,mean=0,sd=lambda)
}
#Now write function that finds the phi_j(x) over the interval
gaussprocess2<-function(from=0,to=1,m=1000,j1,term=function(j,x1){sin(j*pi*x1)}){
  x<-seq(from=from,to=to,length.out=1000)
  fx<-sapply(x,function(x1){
    term(j1,x1)
  })
  return(fx)
}
#Call the function to compute phi_j(x) for lots of j and combine into single vector
xt<-rep(0,1000)
for (j in 1:N){
  fx<-gaussprocess2(from=0,to=1,m=1000,j,term=function(j,x1){sin(j*pi*x1)})
  xt<-xt+cjVec[j]*fx
}
x<-seq(from=0,to=1,length.out=1000)
plot(x,xt,xlab="t",ylab="xt") 