#First use the exact solution for the problem and the MSE
times<-seq(0,15,by=0.01)
J<-seq(1,1001,by=10)
m_u<-(-1.25)*(1/times)*(exp(-0.9*times)-exp(-0.1*times))
UJHolder<-matrix(0,nrow=length(J),ncol=length(times))
MSE<-rep(0,30*length(J))
for (i in 1:length(J)){
  for (k in 1:30){
  lambda<-runif(J[i],0.1,0.9)
  exactSols<-matrix(0,nrow=J[i],ncol=length(times))
  UJ<-rep(0,length(times))
  for (j in 1:J[i]){
    exactSols[j,]<-exp(-lambda[j]*times)
    UJ<-UJ+(1/J[i])*exactSols[j,]
  }
  UJHolder[i,]<-UJ
  error<-(UJHolder[i,]-m_u)^2
  MSE[k+30*(i-1)]<-error[length(times)]
  }
}
MSEs<-rep(0,length(J))
start<-1
end<-30
for (i in 1:length(J)){
  MSEs[i]<-mean(MSE[start:end])
  start<-start+30
  end<-end+30
}
plot(log10(J),log10(MSEs),cex=0.4)
lines(log10(J),log10(MSEs),lty="dashed")

#Now use the implicit Euler method
implicit=function(N,lambda,u0){ #get implicit Euler approximation
  u<-rep(0,N+1)
  u[1]<-u0
  h<-15/N
  for (i in 1:N){
    u[i+1]<-u[i]/(1+h*lambda)
  }
  return(u)
}
MSE<-rep(0,30*length(J))
for (i in 1:length(J)){
  for (k in 1:30){
    N<-round(15*sqrt(J[i])) #should be order 1/J now (follows the true solution)
    time<-seq(0,15,length.out=N+1)
    UJHolder<-matrix(0,nrow=length(J),ncol=N+1)
    m_uDiscretise<-(-1.25)*(1/time)*(exp(-0.9*time)-exp(-0.1*time))
    lambda<-runif(J[i],0.1,0.9)
    exactSols<-matrix(0,nrow=J[i],ncol=N+1)
    UJ<-rep(0,N+1)
    for (j in 1:J[i]){
      exactSols[j,]<-implicit(N,lambda[j],1)
      UJ<-UJ+(1/J[i])*exactSols[j,]
    }
    UJHolder[i,]<-UJ
    error<-(UJHolder[i,]-m_uDiscretise)^2
    MSE[k+30*(i-1)]<-error[N+1]
  }
}
MSEs<-rep(0,length(J))
start<-1
end<-30
for (i in 1:length(J)){
  MSEs[i]<-mean(MSE[start:end])
  start<-start+30
  end<-end+30
}
lines(log10(J),log10(MSEs),lty="dashed",col="red")