##### Question 1 - Vectorizing a function #####
# Give example of vector theta1/theta2 not working
dat<-read.csv("https://people.bath.ac.uk/kai21/COVID19_MARCH_JUNE.csv",header=TRUE)
a<-dat$age[1:17] # get vector of age midpoints
k<-(dat$population[1:17]+dat$population[18:34])/100000 # sizes
theta1=seq(-5,-2,length.out=34) # theta values to test
theta2=seq(0.05,0.2,length.out=34)
y_male   <- subset(dat$deaths_COVID,dat$sex=="male")
y_female <- subset(dat$deaths_COVID,dat$sex=="female")
y        <- y_male+y_female
#The below is a non-vectorised function
nll<-function(theta1,theta2,y,a,k){
  
  lin_pred     <- theta1+theta2*a
  
  mu           <- k*exp(lin_pred)
  
  sum_log_dens <- sum(dpois(x    = y,
                            lambda = mu,
                            log  = TRUE))
  
  -sum_log_dens
}
nll(theta1,theta2,y,a,k) #Should give 34 values but it instead only gives one!

# Need to fix the negative log-likelihood function
nll<-function(theta1,theta2,y,a,k){
  
  lin_pred     <- theta1+theta2*a
  
  mu           <- k*exp(lin_pred)
  
  sum_log_dens <- sum(dpois(x    = y,
                            lambda = mu,
                            log  = TRUE))
  
  -sum_log_dens
}
nll<-Vectorize(FUN=nll,vectorize.args = c("theta1","theta2")) # Vectorize
nll(theta1,theta2,y,a,k) # This is now a vector with 34 values

##### Question 2 - Plot data/fit linear model #####
dat$rate<-dat$deaths_COVID/(dat$population/100000) # Add new variable
par(mfrow=c(1,1))
plot(dat$age[1:17],log(dat$rate[1:17]),xlab='Age',ylab='log(Death rate)',main='Male (black), female (blue)')
points(dat$age[18:34],log(dat$rate[18:34]),xlab='Age',ylab='log(Death rate)',main='Female',col="blue")
#Fit a linear model
fitted<-lm(log(rate)~age,data=dat[1:17,])
fitted$coefficients

##### Question 3 - Produce contour plot #####
theta1_vals<-seq(-5,-2,length.out=100)
theta2_vals<-seq(0.05,0.2,length.out=100)
nll_grid_theta1_theta2<-outer(X=theta1_vals,
                              Y=theta2_vals,
                              FUN=nll,
                              y=y,
                              a=a,
                              k=k
                              )
levels_nll<-quantile(x=nll_grid_theta1_theta2,
                     probs=seq(0,0.3,length.out=10)) # want to minimise! 
contour(x=theta1_vals,
        y=theta2_vals,
        z=nll_grid_theta1_theta2,
        levels=levels_nll,
        xlab="Theta_1",
        ylab="Theta_2",
        main="Contour plot")
points(x=fitted$coefficients[1],
       y=fitted$coefficients[2],
       cex=1,
       col="red",
       pch=16)