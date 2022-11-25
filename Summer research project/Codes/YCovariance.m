function [TrueCov,X1Cov,X2Cov,Q] = YCovariance(X,l1,thetaX1,l2,thetaX2,l3,q,sigma)
%YCovariance.m: Calculate the true covariance of the model Y, given a 
%vector of lags and correlation lengths. 
TrueCov=zeros(1,length(X));
X1Cov=zeros(1,length(X));
X2Cov=zeros(1,length(X));
X3Cov=zeros(1,length(X));
X1Cov(1)=1;
X2Cov(1)=1;
m1=0.3226*1e-10; %Next few lines define the relevant quantities from the dissertation
m2=0.3714*1e-10;
sigma1=6.14*1e-12;
sigma2=1.4707*1e-11;
m=-0.6283*sqrt(1+sigma^2);
Sig=[1+sigma^2 1;1 1+sigma^2];
Efsquare=mvncdf([-m,-m],[0,0],Sig);
TrueCov(1)=(m1^2+sigma1^2)*Efsquare+m1*m2*(1.4702-2*Efsquare)+(m2^2+sigma2^2)*(1-1.4702+Efsquare); %This is E[Y(0)^2]
Bigfun=@(x,y,avg,std,Sigma) normcdf(x,avg,std).*normcdf(y,avg,std).*mvnpdf([x,y],[0,0],Sigma); %This is the integrand for the below, sigma is used here implicitly (mean of cdf)
fun=@(x,y,Sigma) Bigfun(x,y,m,sigma,Sigma); 
%Convert lags into distances
X=5*X; %Convert into pixels
X=0.0353*X; %Convert pixels to distances (KEEP?)
%Compute true covariance
%v=3; %Define smoothing parameter (if Matern used)
if q(1)==0 %First time code ran if this satisfied
    I=0;
else
    I=1;
end
for i=2:length(X)-1
    if thetaX1~=0 && thetaX2~=0 %If these are set to zero then a Matern 5/2 function is used for X1/X2
       X1Cov(i)=((1+(X(i)/thetaX1))^(-10))*cos(X(i)/l1);%First 3 lines of the for loop find covariance for this distance X(i)
       X2Cov(i)=((1+(X(i))/thetaX2)^(-2))*cos(X(i)/l2); %The first and second processes have covariance functions with oscillating terms to try and get a negative covariance
    else
       X1Cov(i)=(1+(sqrt(5)*X(i))/l1+(5*(X(i)^2))/(3*l1^2))*exp(-(sqrt(5)*X(i))/l1); %Matern 5/2 covariance functions
       X2Cov(i)=(1+(sqrt(5)*X(i))/l2+(5*(X(i)^2))/(3*l2^2))*exp(-(sqrt(5)*X(i))/l2);
    end
    %argX1=(sqrt(2*v)/l1)*X(i);
    %argX2=(sqrt(2*v)/l2)*X(i);
    %X1Cov(i)=((2^(1-v))/factorial(v-1))*(argX1^v)*besselk(v,argX1);
    %X2Cov(i)=((2^(1-v))/factorial(v-1))*(argX2^v)*besselk(v,argX2);
    X3Cov(i)=(1+(sqrt(5)*X(i))/l3+(5*(X(i)^2))/(3*l3^2))*exp(-(sqrt(5)*X(i))/l3); %Keep the third process as Matern 5/2
    if I==0 %This is the first time code ran, when I need to run integral2
        if i==2
            fun2=@(x,y,avg,std) normcdf(x,avg,std).*normcdf(y,avg,std).*(17.46259612).*exp((-6019.31901).*x.^2+(12038.13801).*x.*y-(6019.31901).*y.^2);
            funTwo2=@(x,y) fun2(x,y,m,sigma);
            q(i-1)=integral2(@(x,y) arrayfun(funTwo2,x,y),-Inf,Inf,-Inf,Inf);
        elseif i==3
            fun3=@(x,y,avg,std) normcdf(x,avg,std).*normcdf(y,avg,std).*(8.732374228).*exp((-1505.200728).*x.^2+(3009.901414).*x.*y-(1505.200728).*y.^2);
            funThree3=@(x,y) fun3(x,y,m,sigma);
            q(i-1)=integral2(@(x,y) arrayfun(funThree3,x,y),-Inf,Inf,-Inf,Inf);
        else
        Sigma=[1 X3Cov(i);X3Cov(i) 1]; %Define covariance matrix of bivariate normal for distance X(i)
        funC=@(x,y) fun(x,y,Sigma); %Input this covariance matrix into the function to obtain the integrand
        q(i-1)=integral2(@(x,y) arrayfun(funC,x,y),-Inf,Inf,-Inf,Inf); %Calculate the double integral E[f(X_3(0))f(X_3(x))]
        end
    else
    end
    TrueCov(i)=(m1^2+(sigma1^2)*X1Cov(i))*q(i-1)+m1*m2*(1.4702-2*q(i-1))+(m2^2+(sigma2^2)*X2Cov(i))*(1-1.4702+q(i-1)); %Combine all of the above to obtain Y covariance
end
EY0=0.7351*m1+0.2649*m2; %Define the mean of Y
TrueCov=TrueCov-EY0^2*ones(1,length(X)); %Subtract the squared mean to get the covariance function
Q=q; %The first time this code runs I want q as an output so that I don't have to call integral2 again which takes ages