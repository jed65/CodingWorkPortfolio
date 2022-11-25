function [ICov,X1Cov,X2Cov,q] = TestIndicatorCov(X,l1,thetaX1,l2,thetaX2,a1,a2,L1,L2,l3,I)
ICov=zeros(1,length(X));
X1Cov=zeros(1,length(X));
X2Cov=zeros(1,length(X));
X1Cov(1)=1;
X2Cov(1)=1;
X3Cov=zeros(1,length(X));
m=-0.6283;
m1=0.3226*1e-10; %Next few lines define the relevant quantities from the dissertation
m2=0.3714*1e-10;
sigma1=6.14*1e-12;
sigma2=1.4707*1e-11;
ICov(1)=0.7351*(m1^2+sigma1^2)+0.2649*(m2^2+sigma2^2);
%Convert lags into distances
X=5*X;
X=0.0353*X;
q=zeros(length(X)-2,1);
for i=2:length(X)-1
    if I==0
       X1Cov(i)=(1+(sqrt(5)*X(i))/l1+(5*(X(i)^2))/(3*l1^2))*exp(-(sqrt(5)*X(i))/l1); %First 3 lines of the for loop find covariance for this distance X(i)
       X2Cov(i)=(1+(sqrt(5)*X(i))/l2+(5*(X(i)^2))/(3*l2^2))*exp(-(sqrt(5)*X(i))/l2);
    else
        X1Cov(i)=(1+X(i)/thetaX1)^(-1)*[cos(X(i)/l1)+a1*sin(X(i)/L1)];
        X2Cov(i)=(1+X(i)/thetaX2)^(-1)*[cos(X(i)/l2)+a2*sin(X(i)/L2)];
    end
    X3Cov(i)=(1+(sqrt(5)*X(i))/l3+(5*(X(i)^2))/(3*l3^2))*exp(-(sqrt(5)*X(i))/l3);
    Sigma=[1 X3Cov(i);X3Cov(i) 1]; %Define covariance matrix of bivariate normal for distance X(i)
    p=mvncdf([m m],[0 0],Sigma);
    q(i-1)=1-0.5298+p;
    ICov(i)=(m1^2+sigma1^2*X1Cov(i))*q(i-1)+m1*m2*(1.4702-2*q(i-1))+(m2^2+sigma2^2*X2Cov(i))*(1-1.4702+q(i-1)); %Combine all of the above to obtain Y covariance
end
EY0=0.7351*m1+0.2649*m2; %Define the mean of Y
ICov=ICov-EY0^2*ones(1,length(X)); %Subtract the squared mean to get the covariance function
