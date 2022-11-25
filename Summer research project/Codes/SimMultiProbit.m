%Step 1: Form covariance matrices and simulate the 3 processes on interval
X=1477.5:5:3847.5; %Define position of the different variables at centres of the box (1345 to 33325 is over the whole gap)
D=zeros(length(X),length(X)); %Initialise distance matrix
C=zeros(length(X),length(X)); %Initialise covariance matrix
C2=zeros(length(X),length(X)); %Initialise covariance matrix
C3=zeros(length(X),length(X)); %Initialise covariance matrix (X_3)
X=X*0.0353; %Convert the position in pixels into millimetres
phi1=[40.099, 1]; %Define correlation length and signal s.d. (pi_2)
phi2=[10.025, 1];  
phi3=[25, 1];
for i=1:length(X)
    for j=1:length(X)
        D(i,j)=sqrt((X(i)-X(j))^2); %Calculate Euclidean distance
        C(i,j)=(phi1(2)^2)*(1+(sqrt(5)/phi1(1))*D(i,j)+(5/(3*phi1(1)^2))*(D(i,j)^2))*exp(-(sqrt(5)/phi1(1))*D(i,j)); %Calculate covariance (Matern 5/2)
        C2(i,j)=(phi2(2)^2)*(1+(sqrt(5)/phi2(1))*D(i,j)+(5/(3*phi2(1)^2))*(D(i,j)^2))*exp(-(sqrt(5)/phi2(1))*D(i,j));
        C3(i,j)=(phi3(2)^2)*(1+(sqrt(5)/phi3(1))*D(i,j)+(5/(3*phi3(1)^2))*(D(i,j)^2))*exp(-(sqrt(5)/phi3(1))*D(i,j));
    end
end
mu1=(0.3226*1e-10)*ones(1,length(X)); %Define means of the variables (all same since stationary)
mu2=(0.3714*1e-10)*ones(1,length(X)); %Are these means correct?? 
mu3=zeros(1,length(X));
R1=mvnrnd(mu3,C,2000); %Simulate realisations of the first process, all of the below have mean zero (means added later)
R2=mvnrnd(mu3,C2,2000); %Second process 
R3=mvnrnd(mu3,C3,2000); %Third process 

%Step 2: Define smooth 'indicator' function and combine with the above to
%form the desired simulation
sigma=0.1; %This is the deviation of the indicator function, don't want too low so that it is reasonably differentiable
m=-0.6283*sqrt(1+sigma^2); %Define the mean of the indicator function from sigma, this comes from requiring E[Phi(X_3)]=pi_2
link=zeros(2000,length(X)); %Initialise these two matrices
GMix=zeros(2000,length(X));
for i=1:2000
    link(i,:)=normcdf(R3(i,:),m,sigma); %Use the normal cumulative distribution to find the values using the X3 simulation, this defines the values of the approx. indicator function
    GMix(i,:)=(mu1+6.14*1e-12*R1(i,:)).*link(i,:)+(mu2+1.4707*1e-11*R2(i,:)).*(ones(1,length(X))-link(i,:)); %Use form of model to obtain simulations, each row is one simulation
end



%% Step 3: Produce plots of the CDF and the simulation
figure(1);
subplot(2,1,1);
x=-2:0.001:2;
pcdf=normcdf(x,m,0.1);
plot(x,pcdf,'k-');
xlabel('x in [-2,2]');
ylabel(['CDF (sigma=',num2str(sigma),')']);
title(['Plot of inverse probit with mu=',num2str(m)]);

subplot(2,1,2);
plot(X,GMix(1,:),'-');
hold on
for j=2:5
    plot(X,GMix(j,:),'-');
end
hold off
xlabel('Distance across gap (mm)');
ylabel('Permeability (m^2)');
title('5 permeability simulations across subset of gap');

%% Section 
figure(2);
%subplot(2,1,1);
lag=0:1:400; %Define the vector of lags
acf=zeros(2000,length(lag)); %Initialise matrix of sample autocorrelations
variance=zeros(2000,1);
for k=1:2000 %This for loop finds the sample acf for each simulation in term, assigning each simulation to a row of acf
    mu=mean(GMix(k,:));
    [acf1,var]=sampleacf(GMix(k,:),400,mu);
    acf(k,:)=acf1;
    variance(k)=var;
end
ACFavg2=mean(acf); %Find the average of the sample autocorrelations 
VarAvg=mean(variance);
dist=lag*5*0.0353;
hold on
for l=1:10 %Plot first 10 sample autocorrelations
    plot(dist,acf(l,:),'-','LineWidth',0.1);
end
plot(dist,ACFavg2,'k-','LineWidth',2.5); %Make this line thicker than the rest 
hold off
xlabel('Distance (mm)');
ylabel('Autocorrelation');
title('Simulation Autocorrelations');

%% Section
subplot(2,1,2);
q=zeros(1,length(lag)-2);
TrueCov=YCovariance(lag,40.099,0,10.025,0,25,q);
TrueCorr=TrueCov/TrueCov(1); %Divide TrueCov by the variance to get the correlation
TrueCorr(length(TrueCorr))=0; %The last element is hugely negative for some reason, tried for many lags and always does the same thing 
plot(dist,TrueCorr,'b-');
hold on
ICov=TestIndicatorCov(lag,40.099,10.025,25);
ICorr=ICov/ICov(1);
ICorr(length(ICorr))=0;
plot(dist, ICorr,'r-');
hold off
xlabel('Distance (mm)');
ylabel('Autocorrelation');
title('Comparing Autocorrelations');
legend('Model ACF','Exact Indicator ACF','Location','northeast'); 

