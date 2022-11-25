%% Exercise 2
x=linspace(0,1,20);
y=linspace(0,1,20);
[xgrid,ygrid]=meshgrid(x,y);
p=zeros(length(x),length(y));
for n=1:1000
    for m=1:1000
        f_nm=(2/(n*m*pi^2))*(cos(0.4*n*pi)-cos(0.6*n*pi))*(cos(0.4*m*pi)-cos(0.6*m*pi));
        beta_nmInverse=1/(pi^2*(n^2+m^2));
        for i=1:length(x)
            for j=1:length(y)
                p(i,j)=p(i,j)+f_nm*beta_nmInverse*(2*sin(n*pi.*x(i))*sin(m*pi.*y(j)));
            end
        end
    end
end
surf(x,y,p);
title('p(x,y)');

%% Exercise 4
L=1;
x=linspace(0,L,100);
I=zeros(1,length(x));
for i=1:length(x)
    if x(i)<0.5
        I(i)=x(i);
    else 
        I(i)=L-x(i);
    end
end
kappa=10;
t=1e-3:7.5e-3:3.1e-2;
T=zeros(length(t),length(x));
for j=1:length(t)
for n=1:1000
    for i=1:length(x)
        T(j,i)=T(j,i)+(4*L/(n^2*pi^2))*sin(n*pi/2)*exp(-((n^2*pi^2*kappa)/(L^2))*t(j))*sin((n*pi*x(i))/L);
    end
end
end
plot(x,I,'k-');
hold on
for j=1:length(t)
    plot(x,T(j,1:length(x)));
end
hold off
title('Solution of heat equation with kappa=10');
legend('T(x,0)','T(x,1e-3)','T(x,8.5e-3)','T(x,0.016)','T(x,0.0235)','T(x,0.031)','Location','best');

%% Exercise 5 - Sampling distribution
L=1;
x=linspace(0,L,100);
t=0.03;
kappa=random('Lognormal',log(3),0.25,[500,1]);
T=zeros(length(kappa),length(x));
for j=1:length(kappa)
for n=1:1000
    for i=1:length(x)
        T(j,i)=T(j,i)+(4*L/(n^2*pi^2))*sin(n*pi/2)*exp(-((n^2*pi^2*kappa(j))/(L^2))*t)*sin((n*pi*x(i))/L);
    end
end
end
figure(1);
plot(x,T(1,1:length(x)));
hold on
for j=2:length(kappa)
    plot(x,T(j,1:length(x)));
end
hold off
title('T_k(x,0.03) with random kappa');

mT=zeros(1,length(x));
varT=zeros(1,length(x));
for j=1:length(kappa)
    mT=mT+(1/length(kappa)).*T(j,1:length(x));
end
for j=1:length(kappa)
    varT=varT+(1/(length(kappa)-1)).*(T(j,1:length(x))-mT).^2;
end
figure(2);
plot(x,mT,'r-');
hold on
plot(x,mT-2.*sqrt(varT),'k--');
plot(x,mT+2.*sqrt(varT),'k--');
hold off
title('Uncertainty in T_k(x,0.03)');

figure(3);
histogram(T(1:length(kappa),51),50);
title('Histogram for T_k(0.5051,0.03)');
    
%% Exercise 7
x=linspace(0,1,20);
y=linspace(0,1,20);
[X,Y]=meshgrid(x,y);
CovarianceMatrix=zeros(400,400);
smoothness=0.5;
tau=0.25;
z2=zeros(20,20);
for m=1:100
    rowCounter=1;
for i=1:20
    xFixed=x(i);
    for j=1:20
        yFixed=y(j);
        columnCounter=1;
        for p=1:20
            xPrimeFixed=x(p);
            for q=1:20
                yPrimeFixed=y(q);
                xVector=[xFixed,yFixed];
                xPrimeVector=[xPrimeFixed,yPrimeFixed];
                distance=norm(xVector-xPrimeVector);
                if distance~=0
                    CovarianceMatrix(rowCounter,columnCounter)=(1/(2^(smoothness-1)*gamma(smoothness)))*((distance/tau)^smoothness)*besselk(smoothness,distance/tau);
                else
                    CovarianceMatrix(rowCounter,columnCounter)=1;
                end
                columnCounter=columnCounter+1;
            end
        end
        rowCounter=rowCounter+1;
    end
end
mu=ones(1,400);
z=mvnrnd(mu,CovarianceMatrix,1);
for i=1:20
    z2(i,1:20)=z2(i,1:20)+(1/100).*z((20*(i-1)+1):20*i);
end
end
figure(1);
surf(x,y,z2);
title('E[f(x,y)] from Monte Carlo');
x=linspace(0,1,20);
y=linspace(0,1,20);
[xgrid,ygrid]=meshgrid(x,y);
p=zeros(length(x),length(y));
mean_p_MC=p;
for k=1:100
    evalPoints=rand(2,1);
for n=1:100
    for m=1:100
        beta_nmInverse=1/(pi^2*(n^2+m^2));
        for i=1:length(x)
            for j=1:length(y)
                p(i,j)=p(i,j)+beta_nmInverse*(2*sin(n*pi.*x(i))*sin(m*pi.*y(j)))*(2*sin(n*pi.*evalPoints(1))*sin(m*pi.*evalPoints(2)));
            end
        end
    end
end
mean_p_MC=mean_p_MC+(1/100).*p;
end
figure(2);
surf(x,y,mean_p_MC);
title('E[p(x,y)] from Monte Carlo');

%% Need to do variance
x=linspace(0,1,20);
y=linspace(0,1,20);
G1=zeros(length(x),length(y));
smoothness=0.5;
tau=0.25;
var_p_MC=G1;
G2=G1;
for k=1:100
    evalPoints=rand(4,1);
    xhat=evalPoints(1:2);
    xtilda=evalPoints(3:4);
    distance=norm(xhat-xtilda);
    if distance~=0
        covariance=(1/(2^(smoothness-1)*gamma(smoothness)))*((distance/tau)^smoothness)*besselk(smoothness,distance/tau);
    else
        covariance=1;
    end
for n=1:100
    for m=1:100
        beta_nmInverse=1/(pi^2*(n^2+m^2));
        for i=1:length(x)
            for j=1:length(y)
                G1(i,j)=G1(i,j)+beta_nmInverse*(2*sin(n*pi.*x(i))*sin(m*pi.*y(j)))*(2*sin(n*pi.*evalPoints(1))*sin(m*pi.*evalPoints(2)));
                G2(i,j)=G2(i,j)+beta_nmInverse*(2*sin(n*pi.*x(i))*sin(m*pi.*y(j)))*(2*sin(n*pi.*evalPoints(3))*sin(m*pi.*evalPoints(4)));
            end
        end
    end
end
toAdd=covariance.*G1.*G2;
var_p_MC=var_p_MC+(1/100).*toAdd;
end
figure(3);
surf(x,y,var_p_MC);
title('Var[p(x,y)] from Monte Carlo');

%% Histogram
x=linspace(0,1,19);
y=linspace(0,1,19);
smoothness=0.5;
tau=0.25;
pAtMiddle=zeros(1,200);
for sample=1:200
    z2=zeros(19,19);
    rowCounter=1;
    CovarianceMatrix=zeros(361,361);
for i=1:19
    xFixed=x(i);
    for j=1:19
        yFixed=y(j);
        columnCounter=1;
        for p=1:19
            xPrimeFixed=x(p);
            for q=1:19
                yPrimeFixed=y(q);
                xVector=[xFixed,yFixed];
                xPrimeVector=[xPrimeFixed,yPrimeFixed];
                distance=norm(xVector-xPrimeVector);
                if distance~=0
                    CovarianceMatrix(rowCounter,columnCounter)=(1/(2^(smoothness-1)*gamma(smoothness)))*((distance/tau)^smoothness)*besselk(smoothness,distance/tau);
                else
                    CovarianceMatrix(rowCounter,columnCounter)=1;
                end
                columnCounter=columnCounter+1;
            end
        end
        rowCounter=rowCounter+1;
    end
end
mu=ones(1,361);
z=mvnrnd(mu,CovarianceMatrix,1);
for i=1:19
    z2(i,1:19)=z2(i,1:19)+(1/100).*z((19*(i-1)+1):19*i);
end
randomX=randperm(361);
for MCCounter=1:50
    for n=1:50
    for m=1:50
        beta_nmInverse=1/(pi^2*(n^2+m^2));
        evalPoints=zeros(1,2);
        evalPoints(1)=floor(randomX(MCCounter)/19)*(1/18);
        evalPoints(2)=(1/18)*(randomX(MCCounter)-19*floor(randomX(MCCounter)/19));
        pAtMiddle(sample)=pAtMiddle(sample)+(1/100)*z(randomX(MCCounter))*beta_nmInverse*(2*sin(n*pi*0.5)*sin(m*pi*0.5))*(2*sin(n*pi.*evalPoints(1))*sin(m*pi.*evalPoints(2)));
    end
    end
end
end
histogram(pAtMiddle);
   

%% Exercise 20 - Simulating Gaussian Random Field 
J=300;
h=1;
L=300;
x=0:h:300;
smoothness=1.5;
tau=0.5;
midpoints=zeros(1,length(x)-1);
for i=1:length(x)-1
    midpoints(i)=0.5*(x(i)+x(i+1));
end

C=zeros(length(midpoints),length(midpoints));
for i=1:length(midpoints)
    for j=1:length(midpoints)
        if midpoints(i)~=midpoints(j)
            C(i,j)=0.5*(1/(2^(smoothness-1)*gamma(smoothness)))*(abs((midpoints(i)-midpoints(j))/tau)^smoothness)*besselk(smoothness,abs((midpoints(i)-midpoints(j))/tau));
        else
            C(i,j)=0.5;
        end
    end
end

[V,D]=eig(C);
normalCoefficients=normrnd(0,1,1,J);
GRFSample=zeros(1,length(midpoints));
for i=1:length(midpoints)
    for m=1:J
        GRFSample(i)=GRFSample(i)+sqrt(D(m,m))*sqrt(h)*V(i,m)*normalCoefficients(m);
    end
end
plot(midpoints(1:50),GRFSample(1:50));
xlabel('x');
ylabel('GRF');
title('Sample of Gaussian Random Field with smoothness 1.5');
        
    

    