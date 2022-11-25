x=normrnd(0.5,1,[1,500]);
mu=linspace(-1,1,100);
sum=zeros(1,length(mu));
logL=zeros(1,length(mu));
for i=1:length(mu)
    for k=1:length(x)
        sum(i)=sum(i)+(x(k)-mu(i))^2;
    end
    logL(i)=-10*log(2*pi)-0.5*sum(i);
end
plot(mu,logL,'-b');
xlabel('Mu');
ylabel('logL(mu,1)');
mle=mu(logL==max(logL));
n=length(x);
analyticMle=mean(x);
disp(['When n=',num2str(n),' the maximum of the log-likelihood is ',num2str(mle), ' and the analytical mle (mean) is ',num2str(analyticMle)]);

%% Exercise 3 Part 1
n=linspace(50,10000,996);
error=zeros(1,length(n));
for i=50:10:10000
    p=zeros(1,i);
    for j=1:i
        p(j)=binornd(1,0.6); %Theta specified here
    end
    thetam=mean(p);
    error(i/10-4)=abs(0.6-thetam);
end
plot(n,error,'-b');
xlabel('n');
ylabel('Error');
%% Exercise 3 part 2
nstar=500;
N=100000;
counter=0;
theta=0.6;
for i=1:N
   rng(i);
   p=binornd(1,theta,[1,nstar]);
   [~,pci]=mle(p,'distribution','Bernoulli');
   if theta>pci(1) && theta<pci(2)
        counter=counter+1;
    end
end
disp(counter);
proportion=counter/N;
disp(['The proportion of confidence intervals containing the true value of theta=0.6 is ',num2str(proportion)]);
%% Exercise 4
theta=unifrnd(0,1);
disp(theta);
simN=zeros(1,1000);
for i=1:1000
    simN(i)=binornd(1,theta);
end
k=1000*mean(simN);
beta1=1+k;
beta2=1000-k;
thetaArg=0:0.01:1;
pifunc=betapdf(thetaArg,beta1,beta2);
plot(thetaArg,pifunc,'b');
xlabel('Theta');
ylabel('Pi(Theta given D)');
p=[0.025,0.975];
x=betainv(p,beta1,beta2);
width=diff(x);
disp(['The interval with theta=',num2str(theta),' is']);
disp(x);
disp(['The interval has width ',num2str(width)]);
N=1000;
samp=betarnd(beta1,beta2,1,N);
samp=sort(samp);
a=samp(N*0.025);
b=samp(N*0.975);
MCinterval=[a b];
disp('Interval from Monte Carlo is');
disp(MCinterval);

%% Exercise 5
alpha=5;
beta=2;
tau1=gamrnd(alpha,beta);
disp(tau1);
tau=0:0.01:20;
sigma=sqrt(1/tau1);
mu=3;
for n=20:20:80
    x=normrnd(mu,sigma,[1,n]);
    sum=0;
    for i=1:n
        sum=sum+(x(i)-mu)^2;
    end
    pifunc=(tau./(2*pi)).^(n/2).*exp(-0.5.*sum.*tau).*gampdf(tau,alpha,beta);
    marg=[(beta^alpha)*factorial(alpha+n/2-1)]/[factorial(alpha-1)*((2*pi)^(n/2))*((beta+0.5*sum)^(n+1))];
    pifunc=pifunc/marg;
    pifunc=pifunc/max(pifunc);
    plot(tau,pifunc);
    hold on
end
hold off
xlabel('Tau=1/sigma^2');
ylabel('Posterior pdf');
legend('n=20','n=40','n=60','n=80');

%% Method of composition
mu=1;
alpha=2;
beta=3;
rng('default');
tau=gamrnd(alpha,beta);
for n=20:20:100
    sum=0;
    X=normrnd(mu,1/tau,[1,n]);
    for i=1:n
        sum=sum+(X(i)-mu)^2;
    end
    alphaNew=alpha+n/2;
    betaNew=beta+0.5*sum;
    piFunc=gampdf(0:0.5:200,alphaNew,betaNew);
    plot(0:0.5:200,piFunc);
    hold on
end
legend('n=20','n=40','n=60','n=80','n=100');
xlabel('tau');
ylabel('Posterior density');

    



    