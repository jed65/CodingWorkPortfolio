%% Question 2a)
f=@(x) exp(-x^2); %define the function to integrate
GQEstimate=GaussQuadrature(f); %pass this function and use Gauss Quadrature to integrate
TrueValue=sqrt(pi)*erf(1);
errorInGQ=abs(TrueValue-GQEstimate);
disp(['The estimate from Gauss quadrature with n=10 points is:',num2str(GQEstimate)]);
disp(['The true value of the integral is:',num2str(TrueValue)]);
disp(['The absolute value of the error is:',num2str(errorInGQ)]);

%% Question 2c)
n=[3,5,7]; %vector of n values to pass to function
f=@(x) exp(-x).*cos(pi*x); %function we want to approximate
x=linspace(-1,1,200); %vector of values to evaluate at
q_n=zeros(length(n),length(x)); %initialise matrix of zeroes to store values in approximation
for i=1:length(n)
    c=ComputeLSQCoefficients(f,n(i)); %compute least square coefficients
    for k=0:length(c)-1
        q_n(i,1:length(x))=q_n(i,1:length(x))+c(k+1)*legendreP(k,x); %add on c_k*phi_k for each k to form approximation
    end
    figure(i);
    plot(x,exp(-x).*cos(pi.*x),'b-');
    hold on
    plot(x,q_n(i,1:length(x)));
    hold off
    title(['Least-squares approximation for n=',num2str(n(i))]); 
    legend('f(x)=exp(-x)cos(pi*x)','LSQ approximation','Location','best');
end
figure(4);
plot(x,exp(-x).*cos(pi.*x),'b-');
hold on
for i=1:length(n)
    plot(x,q_n(i,1:length(x)));
end
hold off
title('Least-squares approximation for various n'); 
legend('f(x)=exp(-x)cos(pi*x)','LSQ approximation (n=3)','LSQ approximation (n=5)','LSQ approximation (n=7)','Location','best');
