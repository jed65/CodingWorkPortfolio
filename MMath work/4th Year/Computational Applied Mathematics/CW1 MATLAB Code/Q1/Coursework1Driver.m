%% Check question 1a) function 
l=[0,3,2];
d=[2,4,3];
u=[0.5,1,0];
f=[12.75,33,14.5];
x=SolveTridiagonalSystem(l,d,u,f);
disp('Solution to the problem is:');
disp(x);

%% Check question 1b) and 1c) functions
x=linspace(0,1/2,3);
f=sin(pi*x);
c=FindNaturalSplineCoefficients(f);
[q3,xhat]=EvaluateSpline(c,x,0,1/2);
x2=linspace(-2,2,100);
f2=sin(pi*x2);
plot(x2,f2,'b-');
hold on
plot(xhat,q3,'r-');
hold off

%% Part 1d)
n=16; %set first value of n
nHolder=zeros(1,4); %initialise vectors to hold n values and errors
errorHolder=zeros(1,4);
errorHolder2=zeros(1,4);
figure(1);
for i=1:4
    x=linspace(-1,1,n+1); %get vector of nodes
    f=exp(-x).*cos(6*pi.*x); %evaluate function at nodes
    c=FindNaturalSplineCoefficients(f); %find the natural spline coefficients
    [q3,xhat]=EvaluateSpline(c,x,-1,1); %evaluate the spline within the interval at many uniform points xhat
    
    subplot(2,2,i); %this is so we get a plot for each
    hold on
    plot(xhat,q3);
    xlabel('x');
    ylabel('q3(x)');
    title(['Spline for n=',num2str(n)]);
    
    fTrue=exp(-xhat).*cos(6*pi.*xhat); %evaluate f at points that spline has been evaluated at 
    errorHolder(i)=max(abs(fTrue-q3)); %find the infinity norm of the error and assign
    errorHolder2(i)=max(abs(fTrue(floor(length(xhat)/4):floor(3*length(xhat)/4))-q3(floor(length(xhat)/4):floor(3*length(xhat)/4)))); %find infinity norm of middle 50% of interval
    nHolder(i)=n;
    
    
    
    n=2*n; %update n
end
hold off
figure(2);
plot(log10(nHolder),log10(errorHolder),'b-*');
hold on
n=16:nHolder(4);
plot(log10(nHolder),log10(errorHolder2),'m-*');
plot(log10(n),2*ones(1,length(n))+2*log10(2./n),'g--');
plot(log10(n),3.1*ones(1,length(n))+4*log10(2./n),'k--');
hold off
xlabel('log10(n)');
ylabel('log10(Error)');
title('Log-log plot of infinity norm of error vs. no. of nodes used in spline');
legend('Error across whole interval','Error across middle 50% of interval','Reference line slope -2','Reference line slope -4','Location','northeast');
hold off
VarNames={'n','ErrorAcrossWholeInterval','ErrorAcrossMiddleOfInterval'};
ErrorTable=table(nHolder',errorHolder',errorHolder2','VariableNames',VarNames);
disp(ErrorTable);
