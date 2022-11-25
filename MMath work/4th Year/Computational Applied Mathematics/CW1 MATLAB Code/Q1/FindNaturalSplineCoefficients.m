function c = FindNaturalSplineCoefficients(f)
% FindNaturalSplineCoefficients.m: Function taking function values f at n+1
% equally spaced points, returning the coefficients of the natural cubic
% B-spline interpolating f. The natural spline imposes additional two
% conditions that the second derivative of the spline is equal to zero at
% the endpoints of the interval.

c=zeros(1,length(f)+2); %initialise the vector of coefficients
n=length(f)-1;

if length(f)<=2 %output an error if not enough function values given, 
    disp('Error: Not enough function values given to evaluate cubic spline coefficients');
else 
    c(2)=(1/6)*f(1); %set the value of c_0
    c(n+2)=(1/6)*f(n+1); %set the value of c_n
    
    %Now set up the tridiagonal system to solve
    d=4*ones(1,n-1); %there are n-1 equations for n-1 unknowns
    u=ones(1,n-1); 
    u(n-1)=0; %upper diagonal set up with zero for final element
    l=ones(1,n-1); 
    l(1)=0; %lower diagonal set up with zero for first element
    g=f(2:n); %set up the right-hand side with f(x_i)
    g(1)=g(1)-(1/6)*f(1); %alter first and last elements 
    g(n-1)=g(n-1)-(1/6)*f(n+1);
    
    %Solve system
    c(3:n+1)=SolveTridiagonalSystem(l,d,u,g); 
    
    c(1)=2*c(2)-c(3); %set value of c_(-1)
    c(n+3)=2*c(n+2)-c(n+1); %set value of c_(n+1)
    
end
    

