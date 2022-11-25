function Jf = ComputeJacobian(F,x,epsilon)
%ComputeJacobian.m: Function to compute the Jacobian of a real
%n-dimensional vector-valued function F, which should be passed to this
%function as an array of anonymous functions. Other inputs are x, the
%vector at which you want the Jacobian evaluated, and epsilon which is a
%small number used to approximate derivatives via finite difference. The
%output is the n x n Jacobi matrix at the vector x.


n=length(x); %obtain n by finding length of input x
Jf=zeros(n,n); %initialise Jacobian matrix
for i=1:n %i fixes F_i (component of function to evaluate)
    for j=1:n %j fixes an x_j to differentiate with
        xPlusEps=zeros(n,1); 
        xPlusEps(j)=epsilon; %set jth element equal to epsilon, makes canonical basis vector at j (multiplied by epsilon)
        xPlusEps=xPlusEps+x; %add x
        Jf(i,j)=(F{i}(xPlusEps)-F{i}(x))/epsilon; %use finite differences to compute Jacobian
    end
end

    

