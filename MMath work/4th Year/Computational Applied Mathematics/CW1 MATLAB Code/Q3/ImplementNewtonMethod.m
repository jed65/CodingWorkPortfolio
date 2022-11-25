function [x,xAtIteration] = ImplementNewtonMethod(F,x0,tol,epsilon,iterationOfInterest)
%ImplementNewtonMethod.m: Function which takes inputs of the n-dimensional
%function F, initial guess x0, tolerance for convergence tol, a value
%epsilon (which is used within the finite difference approximation of the
%Jacobian matrix, see ComputeJacobian.m) and optional input
%iterationOfInterest. The outputs include the final iteration x and (if
%reached/requested) the specific iteration requested by the user. Messages
%are displayed to the screen in the event that convergence is not achieved, 
%the requested iteration is not reached or the method could not be 
%implemented (due to zero determinant in the Jacobian). The method used is
%the Newton method. 


convergenceCheck=1; %initialise the difference between approximations for while loop
n=length(x0); %get dimension of problem
xkPlus1=x0; %initialise value
k=0; %number to keep track of which iteration we're on

%disp(['Iteration:',num2str(k)]); %IGNORE (displays iterations, used as a check)
%disp('Estimate:');
%disp(xkPlus1);

errorNotice=0;
if nargin==4 %assumes missing input is iterationOfInterest, if this is the case then message displayed
    disp('No specific iteration has been requested by user, make sure you replace second output by ~');
    xAtIteration=x0;
end
while convergenceCheck>tol %ensures we carry on until convergence
    xk=xkPlus1; %update the value of x
    Jacobian_k=ComputeJacobian(F,xk,epsilon); %call function from a) to compute Jacobian at x_k
    minusF_k=zeros(n,1);
    for i=1:n
        minusF_k(i)=-F{i}(xk);  %calculate -F(x_k)
    end
    if det(Jacobian_k)==0 %display error if Jacobian not invertible
        disp(['Error! The Jacobian at iteration ',num2str(k),' has zero determinant.']);
        errorNotice=1;
        break
    else
        delta_x_k=Jacobian_k\minusF_k; %solve the system for this k
        convergenceCheck=max(abs(delta_x_k)); %find the infinity norm of the above vector to check for convergence
        xkPlus1=xk+delta_x_k; %find next estimate for x*
        k=k+1;
        
        %disp(['Iteration:',num2str(k)]); %IGNORE (displays iterations, used as a check)
        %disp('Estimate:');
        %disp(xkPlus1);
    end
    if k>10000 %this condition ensures code won't run forever
        disp('Convergence not achieved even after 10,000 iterations');
        errorNotice=1; 
        break
    end
    if (nargin==5) && (k==iterationOfInterest) %display message if iteration of interest requested by user reached
        xAtIteration=xkPlus1;
        disp('Iteration of interest reached');
    end
end
x=xkPlus1; %finally, after convergence set the value of the output
if errorNotice~=1
    disp(['Convergence achieved after ',num2str(k),' iterations.']); %display message notifying user that convergence was achieved
end
if (nargin==5) && (k<iterationOfInterest) %display message if iteration of interest requested by user not reached
    disp('Iteration of interest not reached before convergence, ignore second output');
    xAtIteration=x; 
end
