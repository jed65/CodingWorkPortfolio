function x = SolveTridiagonalSystem(l,d,u,f)
% SolveTridiagonalSystem.m: A function to solve the problem Ax=f where A is
% a tridiagonal n x n matrix and x/f are n x 1 vectors using Gaussian
% elimination algorithm. Inputs are the entries on the lower diagonal (plus
% l(1)=0), the entries on the main diagonal, the entries on the upper
% diagonal (plus u(n)=0) and the vector f. Output is the solution x. The
% function displays an instructive error if the input vectors are
% insufficient to solve the problem and also in the situation that the
% elimination step attempts to divide by zero. 


n=length(d); %find the dimension of the problem
x=zeros(1,n); %initialise vector to store solution

if (length(l)~=n)||(length(u)~=n)||(length(f)~=n) %display error if lower or upper diagonal not length n
    disp('Error: Length of vectors for lower/upper diagonals of n x n tridiagonal A or f in Ax=f not correct length');
    disp('Tip: Remember to set l(1)=0 and u(n)=0 so that their length is n');
    
elseif (l(1)~=0)||(u(n)~=0) %display error if first element of l or final element of u not equal to zero
    disp('Error: First element of lower diagonal or final element of upper diagonal not equal to zero');
    
elseif d(1)==0 %check that the first element of the diagonal non-zero, as need to divide by it
    disp('Error: First element of diagonal set equal to zero - cannot divide by zero');
    
else %able to at least begin to solve system (IMPLEMENTATION OF ALGORITHM STARTS HERE)
    errorCheck=0; %this keeps tabs on if there has been an error
    for i=2:n %for loop for implementing elimination step
        if d(i-1)~=0
            d(i)=d(i)-u(i-1)*l(i)/d(i-1); %update diagonal element
            f(i)=f(i)-f(i-1)*l(i)/d(i-1); %update f_i
        else
            disp(['Error: Attempt to divide by zero in elimination step, relevant row is ',num2str(i-1)]);
            disp('Tip: Try using a pivoting method to rearrange rows');
            errorCheck=1;
            break 
        end
    end
    if errorCheck==0 %there has been no errors in elimination step if this holds, continue with backward substitution
        x(n)=f(n)/d(n); %find x_n
        for i=1:n-1
            x(n-i)=(f(n-i)-u(n-i)*x(n-i+1))/d(n-i); %find the other values of x 
        end
    end
    %If there has been an error, these x are not calculated and a zero
    %vector is returned!
end




