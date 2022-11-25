function [q3,xhat]=EvaluateSpline(c,x,a,b)
% EvaluateSpline.m: Function taking nodes and coefficients of the
% interpolating spline as well as the endpoints of the interval [a,b], that
% evaluates the interpolating spline at many more uniform points within the
% interval of interest, returning both the values of the spline at the
% evaluation points and the evaluation points themselves (for plotting).

n=length(x)-1; 
h=x(2)-x(1); %the difference between successive points should be the same everywhere in uniform grid
xhat=linspace(a,b,20*n+1); %form vector of uniform evaluation points in [a,b]
q3=zeros(1,length(xhat)); %initialise vector to hold the evaluation

errorCheck=0; %variable that identifies whether the inputs are correct
if (x(1)~=a) || (x(n+1)~=b)
    disp('Error: First node not equal to lower interval limit or final node not equal to upper interval limit');
    errorCheck=1;
end
for i=1:length(x)-1
    stepLength=x(i+1)-x(i);
    if (stepLength<h-0.0001)||(stepLength>h+0.0001) %pass error if grid of nodes not uniform
        disp('Error: Nodes x passed not a uniform grid');
        errorCheck=1;
        break
    end
end
if length(c)~=length(x)+2 %pass error if not enough coefficients or nodes given
    disp('Error: Length of vector of coefficients or nodes insufficient');
    errorCheck=1;
elseif a>=b %pass error if the limits of the interval given not sensible
    disp('Error: Lower limit of interval set larger than upper limit');
    errorCheck=1;
end

if errorCheck==0 %in this case there has been no errors in the input
    x=[a-h,x,b+h]; %add in x_(-1) and x_(n+1) to nodes
    for i=1:length(xhat) %for each value in xhat
        for j=1:length(c)
            q3(i)=q3(i)+c(j)*CubicSplineValue((xhat(i)-x(j))/h);
        end
    end
end

