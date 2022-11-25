function [y] = fresnelS(x,n)
%fresnelS This function finds a Simpson's rule approximation for the area under the graph f(t)=sin(t^2) with varying interval of integration length x and lower limit zero
f=@(t) sin(t.^2);%calls the function sin(t.^2) indirectly with dummy variable t
l=length(x);%counts the number of elements in the input vector x
y=zeros(1,l);%standardises y as a zero vector of length l just as x is a vector of length l
for e=1:l%e tells us which element in x is being calculated, setting the interval for e as integers between and including 1 and l ensures all elements in x are used
    y(e)=simpson(f,0,x(e),n);%calls our simpson function to compute the area approximation under f with the element x(e) as the upper limit and n number of strips, then assigns the corresponding y element with this approximation for the output
end
