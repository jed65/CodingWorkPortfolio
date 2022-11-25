function S = simpson(f, a, b, n)
% SIMPSON   This function outputs the Simpson's rule approximation for a user inputted function f between the limits a and b, with n strips (n must be even)
if mod(n,2) == 1 %if user input of n congruent to 1 modulo 2 (i.e. if n is odd)
    error('The number of strips must be even!');%Simpson's rule approximates using quadratic curves,finding a quadratic equation requires 3 known points therefore to approximate the whole curve n must be even, thus the function produces an error when n is odd  
end %function continues to run when n is even
h = (b-a)/n; %provides the width of the strips between the limits b and a on the horizontal axis

S = f(a)+f(b); %standardises variable S for the for loop with a value that appears in the final approximation equation
for k = 1:2:n %the values of k in this for loop are odd numbers between 1 and n
    S = S + 4*f(a+k*h); %adds four times each f(a+k*h) for odd additions of h to limit a (on horizontal axis)
end
for k=2:2:n-1 %the values of k in this for loop are the even numbers between 1 and n-1
    S = S + 2*f(a+k*h); %adds twice each f(a+k*h) for even additions of h to the limit a
end
S = S*h/3; %multiplies the output variable S by h/3 to give the final numerical approximation for the area under the curve according to Simpsons rule
