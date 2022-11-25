function c = ComputeLSQCoefficients(FuncIn,n)
%ComputeLSQCoefficients.m: Function that calculates the least-squares
%coefficients for the least-square approximation of a function FuncIn 
%(which acts on the domain [-1,1]) using the Legendre polynomials as a 
%basis. The input n is the number of coefficients requested by the user
%(plus one which is the coefficient for the zero polynomial).
%The coefficients themselves are formed using the weighted inner product
%and the fact that the Legendre basis is orthogonal.

c=zeros(1,n+1); %initialise
for k=0:n
    f_phi_product=@(x) FuncIn(x).*legendreP(k,x); %form anonymous function to integrate over (for numerator of coefficient)
    phi_phi_product=@(x) legendreP(k,x).*legendreP(k,x); %form anonymous function to integrate over (for denominator of coefficient)
    c(k+1)=GaussQuadrature(f_phi_product)/GaussQuadrature(phi_phi_product); %use the GaussQuadrature function I wrote in a) to integrate over [-1,1]
end

