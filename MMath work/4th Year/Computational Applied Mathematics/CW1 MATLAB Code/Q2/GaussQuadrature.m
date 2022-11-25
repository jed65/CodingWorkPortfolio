function GQEstimate = GaussQuadrature(FuncIn)
% GaussQuadrature.m: Takes an input of a function FuncIn which takes one
% input and calculates an estimate for its integral over the interval
% [-1,1] using the Gauss Quadrature method with n=10 points.


%First obtain nodes and weights for n=10 points
syms x
nodes=vpasolve(legendreP(10,x)==0); %nodes are roots of Legendre polynomial
syms x
expr=legendreP(10,x);
derivativeLegendre=eval(subs(diff(expr,x,1),x,nodes)); %evaluate derivative of 10th order Legendre polynomial
previousLegendre=eval(subs(legendreP(9,x),x,nodes)); %evaluate previous (9th order) Legendre polynomial
weights=zeros(1,10); %initialise
for i=1:10
    weights(i)=2/(10*previousLegendre(i)*derivativeLegendre(i)); %assign weights
end

%Next form the Gauss Quadrature estimate
GQEstimate=0; %initialise
for i=1:10
    GQEstimate=GQEstimate+weights(i)*FuncIn(nodes(i)); %add w_i*f(x_i)
end
GQEstimate=double(GQEstimate); %convert from sym to double
