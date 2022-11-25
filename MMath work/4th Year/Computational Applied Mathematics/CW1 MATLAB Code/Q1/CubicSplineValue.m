function B=CubicSplineValue(x)
%CubicSplineValue.m: Function to evaluate the standard cubic spline
%function at point x.

if x<=-2
    B=0;
elseif (x>-2)&&(x<=-1)
    B=(x+2)^3;
elseif (x>-1)&&(x<=0)
    B=1+3*(x+1)+3*(x+1)^2-3*(1+x)^3;
elseif (x>0)&&(x<=1)
    B=1+3*(1-x)+3*(1-x)^2-3*(1-x)^3;
elseif (x>1)&&(x<=2)
    B=(2-x)^3;
else
    B=0;
end
