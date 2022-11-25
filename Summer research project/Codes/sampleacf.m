function [acf,denom] = sampleacf(Data,LagMax,mu)
%sampleacf.m: Calculates the sample autocorrelation from the inputted data
%for lags ranging from 0 to the input LagMax.

acf=zeros(1,LagMax+1); %Initialise the output value
denom=0; %Initialise values
T=length(Data); 

for t=1:T %Calculate denominator in sample acf (sum of squares)
    term=(Data(t)-mu)^2;
    denom=denom+term;
end

for k=0:1:LagMax %Calculate the autocorrelation for different lags
    num=0;
    for i=k+1:T
        numterm=(Data(i)-mu)*(Data(i-k)-mu);
        num=num+numterm;
    end
    acf(k+1)=num/denom;
end

