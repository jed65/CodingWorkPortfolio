function [u,v] = CrankNicolsonForParabolicSystem(L,T,N_x,M_t)
%CrankNicolsonForParabolicSystem.m: Function taking inputs of L (the
%maximum of x), T (maximum in time), N_x (number of space-steps) and M_t
%(number of time-steps), which implements the Crank Nicolson method given
%the initial conditions from part b) and using first-order finite
%difference for the Neumann boundary conditions. The outputs are the
%vectors u and v, approximates for the solution of the system.


%Define step-sizes and discretisation grid vectors plus constants
h=L/N_x; %define space step-size
dt=T/M_t; %define time step-size
x=linspace(0,L,N_x+1); %form vector of x values
t=linspace(0,T,M_t+1); %form vector of t values
beta_v=0.1;
beta_u=0.01*beta_v;

%Assign the initial conditions
u=zeros(length(t),length(x)); %initialise matrices for u and v values 
v=zeros(length(t),length(x)); 
for i=1:length(x)
    if x(i)<=0.9 %in this case, value of u is 0.1
        u(1,i)=0.1; %first row in u is initial time (zero)
    else
        u(1,i)=2; %otherwise, value of u is 2
    end
    v(1,i)=1.97; %this is case for all x values (found integral by hand)
end

%Now form matrix K of Crank-Nicolson, identity matrix and multipliers
K=diag(2*ones(1,N_x-1))+diag(-1*ones(1,N_x-2),1)+diag(-1*ones(1,N_x-2),-1);
K(1,1)=1; K(N_x-1,N_x-1)=1; %these come from using the first order finite difference approx (1/h)*(u1-u0)=0 and (1/h)*(uN-u_(N-1))=0 on Neumann BCs
I=eye(N_x-1); %identity matrix
lambda_u=(beta_u*dt)/(h^2); %multipliers here
lambda_v=(beta_v*dt)/(h^2);
%Find inverse of matrix on left-hand side
inverse_for_u=inv(I+(lambda_u/2).*K);
inverse_for_v=inv(I+(lambda_v/2).*K);

%Now perform the method
for i=1:length(t)-1
    f=v(i,2:length(x)-1).*(0.067+(u(i,2:length(x)-1).^2)./(1+u(i,2:length(x)-1).^2))-u(i,2:length(x)-1); %find forcing term for middle N-1 points
    f=dt.*f'; %scale it by a factor of dt as it is in part a)
    u(i+1,2:length(x)-1)=[(I-(lambda_u/2).*K)*u(i,2:length(x)-1)'+f]'; %find u at next time step (split over two lines)
    u(i+1,2:length(x)-1)=[inverse_for_u*u(i+1,2:length(x)-1)']';
    v(i+1,2:length(x)-1)=[(I-(lambda_v/2).*K)*v(i,2:length(x)-1)'-f]'; %find v at next time step
    v(i+1,2:length(x)-1)=[inverse_for_v*v(i+1,2:length(x)-1)']';
end
%Now recall that from Neumann BCs we have u0=u1 (v0=v1) and uN=u_(N-1) (vN=v_(N-1))
u(1:length(t),1)=u(1:length(t),2); %u0=u1
u(1:length(t),N_x+1)=u(1:length(t),N_x); %uN=u_(N-1)
v(1:length(t),1)=v(1:length(t),2); %repeat for v
v(1:length(t),N_x+1)=v(1:length(t),N_x); 
    