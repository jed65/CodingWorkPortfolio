function [x,y,u] = FiniteDifferencePoisson(N,min_x,max_x,min_y,max_y,f,g)
%FiniteDifferencePoisson.m: Taking inputs of the number of space-steps N,
%constants describing the domain, f (function on right of Poisson equation)
%and g (boundary condition function), assigns values within approximate u
%at the boundaries then calls FormFiniteDifferenceMatrix.m and
%FormRHSVector.m to construct the linear system. The linear system is then
%solved for the unknown nodes to give an approximate solution u.

%Find step-size and discretisation vectors
h_x=(max_x-min_x)/N; %step-size to get N+1 points in x-direction
h_y=(max_y-min_y)/N; %step-size to get N+1 points in y-direction
x=linspace(min_x,max_x,N+1); %vector of x values
y=linspace(min_y,max_y,N+1); %vector of y values

%Initialise at the boundaries
u=zeros((N+1)^2,1); %initialise matrix to hold nodes
u(1:N+1)=g(x,y(1)*ones(1,N+1)); %assign boundary condition to bottom horizontal
row_start=N+2;
row_end=2*(N+1);
for i=1:N %add the start and end of each row 
    u(row_start)=g(x(1),y(i+1));
    u(row_end)=g(x(N+1),y(i+1));
    row_start=row_start+N+1;
    row_end=row_end+N+1;
end
u(2+N*(N+1):(N+1)^2-1)=g(x(2:N),y(N+1)*ones(1,N-1)); %assign boundary condition to top horizontal

%Construct and solve the linear system
A=FormFiniteDifferenceMatrix(N,h_x,h_y); %call function to form general A matrix
F=FormRHSVector(x,y,f,g); %call function to form F vector (on the RHS)
uMiddle=A\F; %Multiply by the inverse of A on each side to solve the system and get the (N-1)^2 unknown nodes

%Incorporate this solution (middle (N-1)^2 set of nodes) into the appropriate parts of vector u
counter=1; 
initial=N+3;
for i=1:N-1 
    u(initial:initial+N-2)=uMiddle(counter:counter+N-2); %need to fill N-1 zeroes on each N+1 row, do this here
    initial=initial+N+1; %next two lines update these terms
    counter=counter+N-1;
end





