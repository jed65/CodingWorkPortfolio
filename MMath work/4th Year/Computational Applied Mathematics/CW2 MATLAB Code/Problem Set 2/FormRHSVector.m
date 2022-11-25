function F = FormRHSVector(x,y,f,g)
%FormRHSVector.m: Takes inputs of x and y (which contain both the step-size
%and number of steps), f (function on right of Poisson equation) and g
%(boundary condition), and forms the (N-1)^2x1 vector for the linear system
%which approximates the Poisson equation by central finite difference.

%Use inputs to find quantities we need to form F
h_x=x(2)-x(1); %find step-sizes from x and y vectors
h_y=y(2)-y(1);
N=length(x)-1; %find N from x vector

%Now fill the vector f
F=zeros((N-1)^2,1); %initialise the vector F
counter=1;
for i=1:N-1 %begin by filling F with evaluations of the function f
    F(counter:counter+N-2)=f(x(2:N),y(i+1).*ones(1,N-1));
    counter=counter+N-1;
end
%Add appropriate elements of g 
F(1:N-1)=F(1:N-1)+(1/h_y^2).*g(x(2:N),y(1).*ones(1,N-1)); %add elements to bottom row in mesh
F((N-1)^2-(N-2):(N-1)^2)=F((N-1)^2-(N-2):(N-1)^2)+(1/h_y^2).*g(x(2:N),y(N+1).*ones(1,N-1)); %add elements to top row in mesh
end1=1; %represents the start (one end) of a row in the mesh
end2=N-1; %represents the end (other end) of a row in the mesh
for i=1:N-1 %add appropriate evaluation of g to the ends of each row 
    F(end1)=F(end1)+(1/h_x^2)*g(x(1),y(i+1)); %beginning of row gets g_(0,i)
    F(end2)=F(end2)+(1/h_x^2)*g(x(N+1),y(i+1)); %end of row gets g_(N,i)
    end1=end1+N-1; %update for next row
    end2=end2+N-1;
end

