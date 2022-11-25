function A = FormFiniteDifferenceMatrix(N,h_x,h_y)
%FormFiniteDifferenceMatrix.m: Takes inputs of the number of steps in each
%direction along with the step-size in each respective direction (which may
%not be equal) and forms the block tridiagonal matrix which is used to
%approximate the Poisson equation. The matrix A outputted is in sparse form
%in order to save storage and allow the use of large N.

%Initialise and define the matrices on the block tridiagonal
A=sparse((N-1)^2,(N-1)^2); %initialise sparse matrix 
Amd=gallery('tridiag',N-1,-1/h_x^2,2/(h_x^2)+2/(h_y^2),-1/h_x^2); %this is the block matrix on the diagonal (also sparse)
Asub=sparse(-1/h_y^2.*eye(N-1)); %these are the block matrices above and below the diagonal (also sparse)
Asup=sparse(-1/h_y^2.*eye(N-1));

%Now form the block tridiagonal matrix
A(1:N-1,1:N-1)=Amd; %assign the first row of matrices (main diagonal and the upper diagonal)
A(1:N-1,N:2*N-2)=Asup; 
central=N; %to keep track of which row of matrices we are at
for i=2:N-2 %for loop assigns matrices to all the matrix rows with 3 blocks (middle N-3 rows)
A(central:central+N-2,central-(N-1):central-1)=Asub; %lower diagonal matrix
A(central:central+N-2,central:central+N-2)=Amd; %main diagonal matrix
A(central:central+N-2,central+N-1:central+2*(N-1)-1)=Asup; %upper diagonal matrix
central=central+N-1; %update this to go to the next row of matrices
end
A((N-1)^2-(N-2):(N-1)^2,(N-1)^2-(N-2):(N-1)^2)=Amd; %finally, assign the final row of matrices (lower diagonal and main diagonal)
A((N-1)^2-(N-2):(N-1)^2,(N-1)^2-2*(N-2)-1:(N-1)^2-(N-2)-1)=Asub;


