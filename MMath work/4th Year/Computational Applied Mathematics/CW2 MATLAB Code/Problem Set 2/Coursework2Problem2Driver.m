%% Question 2b) code
f=@(x,y) 32.*pi.^2.*cos(4.*pi.*x).*cos(4.*pi.*y); %define the function f (found by applying minus Laplacian operator to u)
gridFunction2Norm=zeros(1,7); %initialise vectors
N_holder=zeros(1,7);
N=16; %assign first value of N
for i=1:7
[x,y,u]=FiniteDifferencePoisson(N,0,1,0,1,f,@g); %call function to get finite difference approximation
h_x=x(2)-x(1); %not necessary really because we have square domain but for completeness keep them separate
h_y=y(2)-y(1);
[X,Y]=meshgrid(x,y); 
true_u=cos(4.*pi.*X).*cos(4.*pi.*Y); %calculate true solution
true_u=reshape(true_u,(N+1)^2,1); %reshape into vector
gridFunction2Norm(i)=sqrt(h_x*h_y*sum(abs(true_u-u).^2)); %find grid function 2 norm
N_holder(i)=N;
if N==64 %in the case that N=64, save the three vectors from the function call above
    u_N64=u;
    x_N64=x;
    y_N64=y;
end
N=2*N; %update N by doubling
end
figure(1); %first figure is the log-log plot
loglog(N_holder,gridFunction2Norm);
xlabel('N');
ylabel('Error (grid function 2-norm)');
title('Log-log plot of error versus number of nodal points');
figure(2); %second figure is colour plot
[X,Y]=meshgrid(x_N64,y_N64); %form mesh 
U=reshape(u_N64,65,65); %reshape into matrix
pcolor(X,Y,U);
shading flat
colorbar
xlabel('x');
ylabel('y');
title('Colour plot of approximate solution when N=64');
format long g
ErrorTable=table(N_holder',gridFunction2Norm','VariableNames',{'N','Grid_function_2_norm'}); %also form a table from the results above
disp(ErrorTable);


%Define function for g (found by setting x or y equal to 0/1 from true
%solution)
function gValues=g(x,y)
    gValues=zeros(length(x),1);
    for i=1:length(x)
    if (x(i)==0)||(x(i)==1) %in this case, we are on bottom or top horizontal
        gValues(i)=cos(4*pi*y(i));
    elseif (y(i)==0)||(y(i)==1) %in this case, on left or right vertical
        gValues(i)=cos(4*pi*x(i));
    end
    end
end
