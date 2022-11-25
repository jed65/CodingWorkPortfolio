function [x,y,t,u] = Upwind2D(Func_u0,a,N_x,N_y,M_t)
%Upwind2D.m: Taking inputs of Func_u0 (which finds u at initial time), the
%constant vector a=[a(1),a(2)] and step numbers N_x, N_y and M_t, this
%function implements the upwind scheme on the 2-dimensional scalar
%transport equation. The output u is a three-dimensional matrix, with each
%page (2-dimensional matrix) representing the solution at one time-step.

%Define step-sizes and discretisation vectors in each direction
dx=2*pi/N_x;
dy=2*pi/N_y;
dt=1/M_t;
lambda_x=dt/dx;
lambda_y=dt/dy;
x=linspace(0,2*pi,N_x+1);
y=linspace(0,2*pi,N_y+1);
t=linspace(0,1,M_t+1);

%Initialise u and apply function for initial condition
u=zeros(length(x),length(y),length(t)); %initialise u as matrix where each row represents a fixed time
u(:,:,1)=Func_u0(x,y); %call function to calculate u at initial time (this is given in the problem)

%Now apply the method, which depends on the sign of the elements in the
%vector a
if (a(1)>0)&&(a(2)>0) %case 1/4
    for n=1:M_t %fix time
        for k=1:N_y+1 %fix y
            for j=1:N_x+1 %finally fix x
                if (k==1) && (j==1) %deal with the four corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j,k,n)-u(N_x+1-j,k,n))-a(2)*lambda_y*(u(j,k,n)-u(j,N_y+1-k,n));
                elseif (k==1)&& (j==N_x+1)
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j,k,n)-u(j-1,k,n))-a(2)*lambda_y*(u(j,k,n)-u(j,N_y+1-k,n));
                elseif (k==N_y+1) && (j==1)
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j,k,n)-u(N_x+1-j,k,n))-a(2)*lambda_y*(u(j,k,n)-u(j,k-1,n));
                %elseif (k==N_y+1) && (j==N_x+1) NOT necessary
                    %u(j,k,n+1)=u(1,1,n+1);
                elseif (k==1) && ((j~=1)&&(j~=N_x+1)) %bottom horizontal, not corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j,k,n)-u(j-1,k,n))-a(2)*lambda_y*(u(j,k,n)-u(j,N_y+1-k,n));
                %elseif (k==N_y+1) && ((j~=1)&&(j~=N_x+1)) %top horizontal
                    %u(j,k,n+1)=u(j,1,n+1); %use periodic boundary condition
                elseif ((k~=1)&&(k~=N_y+1)) && (j==1) %left vertical, not corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j,k,n)-u(N_x+1-j,k,n))-a(2)*lambda_y*(u(j,k,n)-u(j,k-1,n));
                %elseif ((k~=1)&&(k~=N_y+1)) && (j==N_x+1) %right vertical, not corners
                    %u(j,k,n+1)=u(1,k,n+1); %use periodic boundary condition
                else %finally, consider general case
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j,k,n)-u(j-1,k,n))-a(2)*lambda_y*(u(j,k,n)-u(j,k-1,n));
                end
            end
        end
    end
elseif (a(1)>0)&&(a(2)<0) %case 2/4
       for n=1:M_t %fix time
        for k=1:N_y+1 %fix y
            for j=1:N_x+1 %finally fix x
                if (k==1) && (j==1) %deal with the four corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j,k,n)-u(N_x+1-j,k,n))-a(2)*lambda_y*(u(j,k+1,n)-u(j,k,n));
                elseif (k==1)&& (j==N_x+1)
                    u(j,k,n+1)=u(1,1,n+1); %use periodic boundary conditions for other three corners
                elseif (k==N_y+1) && (j==1)
                    u(j,k,n+1)=u(1,1,n+1);
                elseif (k==N_y+1) && (j==N_x+1)
                    u(j,k,n+1)=u(1,1,n+1);
                elseif (k==1) && ((j~=1)&&(j~=N_x+1)) %bottom horizontal, not corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j,k,n)-u(j-1,k,n))-a(2)*lambda_y*(u(j,k+1,n)-u(j,k,n));
                elseif (k==N_y+1) && ((j~=1)&&(j~=N_x+1)) %top horizontal, not corners
                    u(j,k,n+1)=u(j,1,n+1); %use periodic boundary condition
                elseif ((k~=1)&&(k~=N_y+1)) && (j==1) %left vertical, not corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j,k,n)-u(N_x+1-j,k,n))-a(2)*lambda_y*(u(j,k+1,n)-u(j,k,n));
                elseif ((k~=1)&&(k~=N_y+1)) && (j==N_x+1) %right vertical, not corners
                    u(j,k,n+1)=u(1,k,n+1); %use periodic boundary condition
                else %finally, consider general case
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j,k,n)-u(j-1,k,n))-a(2)*lambda_y*(u(j,k+1,n)-u(j,k,n));
                end
            end
        end
       end
elseif (a(1)<0)&&(a(2)>0) %case 3/4
    for n=1:M_t %fix time
        for k=1:N_y+1 %fix y
            for j=1:N_x+1 %finally fix x
                if (k==1) && (j==1) %deal with the four corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j+1,k,n)-u(j,k,n))-a(2)*lambda_y*(u(j,k,n)-u(j,N_y+1-k,n));
                elseif (k==1)&& (j==N_x+1)
                    u(j,k,n+1)=u(1,1,n+1); %use periodic boundary conditions for other three corners
                elseif (k==N_y+1) && (j==1)
                    u(j,k,n+1)=u(1,1,n+1);
                elseif (k==N_y+1) && (j==N_x+1)
                    u(j,k,n+1)=u(1,1,n+1);
                elseif (k==1) && ((j~=1)&&(j~=N_x+1)) %bottom horizontal, not corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j+1,k,n)-u(j,k,n))-a(2)*lambda_y*(u(j,k,n)-u(j,N_y+1-k,n));
                elseif (k==N_y+1) && ((j~=1)&&(j~=N_x+1)) %top horizontal, not corners
                    u(j,k,n+1)=u(j,1,n+1); %use periodic boundary condition
                elseif ((k~=1)&&(k~=N_y+1)) && (j==1) %left vertical, not corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j+1,k,n)-u(j,k,n))-a(2)*lambda_y*(u(j,k,n)-u(j,k-1,n));
                elseif ((k~=1)&&(k~=N_y+1)) && (j==N_x+1) %right vertical, not corners
                    u(j,k,n+1)=u(1,k,n+1); %use periodic boundary condition
                else %finally, consider general case
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j+1,k,n)-u(j,k,n))-a(2)*lambda_y*(u(j,k,n)-u(j,k-1,n));
                end
            end
        end
    end
else %a(1),a(2)<0 (both negative), case 4/4
    for n=1:M_t %fix time
        for k=1:N_y+1 %fix y
            for j=1:N_x+1 %finally fix x
                if (k==1) && (j==1) %deal with the four corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j+1,k,n)-u(j,k,n))-a(2)*lambda_y*(u(j,k+1,n)-u(j,k,n));
                elseif (k==1)&& (j==N_x+1)
                    u(j,k,n+1)=u(1,1,n+1); %use periodic boundary conditions for other three corners
                elseif (k==N_y+1) && (j==1)
                    u(j,k,n+1)=u(1,1,n+1);
                elseif (k==N_y+1) && (j==N_x+1)
                    u(j,k,n+1)=u(1,1,n+1);
                elseif (k==1) && ((j~=1)&&(j~=N_x+1)) %bottom horizontal, not corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j+1,k,n)-u(j,k,n))-a(2)*lambda_y*(u(j,k+1,n)-u(j,k,n));
                elseif (k==N_y+1) && ((j~=1)&&(j~=N_x+1)) %top horizontal, not corners
                    u(j,k,n+1)=u(j,1,n+1); %use periodic boundary condition
                elseif ((k~=1)&&(k~=N_y+1)) && (j==1) %left vertical, not corners
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j+1,k,n)-u(j,k,n))-a(2)*lambda_y*(u(j,k+1,n)-u(j,k,n));
                elseif ((k~=1)&&(k~=N_y+1)) && (j==N_x+1) %right vertical, not corners
                    u(j,k,n+1)=u(1,k,n+1); %use periodic boundary condition
                else %finally, consider general case
                    u(j,k,n+1)=u(j,k,n)-a(1)*lambda_x*(u(j+1,k,n)-u(j,k,n))-a(2)*lambda_y*(u(j,k+1,n)-u(j,k,n));
                end
            end
        end
    end
end
u=permute(u,[2 1 3]); %tranpose the matrix for each time, necessary so that the maximum ends up in right position (easier to do this than play with indices above)                     
                
                