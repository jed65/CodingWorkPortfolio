%% Question 3c)
%Set up problem
a=1;
Func_u0=@(x) 1-cos(x);
N_x=20;
M_t=10;

%Lax-Friedrich 1D
[x,~,u_LF]=LaxFriedrich1D(Func_u0,a,N_x,M_t);

%Upwind 1D
[~,~,u_UP]=Upwind1D(Func_u0,a,N_x,M_t);

%True solution
true_tHalf=1-cos(x-a*0.5); %true at t=1/2
true_tOne=1-cos(x-a*1); %true at t=1

%Plots for comparison
figure(1);
subplot(2,1,1);
plot(x,true_tHalf,'k'); %plot true solution at t=1/2
hold on
plot(x,u_LF(6,:),'b'); %Lax-Friedrich approx. at t=1/2
plot(x,u_UP(6,:),'r'); %Upwind approx. at t=1/2
hold off
xlim([0,2*pi]);
xlabel('x');
ylabel('u(x,0.5)');
title('Comparison at t=1/2');
legend('True','Lax-Friedrich','Upwind','Location','best');
subplot(2,1,2);
plot(x,true_tOne,'k'); %plot true solution at t=1
hold on
plot(x,u_LF(11,:),'b'); %Lax-Friedrich approx. at t=1
plot(x,u_UP(11,:),'r'); %Upwind approx. at t=1
hold off
xlim([0,2*pi]);
xlabel('x');
ylabel('u(x,1)');
title('Comparison at t=1');
legend('True','Lax-Friedrich','Upwind','Location','best');

%Now find the error at the final time
dx=x(2)-x(1); %define step-size
error_LF=sqrt(dx*sum(abs((true_tOne-u_LF(11,:)).^2))); %find grid function 2-norm at t=1 for Lax-Friedrich
error_UP=sqrt(dx*sum(abs((true_tOne-u_UP(11,:)).^2))); %find grid function 2-norm at t=1 for Upwind
disp(['Error in Lax-Friedrich method at t=1 is: ',num2str(error_LF)]); %display for each
disp(['Error in Upwind method at t=1 is: ',num2str(error_UP)]);

%% Question 3d)
%Set up problem
a=[1 0.5];
N_x=20;
N_y=20;
M_t=10;

%Upwind 2D
[x,y,~,u_UP_2D] = Upwind2D(@u0,a,N_x,N_y,M_t);

%True solution at final time t=1
[X,Y]=meshgrid(x,y);
u_true_2D=(1-cos(X-a(1))).*(1-cos(Y-a(2)));

%Now produce surfaces for comparison
figure(2);
subplot(2,1,1);
surf(x,y,u_UP_2D(:,:,11));
xlabel('x');
ylabel('y');
zlabel('u');
title('Approximate u(x,y,1) by upwind method');
subplot(2,1,2);
surf(x,y,u_true_2D);
xlabel('x');
ylabel('y');
zlabel('u');
title('True u(x,y,1)');

%Now find error at final time
dx=x(2)-x(1); %define step-sizes in each direction
dy=y(2)-y(1);
error_UP_2D=sqrt(dx*dy*sum(abs((u_true_2D-u_UP_2D(:,:,11)).^2),'all')); %find grid function 2-norm at t=1
disp(['Error in Upwind method at t=1 for 2D problem is: ',num2str(error_UP_2D)]);

%Check the maximum is in similar/same place
[row, col] = find(ismember(u_true_2D, max(u_true_2D(:))));
disp(['Maximum of true solution is: ',num2str(u_true_2D(row,col)),' at position (x,y)=(',num2str(x(col)),',',num2str(y(row)),')']);
[row2, col2] = find(ismember(u_UP_2D(:,:,11), max(u_UP_2D(:,:,11),[],'all')));
disp(['Maximum of approximate solution is: ',num2str(u_UP_2D(row2,col2,11)),' at position (x,y)=(',num2str(x(col2)),',',num2str(y(row2)),')']);

function initial_u=u0(x,y) %function to calculate u(x,y,0)
[X,Y]=meshgrid(x,y);
initial_u=(1-cos(X)).*(1-cos(Y));        
end
