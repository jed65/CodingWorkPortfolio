%% Part 1b)
[u,v]=CrankNicolsonForParabolicSystem(1,200,99,999);
x=linspace(0,1,99+1); %form vector of x values
t=linspace(0,200,999+1); %form vector of t values
[X,T]=meshgrid(x,t);
%Heatmaps
figure(1);
pcolor(X,T,u); %for u
shading flat
colorbar
xlabel('x');
ylabel('t');
title('Heatmap for u');
figure(2);
pcolor(X,T,v); %for v
shading flat
colorbar
xlabel('x');
ylabel('t');
title('Heatmap for v');
figure(3);
plot(x,u(1000,1:100),'b-');
hold on
plot(x,v(1000,1:100),'r-');
hold off
xlabel('x');
legend('u','v','Location','best');
title('Plot of u and v versus x at T=200');

%% Part 1c)
e=zeros(1,3);
N_x=[99,99*2,99*4,99*8];
h=1./N_x;
[u1,~]=CrankNicolsonForParabolicSystem(1,200,N_x(1),999); %call function 4 times to get numerical solutions for the different stepsizes
[u2,~]=CrankNicolsonForParabolicSystem(1,200,N_x(2),999);
[u3,~]=CrankNicolsonForParabolicSystem(1,200,N_x(3),999);
[u4,~]=CrankNicolsonForParabolicSystem(1,200,N_x(4),999);
e(1)=abs(h(2)*sum(u2(1000,:))-h(1)*sum(u1(1000,:))); %find e_j for j=1,2,3
e(2)=abs(h(3)*sum(u3(1000,:))-h(2)*sum(u2(1000,:)));
e(3)=abs(h(4)*sum(u4(1000,:))-h(3)*sum(u3(1000,:)));
figure(4);
h=h(2:4);
loglog(h,e);
xlabel('h');
ylabel('e');
title('Plot of e_{j} against h_{j}');

    

