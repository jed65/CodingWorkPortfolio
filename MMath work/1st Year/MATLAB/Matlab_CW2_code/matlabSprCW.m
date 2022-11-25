%% MATLAB Spring Coursework 2018-19

%% Part (d)
%This code plots C(x) and S(x) on the same axis
x=linspace(-8,8,500); %gives us 500 points (as required in the question) in the interval between -8 and 8
y1=fresnelC(x,50);%calls fresnelC to be calculated for the x specified above with 50 strips 
set(groot,'DefaultTextInterpreter','latex');
plot(x,y1,'k-');%plots the graph of C(x)
y2=fresnelS(x,50);%calls fresnelS to be calculated for the x specified above with 50 strips
hold on %holds the above plot command so the second curve can be added
plot(x,y2,'b-');
hold off
xlabel('$x$');
ylabel('$y$');
title('Plots of C($x$) and S($x$)');
axis([-8,8,-1,1]);%sets axis limits to match the x values specified and the limits of each curve
legend('C(x)','S(x)','Location','Best')%adds legends to each curve, placing the key in an area of the plot not occupied by the curves


%% Part (e)
%This code produces a plot of S(x) against C(x)
x=linspace(-8,8,500);%as before gives us 500 points in [-8,8]
y1=fresnelC(x,50);%calls fresnelC to be calculated for these x values 
y2=fresnelS(x,50);%calls fresnelS to be calculated for these x values
set(groot,'DefaultTextInterpreter','latex');
plot(y1,y2,'r-');%plots y1(C(x)) on horizontal axis with y2(S(x)) on vertical axis
xlabel('C($x$)');
ylabel('S($x$)');
title('Plot of S($x$) against C($x$)');

%% Part (f)
% At n=30, there is far greater oscillation at the ends of the curve (-8 to -6 and 6 to 8), whereas the larger n are much more similar (better approximations) converging quickly. The largest differences are at the local extrema at the ends of the curves. In the region of x from -6 to 6 all the approximations are extremely close together (similar)

x=linspace(-8,8,500);%Provides points on the required interval
y1=fresnelS(x,90);%Calls fresnelS to be calculated (90 strips)
plot(x,y1,'color',[0.9100 0.4100 0.1700],'LineWidth',2);%Plot the n=90 approximation with orange colour and thicker line(found this more obvious to see) 
set(groot,'DefaultTextInterpreter','latex');
y2=fresnelS(x,70);%Calls fresnelS to be calculated (70 strips)
y3=fresnelS(x,50);%Calls fresnelS to be calculated (50 strips)
y4=fresnelS(x,30);%Calls fresnelS to be calculated (30 strips)
hold on %allows the 3 other curves to be plotted on the same graph
plot(x,y2,'m-','LineWidth',2);%Plot n=70 approx. in magenta
plot(x,y3,'g-','LineWidth',2);%Plot n=50 approx. in green
plot(x,y4,'c-','LineWidth',2);%Plot n=30 approx. in cyan
hold off
xlabel('$x$');
ylabel('S($x$)');
title('Plots of S($x$) for varying numbers of strips n');
xlim([-8,8]);
legend('n=90','n=70','n=50','n=30','Location','Best')%Places the key for legends in a place where we can still see the curves