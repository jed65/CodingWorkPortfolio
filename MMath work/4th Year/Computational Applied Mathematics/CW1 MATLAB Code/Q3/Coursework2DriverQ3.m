%% Question 3a) (TEST)
F={@(x) x(1)^2-x(2); @(x) x(1)^2+x(2)^2-2};
TestJacobi=ComputeJacobian(F,[2,2],0.005);
disp(TestJacobi); %works well for this test problem!

%% Question 3b) 
format long g
epsilon=[0.005,0.05,0.5];
F={@(x) x(1)^2-x(2); @(x) x(1)^2+x(2)^2-2};
for i=1:length(epsilon)
    disp(['Epsilon now set equal to: ',num2str(epsilon(i))]);
    [x,x5]=ImplementNewtonMethod(F,[2;2],1e-10,epsilon(i),5);
    disp(['Value of x at iteration 5 using epsilon=',num2str(epsilon(i)),' is:']);
    disp(x5);
end

