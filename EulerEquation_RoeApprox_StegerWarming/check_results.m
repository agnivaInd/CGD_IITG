given = readmatrix("StandardSod_0_15_transposed.csv");
x_loc = given(1,:);
rho = given(2,:);
u = given(3,:);
p = given(4,:);
ie = given(5,:);

roe_without = readmatrix("Roe_WITHOUT_Entropy_Fix_Output.csv");
u11 = roe_without(1,:);
u12 = roe_without(2,:);
u13 = roe_without(3,:);
u14 = roe_without(4,:);

roe_with = readmatrix("Roe_WITH_Entropy_Fix_Output.csv");
u21 = roe_with(1,:);
u22 = roe_with(2,:);
u23 = roe_with(3,:);
u24 = roe_with(4,:);

steger = readmatrix("Steger_Warming_output.csv");
u31 = steger(1,:);
u32 = steger(2,:);
u33 = steger(3,:);
u34 = steger(4,:);

count_roe_without = numel(u11);
start_val = 0;
end_val = 1;
stepsize = (end_val - start_val)/count_roe_without;
x_roe_without = (start_val + stepsize/2):stepsize:(end_val - stepsize/2);

count_roe_with = numel(u21);
start_val = 0;
end_val = 1;
stepsize = (end_val - start_val)/count_roe_with;
x_roe_with = (start_val + stepsize/2):stepsize:(end_val - stepsize/2);

count_steger = numel(u31);
start_val = 0;
end_val = 1;
stepsize = (end_val - start_val)/count_steger;
x_steger = (start_val + stepsize/2):stepsize:(end_val - stepsize/2);

figure(1)
subplot(2,2,1)
plot(x_loc,rho)
hold on
plot(x_roe_without,u11,'ro',MarkerSize=2);
hold off
title("Roe WITHOUT Entropy Fix - Density")
legend('Analytical','Numerical','Location','southwest')
subplot(2,2,2)
plot(x_loc,p)
hold on
plot(x_roe_without,u12,'ro',MarkerSize=2);
hold off
title("Roe WITHOUT Entropy Fix - Pressure")
legend('Analytical','Numerical','Location','southwest')
subplot(2,2,3)
plot(x_loc,u)
hold on
plot(x_roe_without,u13,'ro',MarkerSize=2);
hold off
title("Roe WITHOUT Entropy Fix - Velocity")
legend('Analytical','Numerical','Location','southwest')
subplot(2,2,4)
plot(x_loc,ie)
hold on
plot(x_roe_without,u14,'ro',MarkerSize=2);
hold off
title("Roe WITHOUT Entropy Fix - Internal Energy")
legend('Analytical','Numerical','Location','southwest')

figure(2)
subplot(2,2,1)
plot(x_loc,rho)
hold on
plot(x_roe_with,u21,'ro',MarkerSize=2);
hold off
title("Roe WITH Entropy Fix - Density")
legend('Analytical','Numerical','Location','southwest')
subplot(2,2,2)
plot(x_loc,p)
hold on
plot(x_roe_with,u22,'ro',MarkerSize=2);
hold off
title("Roe WITH Entropy Fix - Pressure")
legend('Analytical','Numerical','Location','southwest')
subplot(2,2,3)
plot(x_loc,u)
hold on
plot(x_roe_with,u23,'ro',MarkerSize=2);
hold off
title("Roe WITH Entropy Fix - Velocity")
legend('Analytical','Numerical','Location','southwest')
subplot(2,2,4)
plot(x_loc,ie)
hold on
plot(x_roe_with,u24,'ro',MarkerSize=2);
hold off
title("Roe WITH Entropy Fix - Internal Energy")
legend('Analytical','Numerical','Location','southwest')

figure(3)
subplot(2,2,1)
plot(x_loc,rho)
hold on
plot(x_steger,u31,'ro',MarkerSize=2);
hold off
title("Steger Warming - Density")
legend('Analytical','Numerical','Location','southwest')
subplot(2,2,2)
plot(x_loc,p)
hold on
plot(x_steger,u32,'ro',MarkerSize=2);
hold off
title("Steger Warming - Pressure")
legend('Analytical','Numerical','Location','southwest')
subplot(2,2,3)
plot(x_loc,u)
hold on
plot(x_steger,u33,'ro',MarkerSize=2);
hold off
title("Steger Warming - Velocity")
legend('Analytical','Numerical','Location','southwest')
subplot(2,2,4)
plot(x_loc,ie)
hold on
plot(x_steger,u34,'ro',MarkerSize=2);
hold off
title("Steger Warming - Internal Energy")
legend('Analytical','Numerical','Location','southwest')