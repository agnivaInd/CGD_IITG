u = readmatrix("output.csv");

u1 = u(1,:);
u2 = u(2,:);

x = -9.99:0.02:9.99;

figure(1)
plot(x,u1);
title("u1")

figure(2)
plot(x,u2);
title("u2")