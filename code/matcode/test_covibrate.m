clear;

V = [1,0,1,0,0,0,1,1,0,0] * 200;
% bar(V)

h = -1;
beta = 10;

[U, lam] = eigR1d(V, h);
u1 = getval1d(U(:,1), 50); u1 = my_nmlz(u1);

W = solveR1d(V, h/beta);
w = getval1d(W, 50);

tV = V(8:end);
[tU, tlam] = eig_temp(tV, h);
tu1 = getval1d(tU(:,1), 50); tu1 = my_nmlz(tu1);

figure();
hold on
x = linspace(0,1,length(w));
plot(x, u1)
x = linspace(0.7,1,length(tu1));
plot(x, tu1)

disp(lam(1:2))
disp(tlam(1:2))