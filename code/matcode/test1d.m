clear;
rng(0);

% peremeters
K = 3000;
V =[ 1, 5, 16, 12, 4, 18, 13, 7, 14, 5,...
    11, 17, 4, 19, 5, 13, 8, 10, 6, 2] / 20;
V = [1,0,1,0,0,0,1,1,0,0];
V = V*K;
bar(V)

h = 1;
beta = 100;

% Robin boundary
[U, lam] = eigR1d(V, h, 10, 20);
u1 = getval1d(U(:,1), 50); u1 = my_nmlz(u1);
u2 = getval1d(U(:,2), 50); u2 = my_nmlz(u2);
u3 = getval1d(U(:,3), 50); u3 = my_nmlz(u3);
u4 = getval1d(U(:,4), 50); u4 = my_nmlz(u4);


W = solveR1d(V, h/beta);
w = getval1d(W, 50);

x = linspace(0,1,length(w));
figure();
hold on
plot(x, w, 'k')
plot(x, u1 / (lam(1) + 1/max(w) + beta) )
plot(x, u2 / (lam(2) + 1/max(w) + beta) )
plot(x, u3 / (lam(3) + 1/max(w) + beta) )
plot(x, u4 / (lam(4) + 1/max(w) + beta) )
legend('w','u1','u2','u3','u4')
title(sprintf('K=%g \\beta=%g h=%g no enforce',K, beta, h))

% Dirichlet boundary
[U, lam] = eigD1d(V, 10, 20);
u1 = getval1d(U(:,1), 50); u1 = my_nmlz(u1);
u2 = getval1d(U(:,2), 50); u2 = my_nmlz(u2);
u3 = getval1d(U(:,3), 50); u3 = my_nmlz(u3);
u4 = getval1d(U(:,4), 50); u4 = my_nmlz(u4);

W = solveD1d(V);
w = getval1d(W, 50);

x = linspace(0,1,length(w));
figure();
hold on
plot(x, w, 'k')
plot(x, u1 / lam(1) )
plot(x, u2 / lam(2) )
plot(x, u3 / lam(3) )
plot(x, u4 / lam(4) )
legend('w','u1','u2','u3','u4')
title(sprintf('K=%g Dirichlet boundary',K))