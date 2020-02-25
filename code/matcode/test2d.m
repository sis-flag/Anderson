clear;
rng(0);

% peremeters
K = 1000;
V = rand(20)*K;
h = 0;

x0 = 0.51; y0 = 0.51;
beta = 1;
alpha = 200;

% Dirichlet boundary
W = solveR2d(V, h);
[U, lam] = eigR2d(V, h);

w = getval2d(W, 20);
u1 = getval2d(U(:,1), 20); u1 = my_nmlz(u1);
u2 = getval2d(U(:,2), 20); u2 = my_nmlz(u2);
u3 = getval2d(U(:,3), 20); u3 = my_nmlz(u3);

x = linspace(0,1,length(w));
[x2, x1] = meshgrid(x, x); % caution!

figure(); mesh(x1, x2, w);
figure(); mesh(x1, x2, u1);
figure(); mesh(x1, x2, u2);
figure(); mesh(x1, x2, u3);

% Dirichlet boundary
W = solveD2d(V);
[U, lam] = eigD2d(V);

w = getval2d(W, 20);
u1 = getval2d(U(:,1), 20); u1 = my_nmlz(u1);
u2 = getval2d(U(:,2), 20); u2 = my_nmlz(u2);
u3 = getval2d(U(:,3), 20); u3 = my_nmlz(u3);

x = linspace(0,1,length(w));
[x2, x1] = meshgrid(x, x); % caution!

figure(); mesh(x1, x2, w);
figure(); mesh(x1, x2, u1);
figure(); mesh(x1, x2, u2);
figure(); mesh(x1, x2, u3);
