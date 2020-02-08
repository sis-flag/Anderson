clear;
rng(0);

% peremeters
K = 1000;
V = rand(20)*K;
h = 0;

x0 = 0.51; y0 = 0.51;
beta = 1;
alpha = 200;

% Robin boundary
W0 = solveR2d(V, 0);
W1 = solveR2d(V, 0.1);
W2 = solveR2d(V, -0.1);
W3 = solveR2d(V, 0.5);
Wd = solveD2d(V);
% [U, lam] = eigR2d(V, h);

w1 = getval2d(W1, 20);
w2 = getval2d(W2, 20);
w3 = getval2d(W3, 20);
w0 = getval2d(W0, 20);
wd = getval2d(Wd, 20);
% u1 = getval2d(U(:,1), 20); u1 = my_nmlz(u1);
% u2 = getval2d(U(:,2), 20); u2 = my_nmlz(u2);
% u3 = getval2d(U(:,3), 20); u3 = my_nmlz(u3);

x = linspace(0,1,length(w1));
[x2, x1] = meshgrid(x, x); % caution!

%%

figure(); s = pcolor(x1, x2, w0); s.LineStyle='none'; colorbar; title('h/\beta=0')
figure(); s = pcolor(x1, x2, w1); s.LineStyle='none'; colorbar; title('h/\beta=0.1')
figure(); s = pcolor(x1, x2, w2); s.LineStyle='none'; colorbar; title('h/\beta=-0.1')
figure(); s = pcolor(x1, x2, w3); s.LineStyle='none'; colorbar; title('h/\beta=0.5')
figure(); s = pcolor(x1, x2, wd); s.LineStyle='none'; colorbar; title('Dirichlet')

% Dirichlet boundary
% W = solveD2d(V);
% [U, lam] = eigD2d(V);

% w = getval2d(W, 20);
% u1 = getval2d(U(:,1), 20); u1 = my_nmlz(u1);
% u2 = getval2d(U(:,2), 20); u2 = my_nmlz(u2);
% u3 = getval2d(U(:,3), 20); u3 = my_nmlz(u3);
% 
% x = linspace(0,1,length(w));
% [x2, x1] = meshgrid(x, x); % caution!
% 
% figure(); mesh(x1, x2, w);
% figure(); mesh(x1, x2, u1);
% figure(); mesh(x1, x2, u2);
% figure(); mesh(x1, x2, u3);
