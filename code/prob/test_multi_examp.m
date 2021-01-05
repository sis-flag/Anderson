clear;
rng(0);

% peremeters
K = 5000;
h = 0;

V = K * [1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1];

% plot potential
figure
M = length(V);
cf = bar(((1:M) - 0.5) / M, V);
cf.BarWidth = 1;
set(gcf, 'Position', [500 500 400 300]);

[U, lam] = eig1d(V, h, 2);
W = solve1d(V, h);

figure
hold on
w = getval1d(W);
x = linspace(0, 1, length(w));
plot(x, w, 'k')
set(gcf,'Position',[500 500 400 300]);

figure
hold on
u1 = my_nmlz(getval1d(U(:,1)));
plot(x, u1, 'r.')
u2 = my_nmlz(getval1d(U(:,2)));
plot(x, u2, 'b-')
legend('eigenmode 1', 'eigenmode2')
set(gcf,'Position',[500 500 400 300]);

disp(lam)