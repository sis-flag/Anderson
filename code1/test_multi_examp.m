clear

% peremeters
K = 1e4;
V = K * [1,1,0,1,1,0,0,0,0,1,...
         0,1,1,1,0,1,0,0,1,1,...
         0,1,0,0,1,0,1,1,1,0,...
         1,0,0,0,0,1,1,0,1,1];
h = 0;
N = length(V);

[U, lam] = eig1d(V, h, 2);
W = solve1d(V, h);

figure
hold on
cl = bar(((1:N) - 0.5) / N, V);
cl.BarWidth = 1;
cl.LineStyle = 'None';
cl.FaceColor = [0.8, 0.8, 0.8];
cl = bar(((1:N) - 0.5) / N, -V);
cl.BarWidth = 1;
cl.LineStyle = 'None';
cl.FaceColor = [0.8, 0.8, 0.8];
w = getval1d(W);
x = linspace(0, 1, length(w));
plot(x, w, 'k', 'LineWidth', 1)
u1 = my_nmlz(getval1d(U(:,1)));
u2 = my_nmlz(getval1d(U(:,2)));
nu1 = my_nmlz(0.5 * u1 + 0.5 * u2);
nu2 = my_nmlz(0.5 * u2 - 0.5 * u1);
plot(x, nu1 / lam(1), 'LineWidth', 1)
plot(x, nu2 / lam(2), 'LineWidth', 1)
xlim([0, 1])
ylim([-2e-3, 2e-3])
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 16)

disp(lam)
