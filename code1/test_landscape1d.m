clear;
rng(0);

% peremeters
K = 8000;
N = 30;
h = 0;

V = K * rand(1, N);
% V = K * binornd(1, p, N, 1);

[U, lam] = eig1d(V, h, 4);
W = solve1d(V, h);

figure
hold on
w = getval1d(W);
x = linspace(0, 1, length(w));
plot(x, w, 'k', 'LineWidth', 1)
for k = 1:4
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, uk / lam(k), 'LineWidth', 1)
end
ylim([0, 6e-4])
xlim([0, 1])
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 16)
