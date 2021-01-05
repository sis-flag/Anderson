clear;
rng(0);

% peremeters
K = 8000;
N = 30;
p = 0.5;
h = 10;

V = K * rand(1, N);
% V = K * binornd(1, p, N, 1);

% plot potential
figure
M = length(V);
cf = bar(((1:M) - 0.5) / M, V);
cf.BarWidth = 1;
set(gcf, 'Position', [500 500 400 300])

[U, lam] = eig1d(V, h, 4);
W = solve1d(V, h);

figure
hold on
w = getval1d(W);
x = linspace(0, 1, length(w));
plot(x, w, 'k')
for k = 1:4
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, uk / lam(k))
end
set(gcf, 'Position', [500 500 400 300])