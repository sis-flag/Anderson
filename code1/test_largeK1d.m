clear

% peremeters
K = 1e3;
h = 0;
V = K * [1,0,0,1,0,1,1,0,1,0,...
         1,1,0,0,0,1,0,1,1,0,...
         1,1,0,0,0,0,0,0,0,1,...%%%%
         1,0,1,0,1,0,1,0,0,0,...
         1,1,1,0,1,0,0,0,1,1,...
         1,0,1,0,0,0,1,1,1,0,...
         1,0,1,1,1,0,0,0,0,1,...
         0,0,0,1,0,0,0,0,0,0,...%%%%%
         1,1,0,1,1,1,1,0,1,1,...
         1,1,0,0,1,1,1,0,1,1];
N = length(V);

[U, lam] = eig1d(V, h, 4);
W = solve1d(V, h);

figure
hold on
cl = bar(((1:N) - 0.5) / N, V);
cl.BarWidth = 1;
cl.LineStyle = 'None';
cl.FaceColor = [0.8, 0.8, 0.8];
w = getval1d(W);
x = linspace(0, 1, length(w));
plot(x, w, 'k', 'LineWidth', 1)
for k = 1:2
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, abs(uk) / lam(k), 'LineWidth', 1)
end
xlim([0, 1])
ylim([0, max(w)*1.01])
title(['K = ', num2str(K)])
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 16)
