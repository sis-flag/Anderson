clear
rng(0)

% peremeters
K = 10000;
N = 20;
p = 0.7;
h = 0;

% V = K * rand(N, N);
V = K * binornd(1, p, N, N);

[U, lam] = eig2d(V, h, 4);
W = solve2d(V, h);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!

% plot valley line and potential
figure
hold on
M = size(V, 1);
pV = [V, zeros(M,1); zeros(1,M+1)];
px = linspace(0, 1, M+1);
[px2, px1] = meshgrid(px, px);
cf = pcolor(px1, px2, pV);
cf.LineStyle = 'None';
[vx1, vx2, vw] = valley_line(w);
cf = plot(vx1, vx2, 'r.');
cf.MarkerSize = 2;
caxis([0, K])
axis square
set(gcf, 'Position', [500 500 300 300])

% % plot eigen mode
% figure
% for k = 1:4
%     subplot(2, 2, k)
%     uk = my_nmlz(getval2d(U(:,k)));
%     cf = pcolor(x1, x2, uk);
%     cf.LineStyle = 'None';
%     caxis([-1, 1])
%     colorbar
% end
% set(gcf, 'Position', [500 500 400 300])
