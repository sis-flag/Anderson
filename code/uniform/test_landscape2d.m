clear;
rng(456);

% peremeters
K = 8000;
N = 20;
p = 0.5;
h = 0;

V = K * rand(N, N);
% V = K * binornd(1, p, N, N);

% % plot potential 3d
% figure
% M = size(V, 1);
% px = ((1:M) - 0.5) / M;
% bar3(px, V);
% set(gcf, 'Position', [500 500 400 300]);
% view(-70, 70);

% plot potential 2d
figure
M = size(V, 1);
pV = [V, zeros(M,1); zeros(1,M+1)];
px = linspace(0, 1, M+1);
[px2, px1] = meshgrid(px, px);
cf = pcolor(px1, px2, pV);
cf.LineStyle = 'None';
caxis([0, K])
axis square
colorbar;
set(gcf, 'Position', [500 500 400 300]);

[U, lam] = eig2d(V, h, 4);
W = solve2d(V, h);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!

% plot landscape
figure
hold on
cf = surf(x1, x2, w);
cf.LineStyle = 'None';
cf = surf(x1, x2, zeros(size(x1)));
cf.LineStyle = 'None';
set(gcf, 'Position', [500 500 400 300]);
view(-70, 70)

% plot valley line
figure();
hold on
cf = pcolor(x1, x2, w);
cf.LineStyle = 'none';
axis square
colorbar
[vx1, vx2, vw] = valley_line(w);
cf = plot(vx1, vx2, 'r.');
cf.MarkerSize = 2;
set(gcf, 'Position', [500 500 400 300]);

% plot eigen mode and valley line
figure
hold on

mark = my_watershed(-w);
[vi, vj] = find(isnan(mark));

vu = zeros(size(vi));
for k = 1:4
    uk = my_nmlz(getval2d(U(:,k)));
    cf = surf(x1, x2, uk);
    cf.LineStyle = 'None';
    
    tuk = uk(:);
    vu = max(vu, tuk(vi+vj*size(uk,1)));
end
zlim([0, 1])
cf = plot3(vx1, vx2, vu, 'r.');
cf.MarkerSize = 2;
set(gcf, 'Position', [500 500 400 300]);
view(-70, 70)
