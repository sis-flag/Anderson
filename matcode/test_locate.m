clear;

p = 0.3;
V = rand(20);
V(V<p) = 0; V(V>=p) = 1;

% plot potential
figure();
M = size(V, 1);
pV = [V, zeros(M,1); zeros(1,M+1)];
px = linspace(0, 1, M+1);
[px2, px1] = meshgrid(px, px);
cf = pcolor(px1, px2, pV);
cf.LineStyle = 'none';
caxis([0, 1])
axis square
colorbar;

mark = mark_branch(V);
