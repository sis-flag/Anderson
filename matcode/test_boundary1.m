clear;

%% 1d potential
rng(0);

K = 8000;
V = rand(1, 20);
V = K * V;

% plot potential
figure();
cf = bar(1/40:1/20:1, V);
cf.BarWidth = 1;
saveas(gcf, '../report0/boundary/V1d.png')

%% 1d Robin boundary
h = 1;
beta = 30;

[U, lam] = eigR1d(V, h, 6);
W = solveR1d(V, h/beta);

w = getval1d(W);
x = linspace(0,1,length(w));

figure();
hold on
plot(x, w, 'k')
for k = 1:4
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, uk / (lam(k) + 1/max(w) + beta) )
end
title('Robin boundary')
saveas(gcf, '../report0/boundary/R1d.png')

%% 1d Neumann boundary
[U, lam] = eigR1d(V, 0, 6);
W = solveR1d(V, 0);

w = getval1d(W);
x = linspace(0,1,length(w));

figure();
hold on
plot(x, w, 'k')
for k = 1:4
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, uk / (lam(k) + 1/max(w)) )
end
title('Neumann boundary')
saveas(gcf, '../report0/boundary/N1d.png')

%% 1d Dirichlet boundary
[U, lam] = eigD1d(V, 6);
W = solveD1d(V);

w = getval1d(W);
x = linspace(0,1,length(w));

figure();
hold on
plot(x, w, 'k')
for k = 1:4
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, uk / lam(k) )
end
title('Dirichlet boundary')
saveas(gcf, '../report0/boundary/D1d.png')

%% 2d potential
rng(0);

K = 8000;
V = rand(20);
V = K * V;

% plot potential
figure();
M = size(V, 1);
pV = [V, zeros(M,1); zeros(1,M+1)];
px = linspace(0, 1, M+1);
[px2, px1] = meshgrid(px, px);
cf = pcolor(px1, px2, pV);
cf.LineStyle = 'none';
caxis([0, K])
axis square
colorbar;
saveas(gcf, '../report0/boundary/V2d.png')

%% 2d Robin boundary
h = 1;
beta = 30;

[U, lam] = eigR2d(V, h, 6);
W = solveR2d(V, h/beta);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!

% plot valley line
figure();
hold on
cf = pcolor(x1, x2, w);
cf.LineStyle = 'none';
axis square
colorbar;
[vx1, vx2] = valley_line(w);
cf = plot(vx1, vx2, 'k.');
cf.MarkerSize = 3;
title('landscape on Robin boundary')
saveas(gcf, '../report0/boundary/Rw2d.png')

% plot eigen mode
figure();
for k = 1:4
    subplot(2, 2, k);
    uk = my_nmlz(getval2d(U(:,k)));
    cf = pcolor(x1, x2, uk);
    cf.LineStyle = 'none';
    caxis([-1,1])
    axis square
    colorbar;
    title(sprintf('eigenmode %d', k));
end
saveas(gcf, '../report0/boundary/Ru2d.png')

%% 2d Neumann boundary
[U, lam] = eigR2d(V, 0, 6);
W = solveR2d(V, 0);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!

% plot valley line
figure();
hold on
cf = pcolor(x1, x2, w);
cf.LineStyle = 'none';
axis square
colorbar;
[vx1, vx2] = valley_line(w);
cf = plot(vx1, vx2, 'k.');
cf.MarkerSize = 3;
title('landscape on Neumann boundary')
saveas(gcf, '../report0/boundary/Nw2d.png')

% plot eigen mode
figure();
for k = 1:4
    subplot(2, 2, k);
    uk = my_nmlz(getval2d(U(:,k)));
    cf = pcolor(x1, x2, uk);
    cf.LineStyle = 'none';
    caxis([-1,1])
    axis square
    colorbar;
    title(sprintf('eigenmode %d', k));
end
saveas(gcf, '../report0/boundary/Nu2d.png')

%% 2d Dirichlet boundary
[U, lam] = eigD2d(V, 6);
W = solveD2d(V);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!

% plot valley line
figure();
hold on
cf = pcolor(x1, x2, w);
cf.LineStyle = 'none';
axis square
colorbar;
[vx1, vx2] = valley_line(w);
cf = plot(vx1, vx2, 'k.');
cf.MarkerSize = 3;
title('landscape on Dirichlet boundary')
saveas(gcf, '../report0/boundary/Dw2d.png')

% plot eigen mode
figure();
for k = 1:4
    subplot(2, 2, k);
    uk = my_nmlz(getval2d(U(:,k)));
    cf = pcolor(x1, x2, uk);
    cf.LineStyle = 'none';
    caxis([-1,1])
    axis square
    colorbar;
    title(sprintf('eigenmode %d', k));
end
saveas(gcf, '../report0/boundary/Du2d.png')