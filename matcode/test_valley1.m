clear;

%% poteitial
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
saveas(gcf, '../report0/valley/V.png')

%% simulation
[U, lam] = eigR2d(V, 0, 50);
W = solveN2d(V, 0);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!

[vx1, vx2, vw] = valley_line(w);

%% plot valley line
figure();
hold on
cf = pcolor(x1, x2, w);
cf.LineStyle = 'none';
axis square
colorbar;
cf = plot(vx1, vx2, 'k.');
cf.MarkerSize = 3;
saveas(gcf, '../report0/valley/W.png')

%% plot eigen mode
for k = [1,2,3,5,7,10,20,30,50]
    evx1 = vx1(vw < 1/lam(k) );
    evx2 = vx2(vw < 1/lam(k) );
    
    uk = my_nmlz(getval2d(U(:,k)));
    
    figure();
    hold on
    cf = pcolor(x1, x2, uk);
    cf.LineStyle = 'none';
    caxis([-1,1])
    axis square
    colorbar;
    cf = plot(evx1, evx2, 'k.');
    cf.MarkerSize = 3;
    title(sprintf('eigen mode %d, \\lambda = %g', k, lam(k)))
    saveas(gcf, sprintf('../report0/valley/U(%d).png', k))
end

%% plot black-white
for k = [1,2,3,5,7,10,20,30]
    uk = my_nmlz(getval2d(U(:,k)));
    
    figure();
    subplot(1, 2, 1);
    cf = pcolor(x1, x2, uk);
    cf.LineStyle = 'none';
    caxis([-1,1])
    axis square    
    title(sprintf('eigen mode %d, \\lambda = %g', k, lam(k)))

    subplot(1, 2, 2);
    bw = double(w < 1/lam(k) );
    cf = pcolor(x1, x2, bw);
    cf.LineStyle = 'none';
    caxis([0, 1])
    axis square
    title('subregions')

    set(gcf,'position',[100 100 600 300]);
    saveas(gcf, sprintf('../report0/valley/B(%d).png', k))
end