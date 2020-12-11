clear;

%% example 1
rng(0);
V = rand(20);
K = [1e4, 1e5, 1e6, 1e7];

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
saveas(gcf, '../report0/potential/P0V.png')

for i = 1:length(K)
    KV = K(i) * V;
    
    [U, lam] = eigR2d(KV, 0, 6);
    W = solveN2d(KV, 0);
    
    w = getval2d(W);
    x = linspace(0, 1, size(w,1));
    [x2, x1] = meshgrid(x, x); % caution!
    
    % plot landscape
    figure();
    hold on
    cf = pcolor(x1, x2, w);
    cf.LineStyle = 'none';
    axis square
    colorbar;
    title(sprintf('landscape K = %g', K(i)));
    saveas(gcf, sprintf('../report0/potential/P0K%dW.png', i))
    
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
    saveas(gcf, sprintf('../report0/potential/P0K%dU.png', i))
    
end


%% example 2
rng(123);
p = 0.3;
V = rand(20);
V(V<p) = 0; V(V>=p) = 1;

K = [1e4, 1e5, 1e6, 1e7];

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
saveas(gcf, '../report0/potential/P3V.png')

for i = 1:length(K)
    KV = K(i) * V;
    
    [U, lam] = eigR2d(KV, 0, 6);
    W = solveN2d(KV, 0);
    
    w = getval2d(W);
    x = linspace(0, 1, size(w,1));
    [x2, x1] = meshgrid(x, x); % caution!
    
    % plot landscape
    figure();
    hold on
    cf = pcolor(x1, x2, w);
    cf.LineStyle = 'none';
    axis square
    colorbar;
    title(sprintf('landscape K = %g', K(i)));
    saveas(gcf, sprintf('../report0/potential/P3K%dW.png', i))
    
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
    saveas(gcf, sprintf('../report0/potential/P3K%dU.png', i))
    
end


%% example 3
rng(456);
p = 0.5;
V = rand(20);
V(V<p) = 0; V(V>=p) = 1;

K = [1e4, 1e5, 1e6, 1e7];

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
saveas(gcf, '../report0/potential/P5V.png')

for i = 1:length(K)
    KV = K(i) * V;
    
    [U, lam] = eigR2d(KV, 0, 6);
    W = solveN2d(KV, 0);
    
    w = getval2d(W);
    x = linspace(0, 1, size(w,1));
    [x2, x1] = meshgrid(x, x); % caution!
    
    % plot landscape
    figure();
    hold on
    cf = pcolor(x1, x2, w);
    cf.LineStyle = 'none';
    axis square
    colorbar;
    title(sprintf('landscape K = %g', K(i)));
    saveas(gcf, sprintf('../report0/potential/P5K%dW.png', i))
    
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
    saveas(gcf, sprintf('../report0/potential/P5K%dU.png', i))
    
end


%% example 4
rng(789);
p = 0.7;
V = rand(20);
V(V<p) = 0; V(V>=p) = 1;

K = [1e4, 1e5, 1e6, 1e7];

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
saveas(gcf, '../report0/potential/P7V.png')

for i = 1:length(K)
    KV = K(i) * V;
    
    [U, lam] = eigR2d(KV, 0, 6);
    W = solveN2d(KV, 0);
    
    w = getval2d(W);
    x = linspace(0, 1, size(w,1));
    [x2, x1] = meshgrid(x, x); % caution!
    
    % plot landscape
    figure();
    hold on
    cf = pcolor(x1, x2, w);
    cf.LineStyle = 'none';
    axis square
    colorbar;
    title(sprintf('landscape K = %g', K(i)));
    saveas(gcf, sprintf('../report0/potential/P7K%dW.png', i))
    
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
    saveas(gcf, sprintf('../report0/potential/P7K%dU.png', i))
    
end