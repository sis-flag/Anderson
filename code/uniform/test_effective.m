clear;
rng(456);

% peremeters
K = 8000;
N = 20;
p = 0.5;
h = 0;

V = K * rand(N, N);
% V = K * binornd(1, p, N, N);

[U, lam] = eig2d(V, h, 100);
W = solve2d(V, h);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!
[vx1, vx2, vw] = valley_line(w);

% plot eigen mode and valley line
for k = [1,2,3,5,10,20,30,50,70,100]
    uk = my_nmlz(getval2d(U(:,k)));
    
    figure
    hold on
    
    cf = pcolor(x1, x2, uk);
    axis square
    colorbar;
    cf.LineStyle = 'None';
    caxis([-1, 1])
    
    evx1 = vx1(vw<1/lam(k));
    evx2 = vx2(vw<1/lam(k));
    cf = plot(evx1, evx2, 'k.');
    cf.MarkerSize = 2;
    
    title(sprintf('eigen mode %d', k));
    set(gcf, 'Position', [500 500 400 300]);
end

