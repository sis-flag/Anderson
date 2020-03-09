clear;
rng(123);

% peremeters
h = 0;
p = 0.3;
V = rand(20);
V(V<p) = 0; V(V>=p) = 1; 

for K = [1e2]

    W = solveR2d(V * K, h);
    w = getval2d(W, 20);

    N = size(w,1);
    x = linspace(0, 1, N);
    [x2, x1] = meshgrid(x, x);

    mark = my_watershed(-w);
    [vi, vj] = find(isnan(mark));
    vx1 = (vi-1)/N; vx2 = (vj-1)/N;
    tw = w(:);
    vw = tw(vi+vj*N);
    
    figure()
    subplot(1,2,1)
    hold on
    
    pV = [V*K, zeros(20,1); zeros(1,21)];
    px = linspace(0, 1, 21);
    [px2, px1] = meshgrid(px, px);
    s = pcolor(px1,px2, pV);
    s.LineStyle = 'none';
    colorbar;
    
    fj = plot(vx1, vx2, 'r.');
    fj.MarkerSize = 3;
    
    subplot(1,2,2)
    hold on
    
    s = pcolor(x1,x2, w);
    s.LineStyle = 'none';
    colorbar;
    
    fj = plot(vx1, vx2, 'r.');
    fj.MarkerSize = 3;
end
