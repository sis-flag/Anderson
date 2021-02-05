clear
rng(0)

% peremeters
K = 1e4;
N = 20;
p = 0.8;
h = 0;

V = K * binornd(1, p, N, N);

W = solve2d(V, h);
w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!
[vx1, vx2, vw] = valley_line(w);

% plot valley line and potential
figure
hold on
M = size(V, 1);
pV = [V, zeros(M,1); zeros(1,M+1)];
px = linspace(0, 1, M+1);
[px2, px1] = meshgrid(px, px);
cl = pcolor(px1, px2, pV);
cl.LineStyle = 'None';
plot(vx1, vx2, 'r.', 'MarkerSize', 3);
caxis([0, K])
axis square

temp1 = pcolor(100*ones(2), 100*ones(2), zeros(2));
temp1.LineStyle = 'None';
temp2 = pcolor(100*ones(2), 100*ones(2), K*ones(2));
temp2.LineStyle = 'None';
temp3 = plot(100*ones(2,1), 100*ones(2,1), 'r-');
temp3.LineWidth = 1;
legend([temp1, temp2, temp3], ...
       {'V(x) = 0', 'V(x) = 1', 'valley line'}, ...
       'Location', 'NorthEastOut',...
       'EdgeColor', 'None')
   
xlim([0,1])
ylim([0,1])
set(gcf, 'Position', [300 300 500 300])
set(gca, 'FontSize', 16)
