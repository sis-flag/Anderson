clear;

K = 1e8;
V = K * ones(5);
V(1,2) = 0;
V(1,3) = 0;
V(2,3) = 0;
V(1,4) = 0;

figure
[UN, lamN] = eig2d(V, 0, 1);
uN = getval2d(UN, 100);
x = linspace(0, 1, size(uN,1));
[x2, x1] = meshgrid(x, x); % caution!
% uN(abs(uN) < 1e-3) = nan;
cl = surfl(x1, x2, -uN, 'light');
cl(1).LineStyle = 'None';
title(['eigen value:', num2str(lamN)])
grid off
box off
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
set(gcf, 'Position', [300 300 400 400])
set(gca, 'FontSize', 16)
set(gca, 'ZColor', 'None')

figure
[U2D, lamD] = eig2d(V, inf, 1);
u2D = getval2d(U2D, 100);
x = linspace(0, 1, size(u2D,1));
[x2, x1] = meshgrid(x, x); % caution!
% uD(abs(uD) < 1e-3) = nan;
cl = surfl(x1, x2, -u2D, 'light');
cl(1).LineStyle = 'None';
title(['eigen value:', num2str(lamD)])
grid off
box off
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
set(gcf, 'Position', [300 300 400 400])
set(gca, 'FontSize', 16)
set(gca, 'ZColor', 'None')

K = 1e8;
V = K * ones(5);
V(2:3,2:4) = 0;
V(1,3) = 0;
V(4,3) = 0;

figure
[U2D, lam2D] = eig2d(V, inf, 1);
u2D = getval2d(U2D, 100);
x = linspace(0, 1, size(u2D,1));
[x2, x1] = meshgrid(x, x); % caution!
% uD(abs(uD) < 1e-3) = nan;
cl = surfl(x1, x2, -u2D, 'light');
cl(1).LineStyle = 'None';
title(['eigen value:', num2str(lam2D)])
grid off
box off
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
set(gcf, 'Position', [300 300 400 400])
set(gca, 'FontSize', 16)
set(gca, 'ZColor', 'None')