clear;
rng(1234);

% peremeters
K = 8000;
N = 20;
h = 0;

V = K * rand(N, N);
% V = K * binornd(1, p, N, N);

[U, lam] = eig2d(V, h, 4);
W = solve2d(V, h);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!
[vx1, vx2, vw] = valley_line(w);

% plot landscape and valley line
figure
hold on
surf(x1, x2, w, 'LineStyle', 'None')
surf(x1, x2, zeros(size(x1)), 'LineStyle', 'None')
plot3(vx1, vx2, vw, 'r.', 'MarkerSize', 2)
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 12)
zticks([])
view(-150, 70)

% plot eigen mode and valley line
figure
hold on
mark = my_watershed(-w);
[vi, vj] = find(isnan(mark));
vu = zeros(size(vi));
for k = 1:4
    uk = my_nmlz(getval2d(U(:,k)));
    surf(x1, x2, uk, 'LineStyle', 'None');
    
    tuk = uk(:);
    vu = max(vu, tuk(vi+vj*size(uk,1)));
end
plot3(vx1, vx2, vu, 'r.', 'MarkerSize', 2);
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 12)
zticks([])
view(-150, 70)
