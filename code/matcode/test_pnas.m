clear;
rng(0);

% peremeters
h = 0;
V = rand(20)*3000;

W = solveD2d(V);
w = getval2d(W, 20);

[U, lam] = eigD2d(V, 10, 100);
u1 = getval2d(U(:,1), 20); u1 = my_nmlz(u1);
u2 = getval2d(U(:,2), 20); u2 = my_nmlz(u2);
u3 = getval2d(U(:,3), 20); u3 = my_nmlz(u3);

N = size(w,1);
x = linspace(0, 1, N);
[x2, x1] = meshgrid(x, x);

mark = my_watershed(-w);
[vi, vj] = find(isnan(mark));
vx1 = (vi-1)/N; vx2 = (vj-1)/N;
tw = w(:);
vw = tw(vi+vj*N);

%% fig2-1
figure()
mesh(x1, x2, w);

%% fig2-2
figure();
hold on
s = pcolor(x1, x2, w);
s.LineStyle = 'none';
colorbar;
fj = plot(vx1, vx2, 'k.');
fj.MarkerSize = 3;

%% fig2-3
figure()
hold on
mesh(x1, x2, u1)
mesh(x1, x2, u2)
mesh(x1, x2, u3)
plot3(vx1, vx2, 0.1*ones(size(vx1)), 'r.');
zlim([0,1])

%% fig3
for k = [1, 2, 3]
    uk = getval2d(U(:,k), 20);
    uk = my_nmlz(uk);
    
    evx1 = vx1(vw<1/lam(k));
    evx2 = vx2(vw<1/lam(k));
    
    figure();
%     hold on
%     s = pcolor(x1, x2, uk);
%     s.LineStyle = 'none';
%     colorbar;
%     caxis([-1,1]);
    fj = plot(evx1, evx2, 'k.');
    fj.MarkerSize = 3;
    title(sprintf('eigen mode %d',k));
end
