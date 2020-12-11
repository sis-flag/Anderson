clear;
rng(0);

% peremeters
rng(0);

K = 100;
p = 0.3;
V = rand(10);
V(V<p) = 0; V(V>=p) = 1;
V = K * V;

[U, lam] = eigR2d(V, 0, 100);
W = solveN2d(V, 0);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!

[vx1, vx2, vw] = valley_line(w);
%% fig1-1
figure()
tx = linspace(0, 1, size(V,1));
cf = bar3(tx, V);
for n=1:numel(cf)
    cdata=get(cf(n),'zdata');
    set(cf(n),'cdata',cdata,'facecolor','interp')
end

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
for k = 1:5
    uk = my_nmlz(getval2d(U(:,k)));
    mesh(x1, x2, uk);
end
plot3(vx1, vx2, 0.001*ones(length(vx1),1), 'r.');
zlim([0,1])

%% fig3
for k = [1, 2, 3, 10, 20, 30, 50, 75, 100]
    
    uk = my_nmlz(getval2d(U(:,k)));
    
    evx1 = vx1(vw<1/(lam(k)+1/max(max(w))));
    evx2 = vx2(vw<1/(lam(k)+1/max(max(w))));
    
    figure();
    hold on
    s = pcolor(x1, x2, uk);
    s.LineStyle = 'none';
    colorbar;
    caxis([-1,1]);
    fj = plot(evx1, evx2, 'k.');
    fj.MarkerSize = 3;
    title(sprintf('eigen mode %d \\lambda=%g',k, lam(k)));
end
