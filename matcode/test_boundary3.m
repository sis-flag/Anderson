clear;

%% example 1
rng(0);

K = 8000;
V = rand(1, 20);
V = K * V;

% plot potential
figure();
cf = bar(1/40:1/20:1, V);
cf.BarWidth = 1;
saveas(gcf, '../report0/boundary/V1.png')

[UR, lamR] = eigR1d(V, 0, 50);
[UD, lamD] = eigD1d(V, 50);

% plot eigenmodes
for k = [1,2,3,5,7,10,20]
    uRk = my_nmlz(getval1d(UR(:,k)));
    uDk = my_nmlz(getval1d(UD(:,k)));
    x = linspace(0,1,length(uRk));
    
    figure();
    hold on
    plot(x, uRk);
    plot(x, uDk);
    title(sprintf('eigen mode %d, \\lambda_D = %g, \\lambda_N = %g',...
                  k, lamD(k), lamR(k)))
    legend('Neumann boundary', 'Dirichlet boundary')
    saveas(gcf, sprintf('../report0/boundary/E1(%d).png', k))
end

figure();
hold on
plot(lamR, 'bo-');
plot(lamD, 'r*-');
xlabel('k')
ylabel('\lambda_k')
legend('Neumann boundary', 'Dirichlet boundary')
saveas(gcf, '../report0/boundary/lam1.png')


%% example 2
rng(0);

K = 8000;
V =[ 1, 5, 16, 12, 4, 18, 13, 7, 14, 5,...
    11, 17, 4, 19, 5, 13, 8, 10, 6, 2] / 20;
V = K * V;

% plot potential
figure();
cf = bar(1/40:1/20:1, V);
cf.BarWidth = 1;
saveas(gcf, '../report0/boundary/V2.png')

[UR, lamR] = eigR1d(V, 0, 50);
[UD, lamD] = eigD1d(V, 50);

% plot eigenmodes
for k = [1,2,3,5,7,10,20]
    uRk = my_nmlz(getval1d(UR(:,k)));
    uDk = my_nmlz(getval1d(UD(:,k)));
    x = linspace(0,1,length(uRk));
    
    figure();
    hold on
    plot(x, uRk);
    plot(x, uDk);
    title(sprintf('eigen mode %d, \\lambda_D = %g, \\lambda_N = %g',...
                  k, lamD(k), lamR(k)))
    legend('Neumann boundary', 'Dirichlet boundary')
    saveas(gcf, sprintf('../report0/boundary/E2(%d).png', k))
end

figure();
hold on
plot(lamR, 'bo-');
plot(lamD, 'r*-');
xlabel('k')
ylabel('\lambda_k')
legend('Neumann boundary', 'Dirichlet boundary')
saveas(gcf, '../report0/boundary/lam2.png')

%% example 3
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
saveas(gcf, '../report0/boundary/V3.png')

[UR, lamR] = eigR2d(V, 0, 50);
[UD, lamD] = eigD2d(V, 50);

u1 = getval2d(UR(:,1));
x = linspace(0, 1, size(u1,1));
[x2, x1] = meshgrid(x, x); % caution!

for k = [1,2,5,7,10,20]
    uRk = my_nmlz(getval2d(UR(:,k)));
    uDk = my_nmlz(getval2d(UD(:,k)));
    
    figure();
    subplot(1,2,1);
    cf = pcolor(x1, x2, uRk);
    cf.LineStyle = 'none';
    caxis([-1,1])
    axis square
    colorbar;
    title(sprintf('Neummann \\lambda_{%d} = %g', k, lamR(k)));
    
    subplot(1,2,2);
    cf = pcolor(x1, x2, uDk);
    cf.LineStyle = 'none';
    caxis([-1,1])
    axis square
    colorbar;
    title(sprintf('Dirichlet \\lambda_{%d} = %g', k, lamD(k)));
    
    set(gcf,'position',[100 100 600 300]);
    saveas(gcf, sprintf('../report0/boundary/E3(%d).png', k))
end

figure();
hold on
plot(lamR, 'bo-');
plot(lamD, 'r*-');
xlabel('k')
ylabel('\lambda_k')
legend('Neumann boundary', 'Dirichlet boundary')
saveas(gcf, '../report0/boundary/lam3.png')
