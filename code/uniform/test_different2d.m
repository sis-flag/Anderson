clear;
rng(0);

% peremeters
K = 8000;
N = 20;
p = 0.5;
h = 0;

V = K * rand(N, N);
% V = K * binornd(1, p, N, N);

[UD, lamD] = eig2d(V, inf, 50);
[UN, lamN] = eig2d(V, 0, 50);

x = linspace(0, 1, 20*N+1);
[x2, x1] = meshgrid(x, x); % caution!

for k = 1:4
    figure
    
    subplot(1, 2, 1)
    uk = my_nmlz(getval2d(UD(:,k)));
    cf = pcolor(x1, x2, uk);
    cf.LineStyle = 'None';
    caxis([-1, 1])
    axis square
    colorbar
    title(sprintf('eigenmode %d (Dirichlet), \\lambda_D = %g', k, lamD(k)))
    
    subplot(1, 2, 2)
    uk = my_nmlz(getval2d(UN(:,k)));
    cf = pcolor(x1, x2, uk);
    cf.LineStyle = 'None';
    caxis([-1, 1])
    axis square
    colorbar
    title(sprintf('eigenmode %d (Neumann), \\lambda_N = %g', k, lamN(k)))
    
    set(gcf,'Position',[500 500 900 300]);
end

figure
hold on
plot(lamD, 'o-')
plot(lamN, '*-')
legend('Dirichlet BC', 'Neumann BC')
xlabel('k')
ylabel('\lambda_k')
set(gcf,'Position',[500 500 400 300])
