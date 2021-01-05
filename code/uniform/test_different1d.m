clear;
rng(0);

% peremeters
K = 8000;
N = 30;
p = 0.5;

V = K * rand(1, N);
% V = K * binornd(1, p, N, 1);

% plot potential
figure
M = length(V);
cf = bar(((1:M) - 0.5) / M, V);
cf.BarWidth = 1;
set(gcf, 'Position', [500 500 400 300])

[UD, lamD] = eig1d(V, inf, 50);
[UN, lamN] = eig1d(V, 0, 50);

for k = [1, 3, 10, 15]
    figure
    hold on

    uDk = my_nmlz(getval1d(UD(:,k)));
    uNk = my_nmlz(getval1d(UN(:,k)));
    x = linspace(0, 1, length(uDk));
    
    plot(x, uDk, x, uNk)
    legend('Dirichlet BC', 'Neumann BC')
    title(sprintf('eigenmode %d, \\lambda_D = %g, \\lambda_N = %g', k, lamD(k), lamN(k)))
    set(gcf,'Position',[500 500 400 300])
end

figure
hold on
plot(lamD, 'o-')
plot(lamN, '*-')
legend('Dirichlet BC', 'Neumann BC')
xlabel('k')
ylabel('\lambda_k')
set(gcf,'Position',[500 500 400 300])
