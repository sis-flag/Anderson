clear;

%% example 3
rng(0);

% peremeters
K = 8000;
V = rand(1, 20);
V = K * V;

h = 10;

% plot potential
figure();
cf = bar((1:20)-0.5, V);
cf.BarWidth = 1;
saveas(gcf, '../report0/parameter/V3.png')

[U, lam] = eigR1d(V, h, 6);

L2e = zeros(9, 6);
all_beta = [100, 180, 300, 560, 1000, 1800, 3000, 5600, 10000];
for kb = 1:length(all_beta)
    
    beta = all_beta(kb);
    
    w = getval1d(solveR1d(V, h/beta));
    x = linspace(0, 1, length(w));

    figure();
    hold on
    plot(x, w, 'k')
    for k = 1:6
        uk = my_nmlz(getval1d(U(:,k)));
        plot(x, uk / (lam(k) + beta) )
        
        te = w - uk / (lam(k) + beta);
        L2e(kb, k) = sqrt(mean(mean(te.^2))) / sqrt(mean(mean(w.^2)));
    end
    title(sprintf('\\beta = %g', beta))
    saveas(gcf, sprintf('../report0/parameter/W3B%d.png', beta))
    close()
end

figure()
semilogx(all_beta, L2e, '-o')
xlabel('\beta')
ylabel('gap')
saveas(gcf, '../report0/parameter/gap3.png')

%% example 4
rng(0);

% peremeters
K = 8000;
V = rand(1, 20);
p = 0.3;
V(V<p) = 0; V(V>=p) = 1;
V = K * V;

h = 10;

% plot potential
figure();
cf = bar((1:20)-0.5, V);
cf.BarWidth = 1;
saveas(gcf, '../report0/parameter/V4.png')

[U, lam] = eigR1d(V, h, 6);

L2e = zeros(9, 6);
all_beta = [100, 180, 300, 560, 1000, 1800, 3000, 5600, 10000];
for kb = 1:length(all_beta)
    
    beta = all_beta(kb);
    
    w = getval1d(solveR1d(V, h/beta));
    x = linspace(0, 1, length(w));

    figure();
    hold on
    plot(x, w, 'k')
    for k = 1:6
        uk = my_nmlz(getval1d(U(:,k)));
        plot(x, uk / (lam(k) + beta) )
        
        te = w - uk / (lam(k) + beta);
        L2e(kb, k) = sqrt(mean(mean(te.^2))) / sqrt(mean(mean(w.^2)));
    end
    title(sprintf('\\beta = %g', beta))
    saveas(gcf, sprintf('../report0/parameter/W4B%d.png', beta))
    close()
end

figure()
semilogx(all_beta, L2e, '-o')
xlabel('\beta')
ylabel('gap')
saveas(gcf, '../report0/parameter/gap4.png')