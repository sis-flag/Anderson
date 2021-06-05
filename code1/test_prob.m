clear
rng(0)
tic

%% peremeters
N_samp = 1000;
N1d = 50;
N2d = 15;

thes = 0.5;

%% relation to h
p = 0.5;
K = 1e4;
all_h = [0, 1e-1, 3e-1, 1e0, 3e0, 1e1, 3e1, 1e2, 3e2, 1e3, inf];

Pb_h = zeros(length(all_h), N_samp, 'logical');
for jh = 1:length(all_h)
    h = all_h(jh);
    for js = 1:N_samp
        V = K * binornd(1, p, N1d, 1);
        U = eig1d(V, h, 1);
        u = abs(getval1d(U));
        tPb = max(u(1), u(end)) / max(u) > thes;
        Pb_h(jh, js) = tPb;
    end
    
    fprintf('1d h %d / %d finished ', jh, length(all_h));
    toc
end

Pe_h = zeros(length(all_h), N_samp, 'logical');
Pc_h = zeros(length(all_h), N_samp, 'logical');
for jh = 1:length(all_h)
    h = all_h(jh);
    for js = 1:N_samp
        V = K * binornd(1, p, N2d, N2d);
        U = eig2d(V, h, 1);
        u = abs(getval2d(U));
        tPe = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(u(:)) > thes;
        tPc = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(u(:)) > thes;
        Pe_h(jh, js) = tPe;
        Pc_h(jh, js) = tPc;
    end
    
    fprintf('2d h %d / %d finished ', jh, length(all_h));
    toc
end

%% relation to p
all_p = 0.1:0.1:0.9;
K = 1e4;
h = 0.01;

Pb_p = zeros(length(all_p), N_samp, 'logical');
for jp = 1:length(all_p)
    p = all_p(jp);
    for js = 1:N_samp
        V = K * binornd(1, p, N1d, 1);
        U = eig1d(V, h, 1);
        u = abs(getval1d(U));
        tPb = max(u(1), u(end)) / max(u) > thes;
        Pb_p(jp, js) = tPb;
    end
    
    fprintf('1d p %d / %d finished ', jp, length(all_p));
    toc
end

Pe_p = zeros(length(all_p), N_samp, 'logical');
Pc_p = zeros(length(all_p), N_samp, 'logical');
for jp = 1:length(all_p)
    p = all_p(jp);
    for js = 1:N_samp
        V = K * binornd(1, p, N2d, N2d);
        U = eig2d(V, h, 1);
        u = abs(getval2d(U));
        tPe = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(u(:)) > thes;
        tPc = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(u(:)) > thes;
        Pe_p(jp, js) = tPe;
        Pc_p(jp, js) = tPc;
    end
    
    fprintf('2d p %d / %d finished ', jp, length(all_p));
    toc
end

%% relation to K
p = 0.5;
all_K = [1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5];
h = 0.01;

Pb_K = zeros(length(all_K), N_samp, 'logical');
for jk = 1:length(all_K)
    K = all_K(jk);
    for js = 1:N_samp
        V = K * binornd(1, p, N1d, 1);
        U = eig1d(V, h, 1);
        u = abs(getval1d(U));
        tPb = max(u(1), u(end)) / max(u) > thes;
        Pb_K(jk, js) = tPb;
    end
    
    fprintf('K-1d %d / %d finished ', jk, length(all_K));
    toc
end

Pe_K = zeros(length(all_K), N_samp, 'logical');
Pc_K = zeros(length(all_K), N_samp, 'logical');
for jk = 1:length(all_K)
    K = all_K(jk);
    for js = 1:N_samp
        V = K * binornd(1, p, N2d, N2d);
        U = eig2d(V, h, 1);
        u = abs(getval2d(U));
        tPe = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(u(:)) > thes;
        tPc = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(u(:)) > thes;
        Pe_K(jk, js) = tPe;
        Pc_K(jk, js) = tPc;
    end
    
    fprintf('2d K %d / %d finished ', jk, length(all_K));
    toc
end

%% save
save('Prob.mat', 'all_p', 'all_K', 'all_h', 'N_samp', 'N1d', 'N2d',...
    'Pb_h', 'Pb_p', 'Pb_K','Pe_h', 'Pe_p', 'Pe_K', 'Pc_h', 'Pc_p', 'Pc_K')

%% plot
% clear
load Prob

figure
hold on
plot(all_h(2:end-1), mean(Pb_h(2:end-1,:), 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h(2:end-1), mean(Pe_h(2:end-1,:), 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h(2:end-1), mean(Pc_h(2:end-1,:), 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
legend('P_b', 'P_e', 'P_c')
ylim([0, 1])
xlim([3e-2, 3e3])
xticks(all_h(2:2:end))
xlabel('h')
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

figure
hold on
plot(all_p, mean(Pb_p, 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_p, mean(Pe_p, 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_p, mean(Pc_p, 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
legend('P_b', 'P_e', 'P_c')
ylim([0, 1])
xticks(all_p(1:2:end))
xlabel('p')
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)

figure
hold on
plot(all_K, mean(Pb_K, 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_K, mean(Pe_K, 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_K, mean(Pc_K, 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
legend('P_b', 'P_e', 'P_c')
ylim([0, 1])
xlim([3, 3e5])
xticks(all_K(1:2:end))
xlabel('K')
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')
