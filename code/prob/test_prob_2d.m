clear
rng(0)
tic

%% peremeters
N_samp = 1000;
N = 15;

%% relation to h
p = 0.5;
K = 1e3;
all_h = [0, 1e-1, 3e-1, 1e0, 3e0, 1e1, 3e1, 1e2, 3e2, 1e3, inf];

Pe_h = zeros(length(all_h), N_samp);
Pc_h = zeros(length(all_h), N_samp);
for jh = 1:length(all_h)
    h = all_h(jh);
    for js = 1:N_samp
        V = K * binornd(1, p, N, N);
        U = eig2d(V, h, 1);
        u = abs(getval2d(U));
        tPe = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(u(:));
        tPc = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(u(:));
        Pe_h(jh, js) = tPe;
        Pc_h(jh, js) = tPc;
    end
    
    fprintf('2d h %d / %d finished ', jh, length(all_h));
    toc
end

%% relation to p
all_p = 0.1:0.1:0.9;
K = 1e3;
h = 1;

Pe_p = zeros(length(all_p), N_samp);
Pc_p = zeros(length(all_p), N_samp);
for jp = 1:length(all_p)
    p = all_p(jp);
    for js = 1:N_samp
        V = K * binornd(1, p, N, N);
        U = eig2d(V, h, 1);
        u = abs(getval2d(U));
        tPe = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(u(:));
        tPc = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(u(:));
        Pe_p(jp, js) = tPe;
        Pc_p(jp, js) = tPc;
    end
    
    fprintf('2d p %d / %d finished ', jp, length(all_p));
    toc
end

%% relation to K
p = 0.5;
all_K = [1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5];
h = 1;

Pe_K = zeros(length(all_K), N_samp);
Pc_K = zeros(length(all_K), N_samp);
for jk = 1:length(all_K)
    K = all_K(jk);
    for js = 1:N_samp
        V = K * binornd(1, p, N, N);
        U = eig2d(V, h, 1);
        u = abs(getval2d(U));
        tPe = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(u(:));
        tPc = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(u(:));
        Pe_K(jk, js) = tPe;
        Pc_K(jk, js) = tPc;
    end
    
    fprintf('2d K %d / %d finished ', jk, length(all_K));
    toc
end

%% save
save('Prob2d.mat', 'all_p', 'all_K', 'all_h', 'N_samp', 'N',...
    'Pe_h', 'Pe_p', 'Pe_K', 'Pc_h', 'Pc_p', 'Pc_K')

%% plot
load Prob2d

figure
semilogx(all_h(2:end-1), mean(Pe_h(2:end-1,:), 2), '*-')
ylim([0, 1])
xlabel('h')
ylabel('P_e')
set(gcf, 'Position', [500 500 400 300])

figure
semilogx(all_h(2:end-1), mean(Pc_h(2:end-1,:), 2), '*-')
ylim([0, 1])
xlabel('h')
ylabel('P_c')
set(gcf, 'Position', [500 500 400 300])

figure
plot(all_p, mean(Pe_p, 2), '*-')
ylim([0, 1])
xlabel('p')
ylabel('P_e')
set(gcf, 'Position', [500 500 400 300])

figure
plot(all_p, mean(Pc_p, 2), '*-')
ylim([0, 1])
xlabel('p')
ylabel('P_c')
set(gcf, 'Position', [500 500 400 300])

figure
semilogx(all_K, mean(Pe_K, 2), '*-')
ylim([0, 1])
xlabel('K')
ylabel('P_e')
set(gcf, 'Position', [500 500 400 300])

figure
semilogx(all_K, mean(Pc_K, 2), '*-')
ylim([0, 1])
xlabel('K')
ylabel('P_c')
set(gcf, 'Position', [500 500 400 300])
