clear
rng(0)
tic

%% peremeters
N_samp = 1000;
N = 50;

%% relation to h
p = 0.5;
K = 1e3;
all_h = [0, 1e-1, 3e-1, 1e0, 3e0, 1e1, 3e1, 1e2, 3e2, 1e3, inf];

Pb_h = zeros(length(all_h), N_samp);
for jh = 1:length(all_h)
    h = all_h(jh);
    for js = 1:N_samp
        V = K * binornd(1, p, N, 1);
        U = eig1d(V, h, 1);
        u = abs(getval1d(U));
        tPb = max(u(1), u(end)) / max(u);
        Pb_h(jh, js) = tPb;
    end
    
    fprintf('1d h %d / %d finished ', jh, length(all_h));
    toc
end

%% relation to p
all_p = 0.1:0.1:0.9;
K = 1e3;
h = 1;

Pb_p = zeros(length(all_p), N_samp);
for jp = 1:length(all_p)
    p = all_p(jp);
    for js = 1:N_samp
        V = K * binornd(1, p, N, 1);
        U = eig1d(V, h, 1);
        u = abs(getval1d(U));
        tPb = max(u(1), u(end)) / max(u);
        Pb_p(jp, js) = tPb;
    end
    
    fprintf('1d p %d / %d finished ', jp, length(all_p));
    toc
end

%% relation to K
p = 0.5;
all_K = [1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5];
h = 1;

Pb_K = zeros(length(all_K), N_samp);
for jk = 1:length(all_K)
    K = all_K(jk);
    for js = 1:N_samp
        V = K * binornd(1, p, N, 1);
        U = eig1d(V, h, 1);
        u = abs(getval1d(U));
        tPb = max(u(1), u(end)) / max(u);
        Pb_K(jk, js) = tPb;
    end
    
    fprintf('K-1d %d / %d finished ', jk, length(all_K));
    toc
end

%% save
save('Prob1d.mat', 'all_p', 'all_K', 'all_h', 'N_samp', 'N',...
    'Pb_h', 'Pb_p', 'Pb_K')

%% plot
load Prob1d

figure
semilogx(all_h(2:end-1), mean(Pb_h(2:end-1,:), 2), '*-')
ylim([0, 1])
xlabel('h')
ylabel('P_b')
set(gcf, 'Position', [500 500 400 300])

figure
plot(all_p, mean(Pb_p, 2), '*-')
ylim([0, 1])
xlabel('p')
ylabel('P_b')
set(gcf, 'Position', [500 500 400 300])

figure
semilogx(all_K, mean(Pb_K, 2), '*-')
ylim([0, 1])
xlabel('K')
ylabel('P_b')
set(gcf, 'Position', [500 500 400 300])
