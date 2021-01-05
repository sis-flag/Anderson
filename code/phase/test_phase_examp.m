clear

L = [5, 24, 3, 1];
all_K = 680:20:780;

Kc = getCritical(L);
h = [L(2)/2; L(1); L(2); L(3); L(4); L(3); L(2)/2];
h = h / sum(h); % normalize

all_F0 = zeros(1, length(all_K));
all_F1 = zeros(1, length(all_K));
all_F2 = zeros(1, length(all_K));
all_lam = zeros(2, length(all_K));
labelK = [];
for jk = 1:length(all_K)
    K = all_K(jk);
    
    [U, lam] = eigPhase(K, h, 2);
    W = solvePhase(K, h);
    
    all_lam(:,jk) = lam;
    
    [u1, ~, hp1, hp2] = getvalPhase(U(:,1), h);
    all_F1(jk) = hp1 / (hp1 + hp2);
    
    [u2, ~, hp1, hp2] = getvalPhase(U(:,2), h);
    all_F2(jk) = hp1 / (hp1 + hp2);
    
    [w, x, hp1, hp2] = getvalPhase(W, h);
    all_F0(jk) = hp1 / (hp1 + hp2);
    
    figure(1)
    hold on
    plot(x, my_nmlz(u1))
    xlim([0, 1])
    set(gcf, 'Position', [500 500 400 300])
    
    figure(2)
    hold on
    u2 = my_nmlz(u2);
    plot(x, u2*sign(u2(75)))
    xlim([0, 1])
    set(gcf, 'Position', [500 500 400 300])
    
    figure(3)
    hold on
    plot(x, w)
    xlim([0, 1])
    set(gcf, 'Position', [500 500 400 300])

    labelK = [labelK; 'K = ', int2str(K)];
end

figure(1)
legend(labelK)

figure(2)
legend(labelK)

figure(3)
legend(labelK)

figure(4)
hold on
plot(all_K, all_lam, '*-')
legend('\lambda_1', '\lambda_2')
xlabel('K')
ylabel('\lambda')
set(gcf, 'Position', [500 500 400 300])

figure(5)
plot(all_K, all_F0, '*-')
ylim([0.5, 0.55])
xlabel('K')
ylabel('F')
set(gcf, 'Position', [500 500 400 300])

figure(6)
plot(all_K, all_F1, '*-')
ylim([0, 1])
xlabel('K')
ylabel('F')
set(gcf, 'Position', [500 500 400 300])

figure(7)
plot(all_K, all_F2, '*-')
ylim([0, 1])
xlabel('K')
ylabel('F')
set(gcf, 'Position', [500 500 400 300])