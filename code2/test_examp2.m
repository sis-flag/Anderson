clear

L = [5, 24, 3, 1];
h = [L(2)/2; L(1); L(2); L(3); L(4); L(3); L(2)/2];
h = h / sum(h); % normalize

all_K = 700:1:750;

all_F1 = zeros(length(all_K), 1);
all_F2 = zeros(length(all_K), 1);
for jk = 1:length(all_K)
    K = all_K(jk);
    [U, lam] = eigPhase(K, h, 2);
    
    [~, ~, hp11, hp12] = getvalPhase(U(:,1), h);
    [~, ~, hp21, hp22] = getvalPhase(U(:,2), h);
    
    all_F1(jk) = hp11 / (hp11 + hp12);
    all_F2(jk) = hp21 / (hp21 + hp22);
end

figure
hold on
plot(all_K, all_F1, '-', 'LineWidth', 1)
plot([700, 750], [0.5, 0.5], 'k--', 'LineWidth', 0.5)
plot(723.5, 0.5, 'ko')
text(723.5, 0.5, 'crtitcal point')
xlabel('K')
ylabel('F')
xlim([700, 750])
set(gcf, 'Position', [300 300 350 300])
set(gca, 'FontSize', 16)