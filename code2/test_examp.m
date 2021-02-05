clear

L = [5, 24, 3, 1];
h = [L(2)/2; L(1); L(2); L(3); L(4); L(3); L(2)/2];
h = h / sum(h); % normalize

all_K = 700:25:750;

for jk = 1:length(all_K)
    K = all_K(jk);
    [U, lam] = eigPhase(K, h, 2);
    
    [u1, ~, hp11, hp12] = getvalPhase(U(:,1), h);
    [u2, x, hp21, hp22] = getvalPhase(U(:,2), h);
    
    figure
    hold on
    plot(x, u1 / (hp11 + hp12), 'LineWidth', 1)
    plot(x, u2 / (hp21 + hp22), 'LineWidth', 1)
    legend(['u_1'; 'u_2'], 'Location', 'NorthEastOut')
    xlim([0, 1])
    title(['K = ', int2str(K)])
    set(gcf, 'Position', [300 300 350 200])
    set(gca, 'FontSize', 14)
end
