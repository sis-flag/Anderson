clear

L = [5, 24, 3, 1];
h = [L(2)/2; L(1); L(2); L(3); L(4); L(3); L(2)/2];
h = h / sum(h); % normalize

all_K = [700, 723.5, 750];

for jk = 1:length(all_K)
    K = all_K(jk);
    [U, lam] = eigPhase(K, h, 2);
    
    [u1, x, hp11, hp12] = getvalPhase(U(:,1), h);
    
    figure
    hold on
    cl = bar([0; cumsum(h)], [1,0,1,0,1,0,1,0],'histc');
    cl.LineStyle = 'None';
    cl.FaceColor = [0.8, 0.8, 0.8];
    plot(x, u1 / (hp11 + hp12), 'LineWidth', 1)
    xlim([0, 1])
    title(['K = ', num2str(K)])
    set(gcf, 'Position', [300 300 350 300])
    set(gca, 'FontSize', 16)
end
