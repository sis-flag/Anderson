clear

L = [5, 24, 3, 1];
h = [L(2)/2; L(1); L(2); L(3); L(4); L(3); L(2)/2];
h = h / sum(h);

h1 = [L(2)/2; L(1); L(2)+L(3)+L(4)+L(3)+L(2)/2];
h1 = h1 / sum(h1);

h2 = [L(2)/2+L(1)+L(2); L(3); L(4); L(3); L(2)/2];
h2 = h2 / sum(h2);

%%
K = 800;

[U0, lam0] = eigPhase(K, h, 2);
[U1, lam1] = eigPhase(K, h1, 1);
[U2, lam2] = eigPhase(K, h2, 1);

[u01, x01] = getvalPhase(U0(:,1), h);
[u02, x02] = getvalPhase(U0(:,2), h);
[u1, x1] = getvalPhase(U1(:,1), h1);
[u2, x2] = getvalPhase(U2(:,1), h2);

px = linspace(0, 1, 201);
pu01 = interp1(x01, my_nmlz(u01), px);
pu02 = interp1(x02, my_nmlz(u02), px);
pu1 = interp1(x1, my_nmlz(u1), px);
pu2 = interp1(x2, my_nmlz(u2), px);

figure
hold on
plot(px, pu1, '.', 'MarkerSize', 14)
plot(px, pu2, '.', 'MarkerSize', 14)
plot(px, pu01, '-', 'LineWidth', 1)
plot(px, pu02, '-', 'LineWidth', 1)

legend([ 'u_1 of S_1'; 'u_1 of S_2'; 'u_1 of S  '; 'u_2 of S  '],...
    'Location', 'NorthEastOut')
xlim([0, 1])
set(gcf, 'Position', [300 300 350 200])
set(gca, 'FontSize', 14)

disp(lam0)
disp(lam1)
disp(lam2)

%%
all_K = 680:10:780;

lam0 = zeros(2, length(all_K));
lam1 = zeros(1, length(all_K));
lam2 = zeros(1, length(all_K));
for jk = 1:length(all_K)
    K = all_K(jk);
    
    [~, lam0(:,jk)] = eigPhase(K, h, 2);
    [~, lam1(jk)] = eigPhase(K, h1, 1);
    [~, lam2(jk)] = eigPhase(K, h2, 1);
end

figure
hold on
plot(all_K, lam1, '-', 'LineWidth', 1)
plot(all_K, lam2, '-', 'LineWidth', 1)
plot(all_K, lam0(1,:), 'o', 'MarkerSize', 6)
plot(all_K, lam0(2,:), '*', 'MarkerSize', 6)
xlabel('K')
ylabel('\lambda')
xlim([670, 790])
ylim([340, 380])
legend([ '\lambda_1 of S_1'; '\lambda_1 of S_2'; '\lambda_1 of S  '; '\lambda_2 of S  '],...
    'Location', 'NorthEastOut')
set(gcf, 'Position', [300 300 350 200])
set(gca, 'FontSize', 14)