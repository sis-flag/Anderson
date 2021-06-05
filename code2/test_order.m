%%
clear

%% convert L and P
P2L = @(P)...
    [P(1) * P(2);...
    (1-P(1)) / 2;...
    P(1) * (1-P(2)) * (1-P(3)) / 2;...
    P(1) * (1-P(2)) * P(3)];
L2P = @(L)...
    [1 - 2*L(2) / (sum(L) + L(2) + L(3));...
    L(1) / (L(1) + 2*L(3) + L(4));...
    L(4) / (2*L(3) + L(4))];

getCriticalP = @(P) getCritical(P2L(P));

%% fit P1 and Kc
P1 = linspace(0.20, 0.30, 12)';
P2 = 0.4;
P3 = 0.1;
Kc = zeros(length(P1), 1);
for jn = 1:length(P1)
    Kc(jn) = getCriticalP([P1(jn); P2; P3]);
end

xx = log(P1);
yy = log(Kc);
coef = polyfit(xx, yy, 1);
yyf = polyval(coef, xx);

disp('P1-Kc')
disp(coef)
disp('R2')
disp(sum((yy-yyf).^2) / sum((yy-mean(yy)).^2))

%% plot
figure
hold on
plot([-1.55, -1.35], [7.4, 7.0], 'k--', 'LineWidth', 1)
plot(xx, yy, 'bo', 'MarkerSize', 5)
plot(xx, yyf, 'r-', 'LineWidth', 1)
xlabel('log(P_1)')
ylabel('log(K_c)')
legend('Slop = -2')
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 16)

%% fit P2 and Kc
P1 = 0.25;
P2 = linspace(0.33, 0.44, 12)';
P3 = 0.1;
Kc = zeros(length(P2), 1);
for jn = 1:length(P2)
    Kc(jn) = getCriticalP([P1; P2(jn); P3]);
end

xx = P2;
yy = log(Kc);
coef = polyfit(xx, yy, 1);
yyf = polyval(coef, xx);

disp('P2-Kc')
disp(coef)
disp('R2')
disp(sum((yy-yyf).^2) / sum((yy-mean(yy)).^2))

%% plot
figure
hold on
plot(xx, yy, 'bo', 'MarkerSize', 5)
plot(xx, yyf, 'r-', 'LineWidth', 1)
xlabel('P_2')
ylabel('log(K_c)')
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 16)

%% fit P3 and Kc
P1 = 0.25;
P2 = 0.4;
P3 = linspace(0.05, 0.15, 12)';
Kc = zeros(length(P3), 1);
for jn = 1:length(P3)
    Kc(jn) = getCriticalP([P1; P2; P3(jn)]);
end

xx = log(P3);
yy = log(Kc);
coef = polyfit(xx, yy, 1);
yyf = polyval(coef, xx);

disp('P3-Kc')
disp(coef)
disp('R2')
disp(sum((yy-yyf).^2) / sum((yy-mean(yy)).^2))

%% plot
figure(3)
hold on
plot([-2.8, -2.2], [7.5, 6.48], 'k--', 'LineWidth', 1)
plot(xx, yy, 'bo', 'MarkerSize', 5)
plot(xx, yyf, 'r-', 'LineWidth', 1)
xlabel('log(P_3)')
ylabel('log(K_c)')
ylim([6.3, 8.7])
legend('Slope = -1.7')
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 16)
