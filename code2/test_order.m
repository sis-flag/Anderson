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
P1 = 10.^linspace(-0.7, -0.5, 30)';
P2 = 0.4;
P3 = 0.1;
Kc = zeros(length(P1), 1);
for jn = 1:length(P1)
    Kc(jn) = getCriticalP([P1(jn); P2; P3]);
end

xx = log10(P1);
yy = log10(Kc);
coef = polyfit(xx, yy, 1);
yyf = polyval(coef, xx);
Kcf = 10.^yyf;
R2 = sum((Kc-Kcf).^2) / sum((Kc-mean(Kc)).^2);

disp('P1-Kc')
disp(coef)
disp('R2')
disp(R2)

%% plot
figure
hold on
plot(xx, yy, 'bo', 'MarkerSize', 4)
plot(xx, yyf, 'r-', 'LineWidth', 1)
xlabel('log_{10}(P_1)')
ylabel('log_{10}(K_c)')
legend('Slop = -2','R^2 = 1 - 1.2 \times 10^{-8}')
set(gcf, 'Position', [300 300 400 400])
set(gca, 'FontSize', 16)

%% fit P2 and Kc
P1 = 0.25;
P2 = linspace(0.38, 0.42, 30)';
P3 = 0.1;
Kc = zeros(length(P2), 1);
for jn = 1:length(P2)
    Kc(jn) = getCriticalP([P1; P2(jn); P3]);
end

xx = P2;
yy = log10(Kc);
coef = polyfit(xx, yy, 1);
yyf = polyval(coef, xx);
Kcf = 10.^yyf;
R2 = sum((Kc-Kcf).^2) / sum((Kc-mean(Kc)).^2);

disp('P2-Kc')
disp(coef)
disp('R2')
disp(R2)

%% plot
figure
hold on
plot(xx, yy, 'bo', 'MarkerSize', 4)
plot(xx, yyf, 'r-', 'LineWidth', 1)
xlabel('P_2')
ylabel('log_{10}(K_c)')
legend('Slop = -23.2','R^2 = 1 - 1.1 \times 10^{-3}')
set(gcf, 'Position', [300 300 400 400])
set(gca, 'FontSize', 16)

%% fit P3 and Kc
P1 = 0.25;
P2 = 0.4;
P3 = 10.^linspace(-1.1, -0.9, 30)';
Kc = zeros(length(P3), 1);
for jn = 1:length(P3)
    Kc(jn) = getCriticalP([P1; P2; P3(jn)]);
end

xx = log10(P3);
yy = log10(Kc);
coef = polyfit(xx, yy, 1);
yyf = polyval(coef, xx);
Kcf = 10.^yyf;
R2 = sum((Kc-Kcf).^2) / sum((Kc-mean(Kc)).^2);

disp('P3-Kc')
disp(coef)
disp('R2')
disp(R2)

%% plot
figure
hold on
plot(xx, yy, 'bo', 'MarkerSize', 4)
plot(xx, yyf, 'r-', 'LineWidth', 1)
xlabel('log_{10}(P_3)')
ylabel('log_{10}(K_c)')
legend('Slop = -1.7','R^2 = 1 - 6.8 \times 10^{-4}')
set(gcf, 'Position', [300 300 400 400])
set(gca, 'FontSize', 16)
