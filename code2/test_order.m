%%
clear

%% fit P1 and Kc
tn = linspace(20, 40, 10)';
Kc = zeros(length(tn), 1);
P1 = zeros(length(tn), 1);
P2 = zeros(length(tn), 1);
P3 = zeros(length(tn), 1);
for jn = 1:length(tn)
    L = [7, tn(jn), 5, 1];
    L = L / (sum(L) + L(2) + L(3));
    P1(jn) = 1 - 2*L(2);
    P2(jn) = L(1) / (L(1) + 2*L(3) + L(4));
    P3(jn) = L(4) / (2*L(3) + L(4));
    Kc(jn) = getCritical(L);
end

xx = P1.^2;
yy = 1./Kc;
coef = sum(yy.*xx) / sum(xx.*xx);
yyf = coef * xx;

figure(1)
hold on
plot(xx, yy, '+', 'MarkerSize', 7)
plot(xx, yyf, '-', 'LineWidth', 1)
xlabel('{P_1}^2')
ylabel('1/K_c')
xlim([0.03, 0.1])
ylim([2e-4, 8e-4])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 16)

%% fit P2 and Kc
tn = linspace(5.5, 8.5, 10)';
Kc = zeros(length(tn), 1);
P1 = zeros(length(tn), 1);
P2 = zeros(length(tn), 1);
P3 = zeros(length(tn), 1);
for jn = 1:length(tn)
    L = [tn(jn), 2*(11+tn(jn)), 5, 1];
    L = L / (sum(L) + L(2) + L(3));
    P1(jn) = 1 - 2*L(2);
    P2(jn) = L(1) / (L(1) + 2*L(3) + L(4));
    P3(jn) = L(4) / (2*L(3) + L(4));
    Kc(jn) = getCritical(L);
end

xx = P2;
yy = log(Kc);
coef = polyfit(xx, yy, 1);
yyf = polyval(coef, xx);

figure(2)
hold on
plot(xx, yy, '+', 'MarkerSize', 7)
plot(xx, yyf, '-', 'LineWidth', 1)
xlabel('P_2')
ylabel('log(K_c)')
xlim([0.33, 0.44])
ylim([6.5, 9.5])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 16)

%% fit P3 and Ac
tn = linspace(0.4, 1.6, 10)';
Kc = zeros(length(tn), 1);
P1 = zeros(length(tn), 1);
P2 = zeros(length(tn), 1);
P3 = zeros(length(tn), 1);
for jn = 1:length(tn)
    L = [7, 30, 5-tn(jn)/2, tn(jn)];
    L = L / (sum(L) + L(2) + L(3));
    P1(jn) = 1 - 2*L(2);
    P2(jn) = L(1) / (L(1) + 2*L(3) + L(4));
    P3(jn) = L(4) / (2*L(3) + L(4));
    Kc(jn) = getCritical(L);
end

xx = P3.^2;
yy = 1./Kc;
coef = sum(yy.*xx) / sum(xx.*xx);
yyf = coef * xx;

figure(3)
hold on
plot(xx, yy, '+', 'MarkerSize', 7)
plot(xx, yyf, '-', 'LineWidth', 1)
xlabel('{P_3}^2')
ylabel('1/K_c')
xlim([0, 0.03])
ylim([0, 2e-3])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 16)
