clear
load('future.mat')

%% Bernoulli
K = 1e4; p = 1/2;

figure
hold on
plot(all_h, Pb_hp(:,1), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,1), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,1), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
hold on
plot([0,0], [0,1-p], 'k-', 'LineWidth', 1)
plot([K,K], [0,p], 'k-', 'LineWidth', 1)
xticks([0, K])
xticklabels({'0','10^4'})
yticks(p)
yticklabels({'0.5'})
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 1])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Normal
mu = 5e3; sigma = 5e3;

figure
hold on
plot(all_h, Pb_hp(:,2), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,2), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,2), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
x = linspace(-5e3, 15e3, 100);
y = pdf('Normal', x, mu, sigma);
plot(x, y, 'k-', 'LineWidth', 1)
xticks(5e3)
xticklabels({'5\times10^3'})
yticks([])
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 2.6e-4])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Gamma
a = 1; b = 5e3;

figure
hold on
plot(all_h, Pb_hp(:,3), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,3), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,3), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
x = linspace(-5e3, 15e3, 100);
y = pdf('Gamma', x, a, b);
plot(x, y, 'k-', 'LineWidth', 1)
xticks(5e3)
xticklabels({'5\times10^3'})
yticks([])
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 2.6e-4])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Uniform

a = 5e3 - sqrt(3)*5e3; b = 5e3 + sqrt(3)*5e3;

figure
hold on
plot(all_h, Pb_hp(:,4), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,4), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,4), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
plot([a, b], [1,1]/(b-a), 'k-', 'LineWidth', 1)
xticks([a, b])
xticklabels({'-3.67\times10^3','13.67\times10^3'})
yticks(1/(b-a))
yticklabels({'p'})
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 2.6e-4])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Bernoulli
K = 2e4/3; p = 3/4;

figure
hold on
plot(all_h, Pb_hp(:,5), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,5), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,5), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
hold on
plot([0,0], [0,1-p], 'k-', 'LineWidth', 1)
plot([K,K], [0,p], 'k-', 'LineWidth', 1)
xticks([0, K])
xticklabels({'0','2/3\times10^4'})
yticks([1-p, p])
yticklabels({'0.25','0.75'})
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 1])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Normal
mu = 5e3; sigma = 5e3 / sqrt(3);

figure
hold on
plot(all_h, Pb_hp(:,6), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,6), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,6), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
x = linspace(-5e3, 15e3, 100);
y = pdf('Normal', x, mu, sigma);
plot(x, y, 'k-', 'LineWidth', 1)
xticks(5e3)
xticklabels({'5\times10^3'})
yticks([])
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 2.6e-4])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Gamma
a = 3; b = 5e3/3;

figure
hold on
plot(all_h, Pb_hp(:,7), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,7), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,7), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
x = linspace(-5e3, 15e3, 100);
y = pdf('Gamma', x, a, b);
plot(x, y, 'k-', 'LineWidth', 1)
xticks(5e3)
xticklabels({'5\times10^3'})
yticks([])
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 2.6e-4])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Uniform
a = 0; b = 1e4;

figure
hold on
plot(all_h, Pb_hp(:,8), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,8), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,8), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
plot([a, b], [1,1]/(b-a), 'k-', 'LineWidth', 1)
xticks([a, b])
xticklabels({'0','10^4'})
yticks(1/(b-a))
yticklabels({'10^{-4}'})
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 2.6e-4])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Bernoulli
K = 5e4/9; p = 0.9;

figure
hold on
plot(all_h, Pb_hp(:,9), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,9), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,9), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
hold on
plot([0,0], [0,1-p], 'k-', 'LineWidth', 1)
plot([K,K], [0,p], 'k-', 'LineWidth', 1)
xticks([0, K])
xticklabels({'0','5/9\times10^4'})
yticks([1-p, p])
yticklabels({'0.1', '0.9'})
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 1])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Normal
mu = 5e3; sigma = 5e3/3;

figure
hold on
plot(all_h, Pb_hp(:,10), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,10), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,10), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
x = linspace(-5e3, 15e3, 100);
y = pdf('Normal', x, mu, sigma);
plot(x, y, 'k-', 'LineWidth', 1)
xticks(5e3)
xticklabels({'5\times10^3'})
yticks([])
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 2.6e-4])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Gamma
a = 9; b = 5e3/9;

figure
hold on
plot(all_h, Pb_hp(:,11), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,11), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,11), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
x = linspace(-5e3, 15e3, 100);
y = pdf('Gamma', x, a, b);
plot(x, y, 'k-', 'LineWidth', 1)
xticks(5e3)
xticklabels({'5\times10^3'})
yticks([])
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 2.6e-4])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)

%% Uniform
a = (1-1/sqrt(3))*5e3; b = (1+1/sqrt(3))*5e3;

figure
hold on
plot(all_h, Pb_hp(:,12), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pe_hp(:,12), '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h, Pc_hp(:,12), '.-', 'MarkerSize', 15, 'LineWidth', 1)
ylim([0, 1])
xlabel('h')
xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
xlim([1e-1, 1e3])
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')

axes('Position',[0.65, 0.65, 0.3, 0.3]);
plot([a, b], [1,1]/(b-a), 'k-', 'LineWidth', 1)
xticks([a, b])
xticklabels({'2.1\times10^3','7.9\times10^3'})
yticks(1/(b-a))
yticklabels({'p'})
box on
grid on
xlim([-5e3, 15e3])
ylim([0, 2.6e-4])
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 10)
