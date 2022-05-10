%% 1d example
clear;

% peremeters
K = 9000;
N = 30;
h = 0;

% rng(12);
% V = exprnd(K, N, 1);
rng(0);
V = K/3 * randn(N, 1) + K;

[U, lam] = eig1d(V, h, 4);
W = solve1d(V, h);

figure
hold on

yyaxis right
cl = bar(((1:N) - 0.5) / N, V);
cl.BarWidth = 1;
cl.LineStyle = 'None';
cl.FaceColor = [0.8, 0.8, 0.8];
ylabel('potential')

yyaxis left
w = getval1d(W);
x = linspace(0, 1, length(w));
plot(x, w, 'k-', 'LineWidth', 1)
for k = 1:4
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, abs(uk) / lam(k), '-', 'LineWidth', 1)
end
xlim([0, 1])
xlabel('x')
ylabel('u(x) or w(x)')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 14)

%% 2d example
clear;

% peremeters
K = 9000;
N = 20;
h = 0;

% rng(5);
% V = exprnd(K, N, N);
rng(12);
V = K/3 * randn(N, N)' + K;

[U, lam] = eig2d(V, h, 4);
W = solve2d(V, h);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x);

% plot eigen mode and valley line
figure
hold on
for k = 1:4
    uk = my_nmlz(getval2d(U(:,k)));
    surf(x1, x2, uk, 'LineStyle', 'None');
end
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 12)
zticks([])
view(-150, 70)
axis off

%% relation to h
clear
rng(0)
tic

N_samp = 200;
N1d = 50;
N2d = 15;

thes = 0.5;

p = 0.5;
K = 1e4;
all_h = [0, 1e-1, 3e-1, 1e0, 3e0, 1e1, 3e1, 1e2, 3e2, 1e3, inf];

Pb_h = zeros(length(all_h), N_samp, 'logical');
for jh = 1:length(all_h)
    h = all_h(jh);
    for js = 1:N_samp
%         V = exprnd(K, N1d, 1);
        V = K/3 * randn(N1d, 1) + K;
        U = eig1d(V, h, 1);
        u = abs(getval1d(U));
        tPb = max(u(1), u(end)) / max(u) > thes;
        Pb_h(jh, js) = tPb;
    end
    
    fprintf('1d h %d / %d finished ', jh, length(all_h));
    toc
end

Pe_h = zeros(length(all_h), N_samp, 'logical');
Pc_h = zeros(length(all_h), N_samp, 'logical');
for jh = 1:length(all_h)
    h = all_h(jh);
    for js = 1:N_samp
%         V = exprnd(K, N2d, N2d);
        V = K/3 * randn(N2d, N2d) + K;
        U = eig2d(V, h, 1);
        u = abs(getval2d(U));
        tPe = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(u(:)) > thes;
        tPc = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(u(:)) > thes;
        Pe_h(jh, js) = tPe;
        Pc_h(jh, js) = tPc;
    end
    
    fprintf('2d h %d / %d finished ', jh, length(all_h));
    toc
end

%% plot

figure
hold on
plot(all_h(2:end-1), mean(Pb_h(2:end-1,:), 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h(2:end-1), mean(Pe_h(2:end-1,:), 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
plot(all_h(2:end-1), mean(Pc_h(2:end-1,:), 2),...
    '.-', 'MarkerSize', 15, 'LineWidth', 1)
legend('P_b', 'P_e', 'P_c')
ylim([0, 1])
xlim([3e-2, 3e3])
xticks(all_h(2:2:end))
xlabel('h')
ylabel('probability')
set(gcf, 'Position', [300 300 400 300])
set(gca, 'FontSize', 15)
set(gca, 'XScale', 'log')
