clear
rng(0)

%% simulate and predict
Nsamp = 200;

Data = zeros(Nsamp, 8);
for count = 1:Nsamp
    L3 = rand() * 0.025 + 0.03;
    L1 = (rand() * 0.4 + 1.3) * L3;
    L4 = (rand() * 0.2 + 0.2) * L3;
    L2 = (1 - L1 - 2*L3 - L4) / 2;

    L = [L1, L2, L3, L4];
    L = L / (sum(L) + L(2) + L(3));

    [Kc_s, lamc_s] = getCritical(L);

    Df = @(x) Dfunc(x(1), x(2), L);
    options = optimoptions('fsolve','Display','none');
    res = fsolve(Df, [sqrt(lamc_s); sqrt(Kc_s-lamc_s)], options);
    Kc_p = res(1)^2 + res(2)^2;
    lamc_p = res(1)^2;

    fprintf('%d & %g & %g & %g & %g & %g & %g & %g & %g \\\\\n',...
        count, L(1), L(2), L(3), L(4),...
        Kc_s, lamc_s, Kc_p, lamc_p);

    Data(count,:) = [L, Kc_s, lamc_s, Kc_p, lamc_p];
end

save('Phase.mat', 'Data')

%% fit
clear
load('Phase.mat')

% figure
% hold on
% plot(Data(1:7:end,6), Data(1:7:end,8), '.', 'MarkerSize', 8)
% plot([0,max(Data(:,8))], [0,max(Data(:,8))], '-', 'LineWidth', 1)
% xlabel('prediction')
% ylabel('simulation')
% title('\lambda_c')
% set(gcf, 'Position', [300 300 400 300])
% set(gca, 'FontSize', 16)

figure
hold on
plot(Data(:,5)*10, Data(:,7)*10, 'o', 'MarkerSize', 5)
plot([0,max(Data(:,7))]*10, [0,max(Data(:,7))]*10, '-', 'LineWidth', 1)
xlabel('prediction')
ylabel('simulation')
title('K_c')
xlim([0, 5e4])
ylim([0, 5e4])
xticks(0:1e4:5e4)
yticks(0:1e4:5e4)
axis square
set(gcf, 'Position', [300 300 300 300])
set(gca, 'FontSize', 17.6)

disp('mse of Kc:')
disp(mean((Data(:,8)-Data(:,6)).^2))
disp('mse of lam:')
disp(mean((Data(:,5)-Data(:,7)).^2))
