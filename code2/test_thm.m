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
Kc0 = Data(:, 8);
lam0 = Data(:, 7);
Kc1 = Data(:, 6);
lam1 = Data(:, 5);

disp('RMSE of Kc:')
disp(100 * mean(abs(Kc0-Kc1) ./ Kc1))
disp('RMSE of lam:')
disp(100 * mean(abs(lam0-lam1) ./ lam1))
