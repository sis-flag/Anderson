clear 

%% simulate and predict
Data = zeros(1000, 10);
count = 1;

for N3 = linspace(3, 6, 5)
    for N1 = linspace(N3+0.5, 1.5*N3, 5)
        for N4 = linspace(0.5, 0.5*N3, 5)
            for N2 = 2*ceil(N1+N3+N4+N3):2:3*floor(N1+N3+N4+N3)

N = [N1, N2, N3, N4, N3];

[Kc_s, b, lamc_s] = getP_Kc(N);

Df = @(x) Dfunc1(x(1), x(2), N);
options = optimoptions('fsolve','Display','none');
res = fsolve(Df, [sqrt(lamc_s); sqrt(Kc_s-lamc_s)], options);
Kc_p = res(1)^2 + res(2)^2;
lamc_p = res(1)^2;

fprintf('%d & %g & %g & %g & %g & %g & %g & %g & %g & %g \\\\\n',...
    count, N(1), N(2), N(3), N(4),...
    Kc_s, lamc_s, Kc_p, lamc_p, b);

Data(count,:) = [N, Kc_s, lamc_s, Kc_p, lamc_p, b];
count = count + 1;
            end
        end
    end
end

Data = Data(1:count-1,:);
save('DataP', 'Data')

%% fit
clear;
load DataP;

figure
hold on
plot(Data(:,8), Data(:,6), '.')
plot([0,max(Data(:,8))], [0,max(Data(:,8))], '-')

figure
hold on
plot(Data(:,9), Data(:,7), '.')
plot([0,max(Data(:,7))], [0,max(Data(:,7))], '-')

disp('mse of Kc:')
disp(mean((Data(:,8)-Data(:,6)).^2))
disp('mse of lam:')
disp(mean((Data(:,9)-Data(:,7)).^2))

%% fit
clear;
load DataP;

L2 = Data(:,2) ./ sum(Data(:,[1,2,2,3,4,5]), 2);
L4 = Data(:,4) ./ sum(Data(:,[3,5]), 2);
L1 = Data(:,1) ./ sum(Data(:,[3,5]), 2);
Kc = Data(:,6);

Fit = (1-2*L2).^2 .* L4 .* exp(7.1 * (L1-0.5) - 4);

figure
hold on
plot(Fit, 1./Kc, '.')
plot([0, 0.0012], [0,0.0012], '-')

% X = [log(1-2*L2), log(L4), log(L1-0.5), ones(size(L1))];
% Y = log(Kc);
% coef = X \ Y;
% disp(coef)
% mse = mean((X*coef - Y).^2);
