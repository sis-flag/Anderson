%% fit L2 and Kc
tn = linspace(20, 40, 10)';
Kc = zeros(length(tn), 1);
L2 = zeros(length(tn), 1);
for jn = 1:length(tn)
    N = [7, tn(jn), 5, 1, 5];
    L2(jn) = N(2) / (sum(N) + N(2));
    Kc(jn) = getP_Kc(N);
    fprintf('%g\n', Kc(jn));
end

xx = (1-2*L2).^2;
yy = 1./Kc;
coef = sum(yy.*xx) / sum(xx.*xx);
yyf = coef * xx;

figure(1)
hold on
plot(xx, yy, 'b*')
plot(xx, yyf, 'b--')

disp(coef)
disp(sum((yy-yyf).^2) / sum((yy-mean(yy)).^2))

%% fit L4 and Ac
tn = linspace(0.4, 1.6, 10)';
Ac = zeros(length(tn), 1);
L4 = zeros(length(tn), 1);
for jn = 1:length(tn)
    N = [5, 2*(14+tn(jn)), 3, tn(jn), 3];
    L2 = N(2) / (sum(N) + N(2));
    L4(jn) = N(4) / sum(N([3,5]));
    
    Kc = getP_Kc(N);
    Ac(jn) = 1./(Kc * (1-2*L2)^2);
    fprintf('%g\t%g\n', Kc, Ac(jn));
end

xx = L4;
yy = Ac;
coef = sum(yy.*xx) / sum(xx.*xx);
yyf = coef * xx;

figure(1)
hold on
plot(xx, yy, 'm*')
plot(xx, yyf, 'm--')

disp(coef)
disp(sum((yy-yyf).^2) / sum((yy-mean(yy)).^2))

%% fit L1 and Bc
tn = linspace(5.5, 7.5, 10)';
Bc = zeros(length(tn), 1);
L1 = zeros(length(tn), 1);
for jn = 1:length(tn)
    N = [tn(jn), 2*(9+tn(jn)), 5, 1, 5];
    L2 = N(2) / (sum(N) + N(2));
    L4 = N(4) / sum(N([3,5]));
    L1(jn) = N(1) / sum(N([3,5]));
    
    Kc = getP_Kc(N);
    Ac = 1./(Kc * (1-2*L2)^2);
    Bc(jn) = Ac ./ L4;
    fprintf('%g\t%g\t%g\n', Kc, Ac, Bc(jn));
end

xx = L1 - 0.5;
yy = log(Bc);
coef = polyfit(xx, yy, 1);
yyf = polyval(coef, xx);

figure(1)
hold on
plot(xx, yy, 'b*')
plot(xx, yyf, 'b--')

disp(coef)
disp(sum((yy-yyf).^2) / sum((yy-mean(yy)).^2))

% xx = log(1 - L1);
% yy = log(Bc);
% coef = polyfit(xx, yy, 1);
% yyf = polyval(coef, xx);
% 
% figure(2)
% hold on
% plot(xx, yy, 'r*')
% plot(xx, yyf, 'r--')
% 
% disp(coef)
% disp(sum((yy-yyf).^2) / sum((yy-mean(yy)).^2))
