%% fit l2 and Kc
tn = (30:2:50)';
Kc = zeros(length(tn), 1);
L2 = zeros(length(tn), 1);
for jn = 1:length(tn)
    N = [6, tn(jn), 4, 1, 4];
    L = N / (sum(N) + N(2));
    L2(jn) = L(2);
    
    Kc(jn) = get_Kc(N);
    
    disp(Kc(jn))
end

xx = (1-2*L2).^2;
yy = 1./Kc;
coef = sum(yy.*xx) / sum(xx.*xx);
yyf = coef * xx;

figure(1)
hold on
plot(xx, yy, 'r*')
plot(xx, yyf, 'r--')

disp(N)
disp(coef)
disp(sum((yy-yyf).^2) / sum((yy-mean(yy)).^2))

%% fit L4 and Ac
tn = (1:12)';
Ac = zeros(length(tn), 1);
L4 = zeros(length(tn), 1);
for jn = 1:length(tn)
    N = [130, 400, 100-tn(jn), 2*tn(jn), 100-tn(jn)];
    L = N / (sum(N) + N(2));
    L4(jn) = L(4) / sum(L([1,3,4,5]));
    
    Kc = get_Kc(N);
    Ac(jn) = 1 / (Kc * (1-2*L(2))^2);
    
    disp(Kc)
end

xx = L4;
yy = Ac;
coef = polyfit(xx, yy, 2);
yyf = polyval(coef, xx);

figure(1)
hold on
plot(xx, yy, 'b*')
plot(xx, yyf, 'b--')

disp(N)
disp(coef)
disp(sum((yy-yyf).^2) / sum((yy-mean(yy)).^2))

%% fit L1 and Ac
tn = (21:29)';
Ac = zeros(length(tn), 1);
L1 = zeros(length(tn), 1);
for jn = 1:length(tn)
    N = [2*tn(jn), 150, 60-tn(jn), 7, 60-tn(jn)];
    L = N / (sum(N) + N(2));
    L1(jn) = L(1) / sum(L([1,3,4,5]));
    
    Kc = get_Kc(N);
    Ac(jn) = 1 / (Kc * (1-2*L(2))^2);
    
    disp(Kc)
%     disp(Ac(jn))
end

xx = L1;
yy = Ac;
coef = polyfit(xx, yy, 2);
yyf = polyval(coef, xx);

figure(1)
hold on
plot(xx, yy, 'b*')
plot(xx, yyf, 'b--')

disp(N)
disp(coef)
disp(sum((yy-yyf).^2) / sum((yy-mean(yy)).^2))

% %% fit N3 and Bc
% tn = (21:80-21)';
% LL = zeros(length(tn), 1);
% L35 = zeros(length(tn), 1);
% Bc = zeros(length(tn), 1);
% for jn = 1:length(tn)
%     N = [60, 200, tn(jn), 10, 80-tn(jn)];
%     L = N / (sum(N) + N(2));
%     L2 = sum(L([1,3,4,5]));
%     L35(jn) = L(3) - L(5);
%     LL(jn) = (L(3)-L(1))*(L(5)-L(1));
%     
%     Kc = get_Kc(N);
%     Ac = 1 / (Kc * L2^2);
%     Bc(jn) = Ac / L(4);
%     
% %     disp(Kc)
% %     disp(Ac)
% %     disp(Bc(jn))
% end
% 
% figure(2)
% hold on
% plot(L35, Bc, 'r*')
% set(gcf,'position',[300,300,400,300])
% 
% figure(1)
% hold on
% plot(LL, Bc, 'r*')
% coef = polyfit(LL, Bc, 1);
% Bcf = polyval(coef, LL);
% plot(LL, Bcf, 'r--')
% set(gcf,'position',[300,300,400,300])
% 
% disp(N)
% disp(coef)
% disp(sum((Bc-Bcf).^2) / sum((Bc-mean(Bc)).^2))
