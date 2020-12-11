load Data_phase

%% delete error data
Data_n = zeros(1000, 11);
n = 1;
for k = 1:size(Data, 1)
    if Data(k,9) < 1.1 && Data(k,9) > 0.9 &&...
       Data(k,6) > 300
        Data_n(n,:) = Data(k,:);
        n = n+1;
    end
end
Data = Data_n(1:n-1,:);

%% fit Kc and b
Kc = Data(:,7);
b = Data(:,11);
logb = log(b);

% coef = polyfit(kc, logb, 1);
% 
% x = min(kc): max(kc);
% y = exp(polyval(coef, x));

figure()
semilogy(1./Kc, b, '.', x, y, '-')

%% fit N2 and Kc
n2 = (30:2:50)';
Kc = zeros(length(n2), 1);
n2n = zeros(length(n2), 1);
for jn = 1:length(n2)
    N = [7, n2(jn), 5, 1, 5];
    Kc(jn) = get_Kc(N);
    n2n(jn) = n2(jn) / (sum(N) + n2(jn));
    
    disp(Kc(jn))
end

figure(1)
hold on
plot(n2n, Kc, 'b*')
coef = polyfit(n2n-0.5, 1./Kc, 2);
Kcf = 1./polyval(coef, n2n-0.5);
plot(n2n, Kcf, 'b-')
set(gcf,'position',[300,300,400,300])

disp(N)
disp(coef)
disp(sum((Kc-Kcf).^2) / sum((Kc-mean(Kc)).^2))

%% fit N4 and Ac
n4 = (1:20)';
n4n = zeros(length(n4), 1);
Ac = zeros(length(n4), 1);
for jn = 1:length(n4)
    N = [50, 200, 40, n4(jn), 30];
    Kc = get_Kc(N);
    n2n = N(2) / (sum(N) + N(2));
    Ac(jn) = 1 / (Kc * (n2n-0.5)^2);
    n4n(jn) = N(4) / (N(3) + N(5));
    
    disp(Kc)
    disp(Ac(jn))
end

figure(1)
hold on
plot(n4n, Ac, 'b*')
coef = polyfit(log(n4n), log(Ac), 1);
Acf = exp(coef(2)) * n4n.^coef(1);
plot(n4n, Acf, 'b-')
set(gcf,'position',[300,300,400,300])

disp(N)
disp([exp(coef(2)), coef(1)])
disp(sum((Ac-Acf).^2) / sum((Ac-mean(Ac)).^2))

figure(2)
hold on
plot(n4n, Ac, 'b*')
coef = sum(Ac.*n4n) / sum(n4n.*n4n);
Acf = coef * n4n;
plot(n4n, Acf, 'b-')
set(gcf,'position',[300,300,400,300])

disp(coef)
disp(sum((Ac-Acf).^2) / sum((Ac-mean(Ac)).^2))

%% fit N3 and Bc
n3 = (21:49)';
Ac = zeros(length(n3), 1);
Bc = zeros(length(n3), 1);
n3n = zeros(length(n3), 1);
n5n = zeros(length(n3), 1);
n1n = zeros(length(n3), 1);
for jn = 1:length(n3)
    N = [50, 200, n3(jn), 10, 70-n3(jn)];
    Kc = get_Kc(N);
    n2n = N(2) / (sum(N) + N(2));
    Ac(jn) = 1 / (Kc * (n2n-0.5)^2);
    Bc(jn) = Ac(jn) / N(4) * (N(3) + N(5));
    
    n3n(jn) = N(3) / (N(5) + N(3));
    n5n(jn) = N(5) / (N(5) + N(3));
    n1n(jn) = N(1) / (N(5) + N(3));
    
%     disp(Kc)
%     disp(Ac(jn))
    disp(Bc(jn))
end

figure(1)
hold on
nn = (n3n-n1n).*(n5n-n1n);
plot(nn, Ac, 'b*')
coef = polyfit(log(nn), log(Ac), 1);
Bcf = exp(coef(2)) * nn.^coef(1);
plot(nn, Bcf, 'b-')
set(gcf,'position',[300,300,400,300])

disp(N)
disp([exp(coef(2)), coef(1)])
% disp(sum((Bc-Bcf).^2) / sum((Bc-mean(Bc)).^2))

figure(2)
hold on
nn = (n3n-n1n).*(n5n-n1n);
plot(nn, Ac, 'b*')
coef = polyfit(nn, Ac, 1);
Bcf = polyval(coef, nn);
plot(nn, Bcf, 'b-')
set(gcf,'position',[300,300,400,300])

disp(coef)
% disp(sum((Bc-Bcf).^2) / sum((Bc-mean(Bc)).^2))

%% fit N1 and Ac
n1 = (42:75)';
Ac = zeros(length(n1), 1);
n1n = zeros(length(n1), 1);
for jn = 1:length(n1)
    N = [n1(jn), 200, 40, 5, 40];
    Kc = get_Kc(N);
    n2n = N(2) / (sum(N) + N(2));
    Ac(jn) = 1 / (Kc * (n2n-0.5)^2);
    n1n(jn) = N(1) - max(N(3), N(5));
    
    disp(Kc)
    disp(Ac(jn))
end

figure(1)
hold on
plot(n1n, Ac, 'b*')
coef = polyfit(log(n1n), log(Ac), 1);
Acf = exp(coef(2)) * n1n.^coef(1);
plot(n1n, Acf, 'b-')
set(gcf,'position',[300,300,400,300])

disp(N)
disp([exp(coef(2)), coef(1)])

figure(2)
hold on
plot(n1n, Ac, 'b*')
coef = polyfit(n1n, log(Ac), 1);
Acf = exp(polyval(coef, n1n));
plot(n1n, Acf, 'b-')
set(gcf,'position',[300,300,400,300])

disp(coef)