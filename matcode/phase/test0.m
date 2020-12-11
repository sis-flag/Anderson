clear

%%
K = 400;

L = [3, 10, 2, 1, 2];
L = L(:) / (sum(L) + L(2));
x0 = L(1) / 2;
x1 = L(4) / 2;
x2 = x1 + L(5);
x3 = 0.5;

Df1 = @(a,b) a*tan(a*x0) - b*tanh(b*(x3-x0));
Df2 = @(a,b) (a^2 - b^2) * (exp(2*b*x2) + exp(2*b*(x1+x3))) / (-exp(2*b*x2) + exp(2*b*(x1+x3)))...
           + (a^2 + b^2) * (exp(2*b*x3) + exp(2*b*(x1+x2))) / (-exp(2*b*x2) + exp(2*b*(x1+x3)))...
           + 2*a*b * cot(a*(x1-x2));

lam1_p = fsolve(@(lam) Df1(sqrt(lam), sqrt(K-lam)), 300);
lam2_p = fsolve(@(lam) Df2(sqrt(lam), sqrt(K-lam)), 300);

h0 = [L(2)/2; L; L(2)/2];
[U0, lam0] = eigP(K, h0, 2);
[u01, x01] = getvalP(U0(:,1), h0);
u01 = my_nmlz(u01);
[u02, x02] = getvalP(U0(:,2), h0);
u02 = my_nmlz(u02);

h1 = [L(2)/2+L(1)+L(2); L(3); L(4); L(5); L(2)/2];
[U1, lam1] = eigP(K, h1, 1);
[u1, x1] = getvalP(U1, h1);
u1 = my_nmlz(u1);

h2 = [L(2)/2; L(1); L(2)+L(3)+L(4)+L(5)+L(2)/2];
[U2, lam2] = eigP(K, h2, 1);
[u2, x2] = getvalP(U2, h2);
u2 = my_nmlz(u2);

figure
hold on
plot(x01, u01, '.')
plot(x02, u02, '.')
plot(x1, u1)
plot(x2, u2)
legend('eigen 1', 'eigen 2', 'eigen 1 under V1',...
    'eigen 1 under V2')

disp(lam0)
disp(lam1)
disp(lam2)
disp(lam1_p)
disp(lam2_p)

%%
all_K = 1100:10:1500;

L = [7, 20, 5, 1, 5];
L = L(:) / (sum(L) + L(2));
h0 = [L(2)/2; L; L(2)/2];
h1 = [L(2)/2+L(1)+L(2); L(3); L(4); L(5); L(2)/2];
h2 = [L(2)/2; L(1); L(2)+L(3)+L(4)+L(5)+L(2)/2];

lam01 = zeros(1, length(all_K));
lam02 = zeros(1, length(all_K));
lam1 = zeros(1, length(all_K));
lam2 = zeros(1, length(all_K));
for jk = 1:length(all_K)
    K = all_K(jk);
    [~, lam0] = eigP(K, h0, 2);
    lam01(jk) = lam0(1);
    lam02(jk) = lam0(2);
    [~, lam1(jk)] = eigP(K, h1, 1);
    [~, lam2(jk)] = eigP(K, h2, 1);
end
figure
hold on
plot(all_K, lam01, 's')
plot(all_K, lam02, 'o')
plot(all_K, lam1, '-')
plot(all_K, lam2, '-')
legend('eigen 1', 'eigen 2', 'eigen 1 under V1',...
    'eigen 1 under V2')
