clear

%%
K = 500;
L = [7, 20, 5, 1, 5];

L = L(:) / (sum(L) + L(2));
h = [L(2)/2; L; L(2)/2];
[~, lam0] = eigP(K, h, 3);

tlam = lam0(1)-10:lam0(end)+50;
func = @(lam) Dfunc(K, h, lam);
tD = arrayfun(func, tlam);

figure
hold on
plot(tlam, tD)
plot(tlam, zeros(size(tlam)))
plot(lam0, zeros(size(lam0)), 'k*')

disp(lam0);

%%
L = [7, 20, 5, 1, 5];
L = L(:) / (sum(L) + L(2));

h = [L(2)/2; L; L(2)/2];

figure
hold on
for K = 1310:20:1390
    [~, lam0] = eigP(K, h, 2);
    
    tlam = lam0(1)-1:0.1:lam0(2)+1;
    func = @(lam) Dfunc(K, [L; L(2)], lam);
    tD = arrayfun(func, tlam);
    plot(tlam, tD)
end
plot([305, 320], [0, 0], 'k-')
