clear
rng(0)

% paremeters
K = 1;
B = 200;
V = K * ones(1,10);
h = 1;

x0 = [];
alpha = 0;
beta = 10;
gamma = 0;

% % plot potential
% figure
% cf = bar(((1:length(V))-0.5)/length(V), V);
% cf.BarWidth = 1;

mymesh = mesh1d(length(V));

prob = ads_prob1d(V, B, 0, 0);
W1 = solve1d_enf(prob, mymesh, x0, 0);
[x, w1] = getval1d(W1, mymesh);

prob = ads_prob1d(V, B, 0, 0);
W2 = solve1d_enf(prob, mymesh, x0, 1/alpha);
[x, w2] = getval1d(W2, mymesh);

prob = ads_prob1d(V, B, h*beta, h*gamma);
W3 = solve1d_enf(prob, mymesh, x0, 0);
[x, w3] = getval1d(W3, mymesh);

prob = ads_prob1d(V, B, h, 0);
W0 = solve1d_enf(prob, mymesh, x0, 0);
[x, w0] = getval1d(W0, mymesh);

prob = ads_prob1d(V, B, h, 0);
[U, lam] = eig1d(prob, mymesh);

disp(lam)

xx = linspace(0, 1, 1001);
xx = xx(2:end-1);
ax = arrayfun(prob.a, xx);
bx = arrayfun(prob.b, xx);
cx = arrayfun(prob.c, xx);
figure
subplot(1,3,1)
plot(xx, ax)
subplot(1,3,2)
plot(xx, bx)
subplot(1,3,3)
plot(xx, cx)

figure
hold on
plot(x, w1, 'm-')
plot(x, w2, 'r-')
plot(x, w3, 'b-')
plot(x, w0, 'k-')
legend('w1', 'w2', 'w3', 'w(old)')

figure
for k = 1:4
    [x, u] = getval1d(U(:,k), mymesh);
    u = abs(u) / max(abs(u));
    
    landscape = abs(lam(k))*w1 + alpha*w2 + w3/(beta*min(w3)+gamma);
    
    subplot(2,2,k);
    hold on
    plot(x, landscape, '-')
    plot(x, w0 * abs(lam(k)), '-')
    plot(x, u, '--')
    legend('landscape', 'landscape(old)', 'eigen function')
end