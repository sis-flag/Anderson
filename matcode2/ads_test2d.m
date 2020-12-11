clear
rng(0)

% paremeters
K = 0;
B = [100; 200];
V = K * ones(5, 5);
h = 1;

x0 = 0.5;
y0 = 0.5;
alpha = 300;
beta = 0;
gamma = 0.1;

mymesh = mesh2d(size(V,1), size(V,2)); 

prob = ads_prob2d(V, B, 0, 0);
W1 = solve2d_enf(prob, mymesh, x0, y0, 0);
[x, y, w1] = getval2d(W1, mymesh);

prob = ads_prob2d(V, B, 0, 0);
W2 = solve2d_enf(prob, mymesh, x0, y0, 1/alpha);
[~, ~, w2] = getval2d(W2, mymesh);

prob = ads_prob2d(V, B, h*beta, h*gamma);
W3 = solve2d_enf(prob, mymesh, x0, y0, 0);
[~, ~, w3] = getval2d(W3, mymesh);

prob = ads_prob2d(V, B, h, 0);
W0 = solve2d_enf(prob, mymesh);
[x, y, w0] = getval2d(W0, mymesh);

prob = ads_prob2d(V, B, h, 0);
[U, lam] = eig2d(prob, mymesh);

figure
subplot(2,2,1)
cf = pcolor(x, y, w1);
cf.LineStyle = 'None';
colorbar
subplot(2,2,2)
cf = pcolor(x, y, w2);
cf.LineStyle = 'None';
colorbar
subplot(2,2,3)
cf = pcolor(x, y, w3);
cf.LineStyle = 'None';
colorbar
subplot(2,2,4)
cf = pcolor(x, y, w0);
cf.LineStyle = 'None';
colorbar

disp(lam)

figure
for k = 1:4
    [x, y, u] = getval2d(U(:,k), mymesh);
    u = abs(u) / max(max(abs(u)));
    
    land = abs(lam(k))*w1 + alpha*w2 + w3/(beta*min(min(w3))+gamma);
    
    land0 = abs(lam(k)) * w0;
    
    if all(land < u)
        fprintf('eigen %d, new failed\n', k)
    end
    
    if all(land0 < u)
        fprintf('eigen %d, old failed\n', k)
    end
    
    subplot(2,2,k)
    cf = pcolor(x, y, u);
    cf.LineStyle = 'none';
    colorbar
end