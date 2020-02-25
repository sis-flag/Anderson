clear;
rng(123);

% peremeters
h = 0;
V = rand(10)*1000;

W = solveD2d(V);
w = getval2d(W, 50);

x = linspace(0,1,length(w));
[x2, x1] = meshgrid(x, x); % caution!

%% watershed
N = size(w,1);
water = zeros(size(w));
[~, ind] = sort(-w(:));
shedi = []; shedj = [];
for k = 1:length(ind)
    water(ind(k)) = 1;
    
    ti = mod(ind(k)-1, N) +1;
    tj = ceil(ind(k) / N);
    if ti == 1 || ti == N || tj == 1 || tj == N
        continue
    end
    wloc = water(ti-1:ti+1, tj-1:tj+1);
    if wloc(2,1)+wloc(1,2)+wloc(2,3)+wloc(3,2) >= 3
        shedi = [shedi, ti];
        shedj = [shedj, tj];
    end
end

figure();
hold on
s = pcolor(x1,x2,w);
s.LineStyle = 'none';
colorbar;
plot(shedi/N, shedj/N, 'r.')