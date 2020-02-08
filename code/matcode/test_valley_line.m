clear;
rng(123);

% peremeters
h = 0;
V = rand(20)*3000;

W = solveR2d(V, 0);
[U, lam] = eigR2d(V, 0);

w = getval2d(W, 100);

slxx = linspace(0, 1, size(w,1));
[x, y] = meshgrid(slxx, slxx);

%% valley line
[dwx, dwy] = gradient(w, 1/size(w,1));
dw = sqrt(dwx.^2 + dwy.^2);

% sx = []; sy = [];
% for k1 = 2:size(w,1)-1
%     for k2 = 2:size(w,2)-1
%         wloc = w(k1-1:k1+1, k2-1:k2+1);
%         dwloc = dw(k1-1:k1+1, k2-1:k2+1);
%         wn = sum(sum(wloc < w(k1,k2)));
%         dwn = sum(sum(dwloc < dw(k1,k2)));
%         if (2 <= wn) && (wn <= 6) && (dwn <= 0) && (dw(k1,k2)<0.001)
%             sx = [sx, x(k1,k2)];
%             sy = [sy, y(k1,k2)];
%         end
%     end
% end

cprg = 5;
mw = mean(mean(dw));
sx = []; sy = [];
for k1 = 1+cprg:size(w,1)-cprg
    for k2 = 1+cprg:size(w,1)-cprg
        dwloc = dw(k1-cprg:k1+cprg, k2-cprg:k2+cprg);
        dwn = sum(sum(dwloc < dw(k1,k2)));
        if (dwn <= 0) && (dw(k1,k2) < 0.1*mw)
            sx = [sx, x(k1,k2)];
            sy = [sy, y(k1,k2)];
        end
    end
end

% tx1 = stream2(x,y,-dwx,-dwy,sx,sy);
% idx1 = zeros(1,length(tx1),'logical');
% for k = 1:length(tx1)
%     tmp = tx1{k};
%     rg = norm(tmp(end,:) - tmp(1,:));
%     idx1(k) = (rg > 0.01);
% end

% delete local maximum
tx2 = stream2(x,y,dwx,dwy,sx,sy);
idx2 = zeros(1,length(tx2),'logical');
for k = 1:length(tx2)
    tmp = tx2{k};
    rg = norm(tmp(end,:) - tmp(1,:));
    idx2(k) = (rg > 0.01);
end

idx = idx2;
sx = sx(idx); sy = sy(idx);

ssx = [sx+0.002, sx-0.002, sx+0.002, sx-0.002];
ssy = [sy+0.002, sy+0.002, sy-0.002, sy-0.002];

figure();
hold on
s = pcolor(x, y, w);
s.LineStyle = 'none';
colorbar;
% plot(ssx, ssy, 'r.')
slxx = stream2(x,y,-dwx,-dwy,ssx,ssy);
for k = 1:length(slxx)
    tmp = slxx{k};
    xx =  tmp(:,1); yy =  tmp(:,2);
    ww = interp2(x,y,w,xx,yy,'nearest');
    ind = (ww<0.5/lam(1));
    plot(xx(ind),yy(ind),'k.')
%      plot(xx,yy,'k.')
end