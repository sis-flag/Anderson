clear;
rng(123);

% peremeters
h = 0;
V = rand(20)*1000;

W = solveR2d(V, 0);
w = getval2d(W, 20);

slxx = linspace(0, 1, size(w,1));
[x, y] = meshgrid(slxx, slxx);

%% valley line
[dwx, dwy] = gradient(w, 1/size(w,1));
dw = sqrt(dwx.^2 + dwy.^2);

% find local minimum of |grad w|
cprg = round(0.2 * size(dw,1)/size(V,1));
mdw = mean(mean(dw)); maxw = max(max(w));
sx = []; sy = []; sdw = []; sw = [];
for k1 = 1+cprg:size(w,1)-cprg
    for k2 = 1+cprg:size(w,1)-cprg
        dwloc = dw(k1-cprg:k1+cprg, k2-cprg:k2+cprg);
        dwn = sum(sum(dwloc < dw(k1,k2)));
        if (dwn == 0) && (dw(k1,k2) < 0.1*mdw) && (w(k1,k2) < 0.85*maxw)
            sx = [sx, x(k1,k2)];
            sy = [sy, y(k1,k2)];
            sdw = [sdw, dw(k1,k2)];
            sw = [sw, w(k1,k2)];
        end
    end
end

% delete local minimum
tx1 = stream2(x,y,-dwx,-dwy,sx,sy);
idx1 = zeros(1,length(tx1),'logical');
for k = 1:length(tx1)
    tmp = tx1{k};
    rg = norm(tmp(end,:) - tmp(1,:));
    idx1(k) = (rg > 0.01);
end

% delete local maximum
tx2 = stream2(x,y,dwx,dwy,sx,sy);
idx2 = zeros(1,length(tx2),'logical');
for k = 1:length(tx2)
    tmp = tx2{k};
    rg = norm(tmp(end,:) - tmp(1,:));
    idx2(k) = (rg > 0.01);
end

idx = idx2 & idx1;
sx = sx(idx); sy = sy(idx);

% disturbance
NN = length(sx);
ssx = [sx+0.002, sx-0.002, sx, sx];
ssy = [sy, sy, sy-0.002, sy-0.002];

% get streamline
slxx = stream2(x,y,-dwx,-dwy,ssx,ssy);

% delete duplicate
for k = 1:NN
    x1 = slxx{k}(:,1);      y1 = slxx{k}(:,2);
    x2 = slxx{k+NN}(:,1);   y2 = slxx{k+NN}(:,2);
    x3 = slxx{k+2*NN}(:,1); y3 = slxx{k+2*NN}(:,2);
    x4 = slxx{k+3*NN}(:,1); y4 = slxx{k+3*NN}(:,2);
    
    d12 = line_dist(x1,y1,x2,y2);
    d13 = line_dist(x1,y1,x3,y3);
    d14 = line_dist(x1,y1,x4,y4);
    d23 = line_dist(x2,y2,x3,y3);
    d24 = line_dist(x2,y2,x4,y4);
    d34 = line_dist(x3,y3,x4,y4);
    
    rem = 4;
    if d12 < 1e-4
        slxx{k+NN} = [];
        rem = rem - 1;
    end
    if d13 < 1e-4 || d23 < 1e-4
        slxx{k+2*NN} = [];
        rem = rem - 1;
    end
    if d14 < 1e-4 || d24 < 1e-4 || d34 < 1e-4
        slxx{k+3*NN} = [];
        rem = rem - 1;
    end
    
    % if more than 2 lines remained
    if rem ~= 2
        slxx{k} = [];
        slxx{k+NN} = [];
        slxx{k+2*NN} = [];
        slxx{k+3*NN} = [];
    end
end

figure();
hold on
s = pcolor(x, y, w);
s.LineStyle = 'none';
colorbar;
for tmp = slxx
    if isempty(tmp{1})
        continue
    end
    xx =  tmp{1}(:,1); yy =  tmp{1}(:,2);
    plot(xx, yy, 'k-');
end


% maxw = max(sw); minw = min(min(w));
% rgw = maxw - minw;

% figure();
% hold on
% s = pcolor(x, y, w);
% s.LineStyle = 'none';
% colorbar;
% plot(ssx, ssy, 'r.')
% for tmp = slxx
%     if isempty(tmp{1})
%         continue
%     end
%     xx =  tmp{1}(:,1); yy =  tmp{1}(:,2);
%     ww = interp2(x,y,w,xx,yy,'nearest');
%     for kk = 1:floor(length(xx)/30)
%         txx = xx((kk-1)*30+1:kk*30);
%         tyy = yy((kk-1)*30+1:kk*30);
%         mtww = mean(ww((kk-1)*30+1:kk*30));
%         fj = plot(txx, tyy, 'k.');
%         fj.MarkerSize = round( (maxw-mtww)/rgw*8 );
%     end
%     txx = xx((kk-1)*30+1:end);
%     tyy = yy((kk-1)*30+1:end);
%     mtww = mean(ww((kk-1)*30+1:end));
%     fj = plot(txx, tyy, 'k.');
%     fj.MarkerSize = round( (maxw-mtww)/rgw*8 );
% end