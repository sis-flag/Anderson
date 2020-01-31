clear
load 2.mat

xx = linspace(0, 1, size(w,1));
[x, y] = meshgrid(xx, xx);
[dwx, dwy] = gradient(w, 1/size(w,1));
[dwxx, dwyx] = gradient(dwx, 1/size(w,1));
[dwxy, dwyy] = gradient(dwy, 1/size(w,1));

% dw = dwx.^2 + dwy.^2;
% delw = dwxx.*dwyy - dwyx.*dwxy;
% [wi, wj] = find(dw < 1e-5 & delw < -2e-2);
% sx = xx(wj);
% sy = xx(wi);

sxx = linspace(0,1,15);
[sx, sy] = meshgrid(sxx(2:end), sxx(2:end));
sx = reshape(sx, 1, []);
sy = reshape(sy, 1, []);

hold on
quiver(x(1:25:end,1:25:end), y(1:25:end,1:25:end),...
      dwx(1:25:end,1:25:end), dwy(1:25:end,1:25:end));

% plot(sx, sy, 'r.')
% 
% if length(sx) < 3000
%     streamline(x,y,-dwx,-dwy,sx,sy);
% end


