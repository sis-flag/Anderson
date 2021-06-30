function [y0, x0] = my_resamp(y, x, num)

dy = diff(y); dx = diff(x);
ds = sqrt(dy.^2 + dx.^2);
s = [0, cumsum(ds)];

s0 = linspace(0, s(end), num);
x0 = interp1(s, x, s0, 'pchip');
y0 = interp1(x, y, x0, 'pchip');
