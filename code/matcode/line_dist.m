function d = line_dist(x1,y1,x2,y2)
dx = x1' - x2; dy = y1' - y2;
dd = dx.^2 + dy.^2;
d1 = max(min(dd)); d2 = min(max(dd));
d = min(d1,d2);
end