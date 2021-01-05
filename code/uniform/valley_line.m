function [vx, vy, vw] = valley_line(w)
% get coordinate of points on valley line
% input:
%     w(2-d array):   value of function on mesh points
% output:
%     vx, vy(1-d array): x and y coordinate of points on valley line
%     vw(1-d array):     function value on (vx, vy)

mark = my_watershed(-w);
[vi, vj] = find(isnan(mark));
vx = (vi-1) / size(w,1);
vy = (vj-1) / size(w,1);

tw = w(:);
vw = tw(vi+vj*size(w,1));

end