function [n, mu] = nlocmax(u, thes)
n = [];
for k = 2:length(u)-1
    if u(k) > thes && u(k-1) < u(k) && u(k) > u(k+1)
        n = [n; k];
    end
end
if u(1) > thes && u(1) > u(2)
    n = [n; 1];
end
if u(end) > thes && u(end) > u(end-1)
    n = [n; length(u)];
end
[mu, ind] = sort(u(n), 'descend');
n = n(ind);