function [n1, n2, mu] = nlocmax2d(u, thes)
[N1, N2] = size(u);
n1 = [];
n2 = [];
mu = [];
for k1 = 1:N1
    for k2 = 1:N2
        if u(k1,k2) < thes
            continue
        end
        if k1 > 1 && u(k1,k2) < u(k1-1,k2)
            continue
        end
        if k1 < N1 && u(k1,k2) < u(k1+1,k2)
            continue
        end
        if k2 > 1 && u(k1,k2) < u(k1,k2-1)
            continue
        end
        if k2 < N2 && u(k1,k2) < u(k1,k2+1)
            continue
        end
        n1 = [n1; k1];
        n2 = [n2; k2];
        mu = [mu; u(k1,k2)];
    end
end

[mu, ind] = sort(mu, 'descend');
n1 = n1(ind);
n2 = n2(ind);
