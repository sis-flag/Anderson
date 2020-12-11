function [pn] = pred_2peak(p, N)

q = 1 - p;
M = N * p * q;

pn = 0;
for k = 1:1e10
    x = M * q * p^k * (1-p^k)^(M-1);
    pn = pn + x;
    if x < 1e-10 && k > 100
        break
    end
end

pn = 1 - pn;
end