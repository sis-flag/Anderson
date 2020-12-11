function [b1, b2] = pred_prob(p, N)

q = 1 - p;
M = N * p * q;

b1 = 0;
for k = 2:1e10
    n = 1:(k-1);
    pow = 2*max(n, k-n)-1;
    x = p^(k-2) * sum((1-p.^pow).^(M-2)) * q^2;
    b1 = b1 + x;
    if x < 1e-10 && k > 100
        break
    end
end

b2 = 0;
for n = 1:1e10
    x = p^(n-1) * q * (1-p^(2*n-1))^(M-1);
    b2 = b2 + x;
    if x < 1e-10 && n > 100
        break
    end
end
