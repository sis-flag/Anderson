function Pb = pred_prob(p, N)

q = 1 - p;
M = N * p * q;

b1 = 0;
for k = 2:1e10
    n = 1:(k-1);
    pow = 2*max(n, k-n)-1;
    x = q^(k-2) * sum((1-q.^pow).^(M-2)) * p^2;
    b1 = b1 + x;
    if x < 1e-10 && k > 100
        break
    end
end

b2 = 0;
for n = 1:1e10
    x = q^(n-1) * p * (1-q^(2*n-1))^(M-1);
    b2 = b2 + x;
    if x < 1e-10 && n > 100
        break
    end
end

Pb = q*q * b1 + 2*p*q * b2;

end
