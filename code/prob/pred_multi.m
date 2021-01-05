function [PD, PN] = pred_multi(p, N)

q = 1 - p;
M = N * p * q;

PD = 1;
for n = 1:1e10
    x = M * (1-q^(n-1))^(M-1) * q^(n-1) * p;
    PD = PD - x;
    if x < 1e-10 && n > 100
        break
    end
end

PN = 1;
for n = 1:1e10
    x = q^2 * (M-2) * (1-q^floor((n-1)/2))^2 * (1-q^(n-1))^(M-3) * q^(n-1) * p...
      + 2 * q^2 * (1 - q^(2*n-1))^(M-2) * (1-q^(n-1)) * q^(n-1) * p...
      + 2 * p * q * (M-1) * (1-q^floor((n-1)/2)) * (1 - q^(n-1))^(M-2) * q^(n-1) * p...
      + 2 * p * q * (1 - q^(2*n-1))^(M-1) * q^(n-1) * p...
      + p^2 * M * (1-q^(n-1))^(M-1) * q^(n-1) * p;
    PN = PN - x;
    if x < 1e-10 && n > 100
        break
    end
end

end
