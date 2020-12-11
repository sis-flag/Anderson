clear 

N = [7, 20, 5, 1, 5];

Df = @(x) Dfunc1(x(1), x(2), N);
res = fsolve(Df, [17; 32]);
Kc_p = res(1)^2 + res(2)^2;
lamc_p = res(1)^2;

[Kc_s, b, lamc_s] = getP_Kc(N);
