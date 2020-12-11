function [u, x] = getvalP(U, h, Ns)
% 非均匀网格下求函数值

N = 10;
M = length(h);
if nargin < 3
    Ns = 30;
end

xm = [0; cumsum(h)];

U = reshape(U, M*N+1, 1);

xhat = linspace(-1, 1, Ns +1);
yhat = basis(N, xhat);

u = zeros(1, Ns*M+1);
x = zeros(1, Ns*M+1);
for m = 1:M
    Uloc = U((m-1)*N+1: m*N+1);
    u((m-1)*Ns+1: m*Ns+1) = yhat * Uloc;
    x((m-1)*Ns+1: m*Ns+1) = h(m)/2 * (xhat+1) + xm(m);
end
end
