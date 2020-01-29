function [U, lam] = eigR1d(V, h, N, num)

% solve anderson eigen problem
% - u''(x) + V(x) u(x) = lam u(x) for x in [0,1]
% u'(x) + h u(x) = 0 for x=0 or x=1
% V(x) is piecewise constant

if nargin <= 3
    num = 6;
end
if nargin <= 2
    N = 10;
end

M = length(V);
hm = 1 / M;
[Ahat, Bhat, ~, ~, ~, ~, ~] = lgmat(N);

l2g = @(m) (m-1) * N + (1:N+1);

A = sparse(M*N+1, M*N+1);
B = sparse(M*N+1, M*N+1);
for m =1:M
    Ae = 2/hm * Ahat + hm/2 * V(m) * Bhat;
    Be = hm/2 * Bhat;
    
    ind = l2g(m);
    
    A(ind, ind) = A(ind, ind) + Ae;
    B(ind, ind) = B(ind, ind) + Be;
end

A(1, 1) = A(1, 1) + h;
A(end, end) = A(end, end) + h;

[U, lam] = eigs(A, B, num, 0);
lam = diag(lam);
end