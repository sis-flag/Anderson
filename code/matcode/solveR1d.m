function [U] = solveR1d(V, g, x0, u0, N)

% solve anderson source problem
% - u''(x) + V(x) u(x) = 1 for x in [0,1]
% u'(x) = g for x=0 or x=1
% enforce u(x0) = u0
% V(x) is piecewise constant

if nargin <= 2
    N = 5;
end

M = length(V);
hm = 1 / M;
[Ahat, Bhat, ~, ~, ~, ~, ~] = lgmat(N);

l2g = @(m) (m-1) * N + (1:N+1);

A = sparse(M*N+1, M*N+1);
F = sparse(M*N+1, 1);
for m =1:M
    Ae = 2/hm * Ahat + hm/2 * V(m) * Bhat;
    Fe = hm/2 * F_hat;
    
    ind = l2g(m);
    
    A(ind, ind) = A(ind, ind) + Ae;
    F(ind) = F(ind) + Fe;
end

F(1) = F(1) + g;
F(end) = F(end) + g;

if ~isempty(x0)
    m = ceil(x0/hm);
    xx = mod(x0, hm)/hm;
    P = legendre(N, xx);
    A = 666;
    F = 777;
end

U = A \ F;
end