function [U] = solveR1d(V, g, x0, u0, N)

% solve anderson source problem
% - u''(x) + V(x) u(x) = 1 for x in [0,1]
% u'(x) = g for x=0 or x=1
% enforce u(x0) = u0
% V(x) is piecewise constant

if nargin <= 4
    N = 10;
end
if nargin <= 3
    x0 = []; u0 = [];
end

M = length(V);
hm = 1 / M;
[Ahat, Bhat, Fhat] = lgmat(N);

l2g = @(m) (m-1) * N + (1:N+1);

A = sparse(M*N+1, M*N+1);
F = sparse(M*N+1, 1);
for m =1:M
    Ae = 2/hm * Ahat + hm/2 * V(m) * Bhat;
    Fe = hm/2 * Fhat;
    
    ind = l2g(m);
    
    A(ind, ind) = A(ind, ind) + Ae;
    F(ind) = F(ind) + Fe;
end

F(1) = F(1) + g;
F(end) = F(end) + g;

if ~isempty(x0)
    m = ceil(x0/hm);
    mm = (m-1)*N+1;
    A(mm, :) = 0;
    A(mm, mm) = 1;
    F(mm) = u0;
end

U = A \ full(F);
end