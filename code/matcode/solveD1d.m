function [U] = solveD1d(V, N)

% solve anderson source problem
% - u''(x) + V(x) u(x) = 1 for x in [0,1]
% u(x) = 0 for x=0 or x=1
% V(x) is piecewise constant

if nargin <= 1
    N = 10;
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

edge = [1, M*N+1];
A(edge,:) = []; A(:,edge) = [];
F(edge) = [];

U = A \ full(F);
U = [0; U; 0];
end