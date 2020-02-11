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
[iAhat, jAhat, vAhat] = find(Ahat);
[iBhat, jBhat, vBhat] = find(Bhat);
[iFhat, ~, vFhat] = find(Fhat);

l2g = @(m, n) (m-1) * N + n;

iA = zeros(1, M*(4*N+8));
jA = zeros(1, M*(4*N+8));
vA = zeros(1, M*(4*N+8));

kA = 0;
for m =1:M
    iA(kA+1:kA+N+3) = l2g(m, iAhat);
    jA(kA+1:kA+N+3) = l2g(m, jAhat);
    vA(kA+1:kA+N+3) = 2/hm * vAhat;
    kA = kA+N+3;
    
    iA(kA+1:kA+3*N+5) = l2g(m, iBhat);
    jA(kA+1:kA+3*N+5) = l2g(m, jBhat);
    vA(kA+1:kA+3*N+5) = hm/2 * V(m) * vBhat;
    kA = kA+3*N+5;
end

A = sparse(iA, jA, vA, M*N+1, M*N+1);

F = zeros(M*N+1, 1);
for m =1:M
    ind = l2g(m, iFhat);
    F(ind) = F(ind) + hm/2 * vFhat;
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