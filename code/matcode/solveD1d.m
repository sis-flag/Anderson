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

edge = [1, M*N+1];
A(edge,:) = []; A(:,edge) = [];
F(edge) = [];

U = A \ F;
U = [0; U; 0];
end