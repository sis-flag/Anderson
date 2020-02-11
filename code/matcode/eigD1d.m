function [U, lam] = eigD1d(V, N, num)

% solve anderson eigen problem
% - u''(x) + V(x) u(x) = lam u(x) for x in [0,1]
% u(x) = 0 for x=0 or x=1
% V(x) is piecewise constant

if nargin <= 2
    num = 6;
end
if nargin <= 1
    N = 10;
end

M = length(V);
hm = 1 / M;

[Ahat, Bhat] = lgmat(N);
[iAhat, jAhat, vAhat] = find(Ahat);
[iBhat, jBhat, vBhat] = find(Bhat);

l2g = @(m, n) (m-1) * N + n;

iA = zeros(1, M*(4*N+8));
jA = zeros(1, M*(4*N+8));
vA = zeros(1, M*(4*N+8));
iB = zeros(1, M*(3*N+5));
jB = zeros(1, M*(3*N+5));
vB = zeros(1, M*(3*N+5));

kA = 0; kB = 0;
for m =1:M
    iA(kA+1:kA+N+3) = l2g(m, iAhat);
    jA(kA+1:kA+N+3) = l2g(m, jAhat);
    vA(kA+1:kA+N+3) = 2/hm * vAhat;
    kA = kA+N+3;
    
    iA(kA+1:kA+3*N+5) = l2g(m, iBhat);
    jA(kA+1:kA+3*N+5) = l2g(m, jBhat);
    vA(kA+1:kA+3*N+5) = hm/2 * V(m) * vBhat;
    kA = kA+3*N+5;
    
    iB(kB+1:kB+3*N+5) = l2g(m, iBhat);
    jB(kB+1:kB+3*N+5) = l2g(m, jBhat);
    vB(kB+1:kB+3*N+5) = hm/2 * vBhat;
    kB = kB+3*N+5;
end

A = sparse(iA, jA, vA, M*N+1, M*N+1);
B = sparse(iB, jB, vB, M*N+1, M*N+1);

edge = [1, M*N+1];
A(edge,:) = []; A(:,edge) = [];
B(edge,:) = []; B(:,edge) = [];

[U, lam] = eigs(A, B, num, 0);
lam = diag(lam);
U = [zeros(1,num); U; zeros(1,num)];
end