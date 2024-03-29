function [U] = solve1d(V, h, N)
% solve 1-d anderson source problem
% - u''(x) + V(x) u(x) = 1 for x in [0, 1]
% u'(x) + h u(x) = 0 for x = 0 or x = 1
% V(x) is piecewise constant
% input:
%     V(1-d array):   piecewise constant of V(x)
%     h(number):      paremeter in equation, h=inf for Dirichlet boundary
%     N(integer):     degree of polynomials
% output:
%     U(1-d array):   array with size (1, M*N+1)
%                     each column represents projection on polynomial basis
%                     use function getval1d to get value of sulution

% default input
if nargin < 3
    N = 10;
end

M = length(V);
hm = 1 / M;

[Ahat, Bhat, Fhat] = lgmat(N);
[iAhat, jAhat, vAhat] = find(Ahat);
[iBhat, jBhat, vBhat] = find(Bhat);
[iFhat, ~, vFhat] = find(Fhat);

l2g = @(m, n) (m-1) * N + n;

nnzA = length(iAhat); nnzB = length(iBhat);
iA = zeros(1, M*(nnzA+nnzB));
jA = zeros(1, M*(nnzA+nnzB));
vA = zeros(1, M*(nnzA+nnzB));
F = zeros(M*N+1, 1);

kA = 0;
for m =1:M
    iA(kA+1:kA+nnzA) = l2g(m, iAhat);
    jA(kA+1:kA+nnzA) = l2g(m, jAhat);
    vA(kA+1:kA+nnzA) = 2/hm * vAhat;
    kA = kA+nnzA;
    
    iA(kA+1:kA+nnzB) = l2g(m, iBhat);
    jA(kA+1:kA+nnzB) = l2g(m, jBhat);
    vA(kA+1:kA+nnzB) = hm/2 * V(m) * vBhat;
    kA = kA+nnzB;

    ind = l2g(m, iFhat);
    F(ind) = F(ind) + hm/2 * vFhat;
end

A = sparse(iA, jA, vA, M*N+1, M*N+1);

if isinf(h)
    edg = [1, M*N+1];
    A(edg,:) = []; A(:,edg) = [];
    F(edg) = [];

    U = [0; A \ F; 0];
else
    A(1, 1) = A(1, 1) + h;
    A(end, end) = A(end, end) + h;

    U = A \ F;
end