function [U, lam] = eig1d(V, h, num, N)
% solve 1-d anderson eigen problem
% - u''(x) + V(x) u(x) = lam u(x) for x in [0, 1]
% u'(x) + h u(x) = 0 for x = 0 or x = 1
% V(x) is piecewise constant
% input:
%     V(1-d array):   piecewise constant of V(x)
%     h(number):      paremeter in equation, h=inf for Dirichlet boundary
%     num(integer):   number of eigenvalues required
%     N(integer):     degree of polynomials
% output:
%     U(2-d array):   array with size (num, M*N+1)
%                     each column represents projection on polynomial basis
%                     use function getval1d to get value of sulution
%     lam(1-d array): array with size (num, 1), eigenvalues

% When V(x) < 0, 'eigs' may not get the smallest eigenvalue
V = V - min(V(:)); % shift

% default input
if nargin < 4
    N = 10;
end
if nargin < 3
    num = 4;
end

M = length(V);
hm = 1 / M;

[Ahat, Bhat] = lgmat(N);
[iAhat, jAhat, vAhat] = find(Ahat);
[iBhat, jBhat, vBhat] = find(Bhat);

l2g = @(m, n) (m-1) * N + n;

nnzA = length(iAhat); nnzB = length(iBhat);
iA = zeros(1, M*(nnzA+nnzB));
jA = zeros(1, M*(nnzA+nnzB));
vA = zeros(1, M*(nnzA+nnzB));
iB = zeros(1, M*(nnzB));
jB = zeros(1, M*(nnzB));
vB = zeros(1, M*(nnzB));

kA = 0; kB = 0;
for m =1:M
    iA(kA+1:kA+nnzA) = l2g(m, iAhat);
    jA(kA+1:kA+nnzA) = l2g(m, jAhat);
    vA(kA+1:kA+nnzA) = 2/hm * vAhat;
    kA = kA+nnzA;
    
    iA(kA+1:kA+nnzB) = l2g(m, iBhat);
    jA(kA+1:kA+nnzB) = l2g(m, jBhat);
    vA(kA+1:kA+nnzB) = hm/2 * V(m) * vBhat;
    kA = kA+nnzB;
    
    iB(kB+1:kB+nnzB) = l2g(m, iBhat);
    jB(kB+1:kB+nnzB) = l2g(m, jBhat);
    vB(kB+1:kB+nnzB) = hm/2 * vBhat;
    kB = kB+nnzB;
end

A = sparse(iA, jA, vA, M*N+1, M*N+1);
B = sparse(iB, jB, vB, M*N+1, M*N+1);

if isinf(h)
    edg = [1, M*N+1];
    A(edg,:) = []; A(:,edg) = [];
    B(edg,:) = []; B(:,edg) = [];
    
    [U, lam] = eigs(A, B, num, 0);
    lam = diag(lam);
    
    U = [zeros(1,num); U; zeros(1,num)];
else
    A(1, 1) = A(1, 1) + h;
    A(end, end) = A(end, end) + h;
    
    [U, lam] = eigs(A, B, num, 0);
    lam = diag(lam) + min(V(:));
end

end