function [U, lam] = eigD2d(V, N, num)
% solve anderson eigen problem
% - laplace u(x) + V(x) u(x) = lam u(x) for x in [0,1] x [0,1]
% u'(x) + h u(x) = 0 for x on boundary
% V(x) is piecewise constant

if nargin <= 2
    num = 6;
end
if nargin <= 1
    N = 10;
end

M = size(V, 1);
hm = 1 / M;

[KA, KB] = lgmat2d(N);
[iKA, jKA, vKA] = find(KA);
[iKB, jKB, vKB] = find(KB);

    function u = l2g(m1, m2, n)
        n1 = mod(n-1, N+1) +1;
        n2 = ceil(n / (N+1));
        u1 = (m1-1)*N + n1;
        u2 = (m2-1)*N + n2;
        u = (u2-1)*(M*N+1) + u1;
    end

nnzKA = length(iKA); nnzKB = length(iKB);

iA = zeros(1, M*M*(nnzKA+nnzKB));
jA = zeros(1, M*M*(nnzKA+nnzKB));
vA = zeros(1, M*M*(nnzKA+nnzKB));
iB = zeros(1, M*M*(nnzKB));
jB = zeros(1, M*M*(nnzKB));
vB = zeros(1, M*M*(nnzKB));

kA = 0; kB = 0;
for m1 =1:M
    for m2 = 1:M
        iA(kA+1:kA+nnzKA) = l2g(m1, m2, iKA);
        jA(kA+1:kA+nnzKA) = l2g(m1, m2, jKA);
        vA(kA+1:kA+nnzKA) = vKA;
        kA = kA+nnzKA;

        iA(kA+1:kA+nnzKB) = l2g(m1, m2, iKB);
        jA(kA+1:kA+nnzKB) = l2g(m1, m2, jKB);
        vA(kA+1:kA+nnzKB) = hm*hm/4 * V(m1,m2) * vKB;
        kA = kA+nnzKB;

        iB(kB+1:kB+nnzKB) = l2g(m1, m2, iKB);
        jB(kB+1:kB+nnzKB) = l2g(m1, m2, jKB);
        vB(kB+1:kB+nnzKB) = hm*hm/4 * vKB;
        kB = kB+nnzKB;
    end
end

A = sparse(iA, jA, vA, (M*N+1)^2, (M*N+1)^2);
B = sparse(iB, jB, vB, (M*N+1)^2, (M*N+1)^2);

tmp = reshape(1:(M*N+1)^2, (M*N+1), (M*N+1));
edge = [tmp(1,:), tmp(end,:), tmp(:,1)', tmp(:,end)'];
A(edge,:) = []; A(:,edge) = [];
B(edge,:) = []; B(:,edge) = [];

[UU, lam] = eigs(A, B, num, 0);
lam = diag(lam);

U = zeros((M*N+1)^2, num);
ind = setdiff(tmp(:), edge);
U(ind,:) = UU;
end