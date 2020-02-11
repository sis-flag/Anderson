function [U] = solveD2d(V, N)

% solve anderson source problem
% - laplace u(x) + V(x) u(x) = 1 for x in [0,1] x [0,1]
% u'(x) = g for x on boundary
% enforce u(x0, y0) = u0
% V(x) is piecewise constant

if nargin <= 1
    N = 10;
end

M = size(V, 1);
hm = 1 / M;
[KA, KB, KF] = lgmat2d(N);

[iKA, jKA, vKA] = find(KA);
[iKB, jKB, vKB] = find(KB);
[iKF, ~, vKF] = find(KF);

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
    end
end

A = sparse(iA, jA, vA, (M*N+1)^2, (M*N+1)^2);

F = zeros((M*N+1)^2, 1);
for m1 =1:M
    for m2 = 1:M
        ind = l2g(m1, m2, iKF);
        F(ind) = F(ind) + hm/hm/4 * vKF;
    end
end

tmp = reshape(1:(M*N+1)^2, (M*N+1), (M*N+1));
edge = [tmp(1,:), tmp(end,:), tmp(:,1)', tmp(:,end)'];
A(edge,:) = []; A(:,edge) = [];
F(edge) = [];

UU = A \ F;

U = zeros((M*N+1)^2, 1);
ind = setdiff(tmp(:), edge);
U(ind) = UU;
end


