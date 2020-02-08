function [U, lam] = eigD2d(V, N, num)
% solve anderson eigen problem
% - laplace u(x) + V(x) u(x) = lam u(x) for x in [0,1] x [0,1]
% u'(x) + h u(x) = 0 for x on boundary
% V(x) is piecewise constant

if nargin <= 3
    num = 6;
end
if nargin <= 2
    N = 10;
end

M = size(V, 1);
hm = 1 / M;
[KA, KB] = lgmat2d(N);

    function u = l2g(m1, m2)
        [n2, n1] = meshgrid(1:(N+1), 1:(N+1));
        n1 = reshape(n1, 1, (N+1)*(N+1));
        n2 = reshape(n2, 1, (N+1)*(N+1));
        u1 = (m1-1)*N + n1;
        u2 = (m2-1)*N + n2;
        u = (u2-1)*(M*N+1) + u1;
    end

A = sparse((M*N+1)^2, (M*N+1)^2);
B = sparse((M*N+1)^2, (M*N+1)^2);
for m1 =1:M
    for m2 = 1:M
        Ae = KA + hm*hm/4 * V(m1,m2) * KB;
        Be = hm*hm/4 * KB;
        
        ind = l2g(m1, m2);
        
        A(ind, ind) = A(ind, ind) + Ae;
        B(ind, ind) = B(ind, ind) + Be;
    end
end

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