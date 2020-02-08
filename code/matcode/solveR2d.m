function [U] = solveR2d(V, g, x0, y0, u0, N)

% solve anderson source problem
% - laplace u(x) + V(x) u(x) = 1 for x in [0,1] x [0,1]
% u'(x) = g for x on boundary
% enforce u(x0, y0) = u0
% V(x) is piecewise constant

if nargin <= 5
    N = 10;
end
if nargin <= 4
    x0 = []; y0 = []; u0 = [];
end

M = size(V, 1);
hm = 1 / M;
[KA, KB, KF, ~, ~, ~, ~, KGF0, KGF1, KFG0, KFG1] = lgmat2d(N);

    function u = l2g(m1, m2)
        [n2, n1] = meshgrid(1:(N+1), 1:(N+1));
        n1 = reshape(n1, 1, (N+1)*(N+1));
        n2 = reshape(n2, 1, (N+1)*(N+1));
        u1 = (m1-1)*N + n1;
        u2 = (m2-1)*N + n2;
        u = (u2-1)*(M*N+1) + u1;
    end

A = sparse((M*N+1)^2, (M*N+1)^2);
F = sparse((M*N+1)^2, 1);
for m1 =1:M
    for m2 = 1:M
        Ae = KA + hm*hm/4 * V(m1,m2) * KB;
        Fe = hm*hm/4 * KF;
        
        ind = l2g(m1, m2);
        
        A(ind, ind) = A(ind, ind) + Ae;
        F(ind) = F(ind) + Fe;
    end
end

% left boundary
m2 = 1;
for m1 = 1:M
    ind = l2g(m1, m2);
    F(ind) = F(ind) + hm/2 * g * KGF0;
end

% right boundary
m2 = M;
for m1 = 1:M
    ind = l2g(m1, m2);
    F(ind) = F(ind) + hm/2 * g * KGF1;
end

% bottom boundary
m1 = 1;
for m2 = 1:M
    ind = l2g(m1, m2);
    F(ind) = F(ind) + hm/2 * g * KFG0;
end

% top boundary
m1 = M;
for m2 = 1:M
    ind = l2g(m1, m2);
    F(ind) = F(ind) + hm/2 * g * KFG1;
end

if ~isempty(x0)
    m1 = round(x0/hm); m2 = round(y0/hm);
    mm = l2g(m1, m2); mm = mm(1);
    A(mm, :) = 0;
    A(mm, mm) = 1;
    F(mm) = u0;
end

U = A \ full(F);
end


