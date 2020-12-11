function U = solve1d_enf(prob, mesh, x0, u0, N)

% defult input
if nargin < 5
    N = 6;
end
if nargin < 4
    x0 = []; u0 = [];
end
M = mesh.M;
h = mesh.h;

% Gauss nodes
[xn, wn] = lgpoints(N);

% basis on Gauss nodes
[phi, dphi] = basis(N, xn);

% local to global
if prob.bd == 'P'
    l2g = @(m, n) mod((m-1) * N + (n-1), M*N) + 1;
else
    l2g = @(m, n) (m-1) * N + n;
end

A = zeros(M*N+1, M*N+1);
F = zeros(M*N+1, 1);

for m =1:M
    lxn = (xn + 1)/2 * h(m) + mesh.x(m);
    lan = arrayfun(prob.a, lxn);
    lbn = arrayfun(prob.b, lxn);
    lcn = arrayfun(prob.c, lxn);
    lfn = arrayfun(prob.f, lxn);
    
    lA = (2/h(m)) * dphi' * diag(wn.*lan) * dphi ...
       + phi' * diag(wn.*lbn) * dphi ...
       + (h(m)/2) * phi' * diag(wn.*lcn) * phi;
       
    lF = (h(m)/2) * phi' * (lfn.*wn);
        
%     [liA, ljA, lvA] = find(lA);
%     iA(nnzA+1: nnzA+(N+1)^2) = l2g(m, liA);
%     jA(nnzA+1: nnzA+(N+1)^2) = l2g(m, ljA);
%     vA(nnzA+1: nnzA+(N+1)^2) = lvA;
%     nnzA = nnzA + (N+1)^2;
    
    gi = l2g(m, 1:N+1);
    A(gi, gi) = A(gi, gi) + lA;
    F(gi) = F(gi) + lF;
end

% enforce
if ~isempty(x0)
    for m = 1:M-1
        if mesh.x(m) <= x0 && x0 < mesh.x(m+1)
            m0 = m;
            break
        end
    end
    mm = l2g(m0, 1);
    A(mm, :) = 0;
    A(mm, mm) = 1;
    F(mm) = u0;
end

if prob.bd == 'D'
    edge = [1, M*N+1];
    A(edge,:) = []; A(:,edge) = [];
    F(edge) = [];
    
    U = A \ F;
    U = [0; U; 0];
    
elseif prob.bd == 'R'
    A(1,1) = A(1,1) + prob.h(0);
    A(end,end) = A(end,end) + prob.h(1);
    F(1,1) = F(1,1) + prob.g(0);
    F(end,end) = F(end,end) + prob.g(1);
    
    U = A \ F;
    
elseif prob.bd == 'P'
    
    U = A \ F;
    U = [U; U(1)];
end

end
