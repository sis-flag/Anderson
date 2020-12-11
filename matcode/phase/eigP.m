function [U, lam] = eigP(K, h, num)
% 专门用于研究相变的一维非均匀网格代码

N = 10;
M = length(h);
V = zeros(M, 1);
for kk = 1:2:M
    V(kk) = K;
end

[Ahat, Bhat] = lgmat(N);
[iAhat, jAhat, vAhat] = find(Ahat);
[iBhat, jBhat, vBhat] = find(Bhat);

    function ind = l2g(m, n)
        ind = (m-1)*N + n;
        ind(ind==M*N+1) = 1;
    end

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
    vA(kA+1:kA+nnzA) = 2/h(m) * vAhat;
    kA = kA+nnzA;
    
    iA(kA+1:kA+nnzB) = l2g(m, iBhat);
    jA(kA+1:kA+nnzB) = l2g(m, jBhat);
    vA(kA+1:kA+nnzB) = h(m)/2 * V(m) * vBhat;
    kA = kA+nnzB;
    
    iB(kB+1:kB+nnzB) = l2g(m, iBhat);
    jB(kB+1:kB+nnzB) = l2g(m, jBhat);
    vB(kB+1:kB+nnzB) = h(m)/2 * vBhat;
    kB = kB+nnzB;
end

A = sparse(iA, jA, vA, M*N, M*N);
B = sparse(iB, jB, vB, M*N, M*N);

[UU, lam] = eigs(A, B, num, 0);
lam = diag(lam);

U = [UU; UU(1,:)];
end