function du = getdiff1d(U, Ns, N)

if nargin <= 2
    N = 10;
end
M = round(length(U)-1) / N;

    function y = dphi(n, x)
        if n == 1
            y = -1/2;
        elseif n == N+1
            y = 1/2;
        else
            p1 = legendre(n-1,x);
            y = p1(1,:);
        end
    end

xhat = linspace(-1, 1, Ns +1);
dyhat = zeros(N+1, Ns+1);
for n = 1:N+1
    dyhat(n,:) = 2*M * dphi(n, xhat);
end

du = zeros(1, Ns*M+1);
for m = 1:M
    Uloc = U((m-1)*N+1: m*N+1);
    du((m-1)*Ns+1: m*Ns+1) = 0;
    for n = 1:N+1
        du((m-1)*Ns+1: m*Ns+1) = ...
            du((m-1)*Ns+1: m*Ns+1) + Uloc(n) * dyhat(n,:);
    end
end
end
