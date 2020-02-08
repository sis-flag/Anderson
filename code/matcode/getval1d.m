function u = getval1d(U, Ns, N)

if nargin <= 2
    N = 10;
end
M = round(length(U)-1) / N;

    function y = phi(n, x)
        if n == 1
            y = (1-x)/2;
        elseif n == N+1
            y = (1+x)/2;
        else
            p1 = legendre(n,x); p1 = p1(1,:);
            p2 = legendre(n-2,x); p2 = p2(1,:);
            y = (p1 - p2) / sqrt(4*n-2);
        end
    end

xhat = linspace(-1, 1, Ns +1);
yhat = zeros(N+1, Ns+1);
for n = 1:N+1
    yhat(n,:) = phi(n, xhat);
end

u = zeros(1, Ns*M+1);
for m = 1:M
    Uloc = U((m-1)*N+1: m*N+1);
    u((m-1)*Ns+1: m*Ns+1) = (Uloc') * yhat;
end
end
