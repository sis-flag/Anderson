function u = getval2d(U, Ns, N)

if nargin <= 2
    N = 10;
end
M = round(sqrt(length(U))-1) / N;

    function y = phi(n, x)
        if n == 1
            y = (1-x)/2;
        elseif n == N+1
            y = (1+x)/2;
        else
            p1 = legendre(n,x); p1 = p1(1,:,:);
            p2 = legendre(n-2,x); p2 = p2(1,:,:);
            y = (p1 - p2) / sqrt(4*n-2);
        end
    end

xhat = linspace(-1, 1, Ns +1);
yhat = zeros(N+1, Ns+1);
for n = 1:N+1
    yhat(n,:) = phi(n, xhat);
end

U = reshape(U, M*N+1, M*N+1);

u = zeros(1, Ns*M+1);
for m1 = 1:M
    for m2 = 1:M
        Uloc = U((m1-1)*N+1: m1*N+1, (m2-1)*N+1: m2*N+1);
        u((m1-1)*Ns+1: m1*Ns+1, (m2-1)*Ns+1: m2*Ns+1) =...
            (yhat') * Uloc * yhat;
    end
end

end
