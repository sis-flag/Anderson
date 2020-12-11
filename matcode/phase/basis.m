function [phi, dphi] = basis(N, x)

x = x(:); % x should ba a row-vector

% legendre polynomial
Lx = zeros(length(x), N);
Lx(:,1) = x;
Lx(:,2) = 1.5 * x.*x - 0.5;
for k = 3:N
    % k L_k(x) = (2k-1) x L_{k-1}(x) - (k-1) L_{k-2}(x)
    Lx(:,k) = ((2*k-1)*x.*Lx(:,k-1) - (k-1)*Lx(:,k-2)) / k;
end

% basis function
phi = zeros(length(x), N+1);
phi(:,1) = (1-x)/2;
phi(:,2) = (Lx(:,2) - 1) / sqrt(6);
phi(:,N+1) = (1+x)/2;
for n = 3:N
    phi(:,n) = (Lx(:,n) - Lx(:,n-2)) / sqrt(4*n-2);
end

if nargout < 2
    return
end

% derivative of basis function
dphi = zeros(length(x), N+1);
dphi(:,1) = -1/2;
dphi(:,N+1) = 1/2;
for n = 2:N
    dphi(:,n) = Lx(:,n-1) * sqrt(4*n-2) / 2;
end