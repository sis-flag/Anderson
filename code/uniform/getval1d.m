function u = getval1d(U, Ns, N)
% get value from projection on polynomial basis (1-d)
% input:
%     U(1-d array): projection on polynomial basis
%     Ns(integer):  number of sample points in each interval
%     N(integer):   degree of polynomials (default N = 10)
% output:
%     u(1-d array): function values on sample points

% default input
if nargin < 3
    N = 10;
end
if nargin < 2
    Ns = 50;
end

M = round(length(U)-1) / N;

U = reshape(U, M*N+1, 1);

xhat = linspace(-1, 1, Ns +1);
uhat = basis(N, xhat);

u = zeros(1, Ns*M+1);
for m = 1:M
    Uloc = U((m-1)*N+1: m*N+1);
    u((m-1)*Ns+1: m*Ns+1) = uhat * Uloc;
end

end