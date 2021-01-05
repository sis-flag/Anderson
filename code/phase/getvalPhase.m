function [u, x, hp1, hp2] = getvalPhase(U, h, Ns, N)
% get value from projection on polynomial basis (1-d ununiform mesh)
% input:
%     U(1-d array): projection on polynomial basis
%     Ns(integer):  number of sample points in each interval
%     N(integer):   degree of polynomials (default N = 10)
% output:
%     u(1-d array):     function values on sample points
%     x(1-d array):     sample points
%     hp1, hp2(number): height of the peaks


% default input
if nargin < 4
    N = 10;
end
if nargin < 3
    Ns = 50;
end

M = length(h);
xm = [0; cumsum(h)];

U = reshape(U, M*N+1, 1);

xhat = linspace(-1, 1, Ns +1);
uhat = basis(N, xhat);

u = zeros(1, Ns*M+1);
x = zeros(1, Ns*M+1);
for m = 1:M
    Uloc = U((m-1)*N+1: m*N+1);
    u((m-1)*Ns+1: m*Ns+1) = uhat * Uloc;
    x((m-1)*Ns+1: m*Ns+1) = h(m)/2 * (xhat+1) + xm(m);
end

if nargout < 3
    return
end

hp1 = max(abs( u(1*Ns+1: 2*Ns+1) ));
hp2 = max(abs( u(3*Ns+1: 6*Ns+1) ));

end