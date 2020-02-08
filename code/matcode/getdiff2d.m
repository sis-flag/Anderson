% function [du1, du2] = getdiff2d(U, Ns, N)
% 
% if nargin <= 2
%     N = 10;
% end
% M = round(sqrt(length(U))-1) / N;
% 
%     function y = phi(n, x)
%         if n == 1
%             y = (1-x)/2;
%         elseif n == N+1
%             y = (1+x)/2;
%         else
%             p1 = legendre(n,x); p1 = p1(1,:,:);
%             p2 = legendre(n-2,x); p2 = p2(1,:,:);
%             y = (p1 - p2) / sqrt(4*n-2);
%         end
%     end
% 
%     function y = dphi(n, x)
%         if n == 1
%             y = -1/2;
%         elseif n == N+1
%             y = 1/2;
%         else
%             p1 = legendre(n-1,x);
%             y = p1(1,:,:);
%         end
%     end
% 
% xhat = linspace(-1, 1, Ns +1);
% [xhat2, xhat1] = meshgrid(xhat, xhat);
% yhat1 = zeros(Ns+1, Ns+1, N+1);
% yhat2 = zeros(Ns+1, Ns+1, N+1);
% dyhat1 = zeros(Ns+1, Ns+1, N+1);
% dyhat2 = zeros(Ns+1, Ns+1, N+1);
% for n = 1:N+1
%     yhat1(:,:,n) = phi(n, xhat1);
%     yhat2(:,:,n) = phi(n, xhat2);
%     dyhat1(:,:,n) = 2*M * dphi(n, xhat1);
%     dyhat2(:,:,n) = 2*M * dphi(n, xhat2);
% end
% 
% U = reshape(U, M*N+1, M*N+1);
% 
% du1 = zeros(1, Ns*M+1);
% du2 = zeros(1, Ns*M+1);
% for m1 = 1:M
%     for m2 = 1:M
%         Uloc = U((m1-1)*N+1: m1*N+1, (m2-1)*N+1: m2*N+1);
%         du1((m1-1)*Ns+1: m1*Ns+1, (m2-1)*Ns+1: m2*Ns+1) = 0;
%         du2((m1-1)*Ns+1: m1*Ns+1, (m2-1)*Ns+1: m2*Ns+1) = 0;
%         for n1 = 1:N+1
%             for n2 = 1:N+1
%                 du1((m1-1)*Ns+1: m1*Ns+1, (m2-1)*Ns+1: m2*Ns+1) = ...
%                     du1((m1-1)*Ns+1: m1*Ns+1, (m2-1)*Ns+1: m2*Ns+1)...
%                     + Uloc(n1, n2) * dyhat1(:,:,n1) .* yhat2(:,:,n2);
%                 du2((m1-1)*Ns+1: m1*Ns+1, (m2-1)*Ns+1: m2*Ns+1) = ...
%                     du1((m1-1)*Ns+1: m1*Ns+1, (m2-1)*Ns+1: m2*Ns+1)...
%                     + Uloc(n1, n2) * yhat1(:,:,n1) .* dyhat2(:,:,n2);
%             end
%         end
%     end
% end
% 
% end
