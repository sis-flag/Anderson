clear;

%% 1d potential
rng(0);

K = 8000;
% V = rand(1, 20);
V =[ 1, 5, 16, 12, 4, 18, 13, 7, 14, 5,...
    11, 17, 4, 19, 5, 13, 8, 10, 6, 2] / 20;
V = K * V;

% plot potential
% figure();
% cf = bar(1/40:1/20:1, V);
% cf.BarWidth = 1;
% saveas(gcf, '../report0/soul/V1d.png')

%% 1d soul Robin boundary
h = 10;
beta = -200;

[U, lam] = eigD1d(V, 6);
w1 = getval1d(solveD1d(V));
w2 = getval1d(solveD1d(V));
x = linspace(0, 1, length(w1));

% figure();
subplot(1,2,1)
hold on
plot(x, w1, 'k')
for k = 1:4
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, uk / lam(k) )
end
title('landscape 1')

% figure();
subplot(1,2,2)
hold on
plot(x, w2, 'k')
for k = 1:4
    uk = my_nmlz(getval1d(U(:,k)));
    plot(x, uk / (lam(k) + beta) )
end
title(sprintf('landscape 2 (\\beta = %g)', beta))

%% test correction
% p = 0.7;
% for h = [1,10,100,1e3,1e4]
%     for K = [1,10,1e2,1e3,1e4,1e5]
%         V = rand(1, 20);
%         V(V<p) = 0; V(V>p) = 1;
%         V = K * V;
%         
%         [U1, lam] = eigR1d(V, h, 20);
%         w = getval1d(solveR1d(V, h));
%         
%         for k = 1:20
%             u1k = my_nmlz(getval1d(U1(:,k)));
%             if any(w < u1k / lam(k) )
%                 disp('wrong!!!')
%                 disp([h,K,k])
%                 
%                 x = linspace(0,1,length(w));
%                 
%                 figure();
%                 hold on
%                 plot(x, w, 'k')
%                 plot(x, u1k / lam(k) )
%                 return
%             end
%         end
%     end
% end