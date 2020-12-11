clear;

N = [4, 16, 3, 1, 3];
V = [ones(ceil(N(2)/2), 1);
    zeros(N(1),1);...
    ones(N(2),1);...
    zeros(N(3),1);...
    ones(N(4),1);...
    zeros(N(5),1);...
    ones(floor(N(2)/2),1)];

LN = [ceil(N(2)/2), ceil(N(2)/2)+N(1)];
RN = [ceil(N(2)/2)+N(1)+N(2), sum(N)+N(2)-floor(N(2)/2)];

all_K = 400;

% % plot potential
% figure();
% hold on
% xv = ((1:length(V))-0.5) / length(V);
% cf = bar(xv, V);
% cf.BarWidth = 1;
% cf.LineStyle = 'None';

all_F0 = zeros(1, length(all_K));
all_F1 = zeros(1, length(all_K));
all_F2 = zeros(1, length(all_K));
all_lam = zeros(6, length(all_K));

for jk = 1:length(all_K)
    K = all_K(jk);

    [U, lam] = eigP1d(K*V, 6);
    W = solveP1d(K*V);
    
    all_lam(:,jk) = lam;
    
    u1 = my_nmlz(getval1d(U(:,1), 20));
    left1 = max(abs(u1(LN(1)*20+1: LN(2)*20+1)));
    right1 = max(abs(u1(RN(1)*20+1: RN(2)*20+1)));
    all_F1(jk) = left1 / (left1 + right1);
    
    u2 = my_nmlz(getval1d(U(:,2), 20));
    u2 = u2 * sign(u2(RN(2)*20));
    left2 = max(abs(u2(LN(1)*20+1: LN(2)*20+1)));
    right2 = max(abs(u2(RN(1)*20+1: RN(2)*20+1)));
    all_F2(jk) = left2 / (left2 + right2);
    
    w = getval1d(W, 20);
    left0 = max(w(LN(1)*20+1: LN(2)*20+1));
    right0 = max(w(RN(1)*20+1: RN(2)*20+1));
    all_F0(jk) = left0 / (left0 + right0);
    
    x = linspace(0, 1, length(w));
    
    figure(1);
    hold on;
    plot(x, u1);
    set(gcf,'position',[300,300,400,300])
    
    figure(2);
    hold on;
    plot(x, u2);
    set(gcf,'position',[300,300,400,300])
    
    figure(3);
    hold on;
    plot(x, w);
    set(gcf,'position',[300,300,400,300])
end

figure(4)
hold on
plot(all_K, all_lam, '*-')
set(gcf,'position',[300,300,400,300])

figure(5)
plot(all_K, all_F0, '*-')
set(gcf,'position',[300,300,400,300])

figure(6)
plot(all_K, all_F1, '*-')
ylim([0, 1])
set(gcf,'position',[300,300,400,300])

figure(7)
plot(all_K, all_F2, '*-')
ylim([0, 1])
set(gcf,'position',[300,300,400,300])
