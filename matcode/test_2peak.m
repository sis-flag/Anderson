clear;

% Dirichlet boundary
V = 2e4 * [0,0,0,1,0,0,0,1,0,0,1];
[U, lam] = eigD1d(V, 6);
W = solveD1d(V);
% [U, lam] = eigR1d(V, 0, 6);
figure
hold on
cf = bar(((1:length(V))-0.5)/length(V), V/max(V));
cf.BarWidth = 1;
cf.LineStyle = 'None';
cf.FaceColor = 'c';
cf.FaceAlpha = 0.3;
cf = bar(((1:length(V))-0.5)/length(V), -V/max(V));
cf.BarWidth = 1;
cf.LineStyle = 'None';
cf.FaceColor = 'c';
cf.FaceAlpha = 0.3;
for k = 1:2
    uk = my_nmlz(getval1d(U(:,k)));
    x = linspace(0, 1, length(uk));
    plot(x, uk)
end
w = getval1d(W);
plot(x, w * lam(1), 'k--')
x = linspace(0, 1, length(V)+1);
plot(x, zeros(length(x),1), 'k+')
xlim([0,1])
ylim([-1,1])

disp(lam)