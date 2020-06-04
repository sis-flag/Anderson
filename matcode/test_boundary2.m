clear;

%% 1d potential
rng(0);

K = 8000;
V = rand(1, 20);
V = K * V;

%% 1d compare
w0 = getval1d(solveD1d(V));
w1 = getval1d(solveR1d(V, 0.03));
w2 = getval1d(solveR1d(V, 0));
w3 = getval1d(solveR1d(V, -0.03));

x = linspace(0, 1, length(w1));

figure();
hold on
plot(x, w0)
plot(x, w1)
plot(x, w2)
plot(x, w3)
legend('Dirichlet',...
    'h/\beta=0.03',...
    'h/\beta=0',...
    'h/\beta=-0.03'...
)
saveas(gcf, '../report0/boundary/C1d.png')

%% 2d potential
rng(0);

K = 8000;
V = rand(20);
V = K * V;

%% 2d compare
w0 = getval2d(solveD2d(V));
w1 = getval2d(solveR2d(V, 0.03));
w2 = getval2d(solveR2d(V, 0));
w3 = getval2d(solveR2d(V, -0.03));

x = linspace(0, 1, size(w0,1));

for k = 0:5
    figure();
    hold on
    plot(x, w0(:,k*80+1))
    plot(x, w1(:,k*80+1))
    plot(x, w2(:,k*80+1))
    plot(x, w3(:,k*80+1))
    title(sprintf('y = %g', k/5))
    legend('Dirichlet', 'h/\beta=0.03', 'h/\beta=0', 'h/\beta=-0.03')
    saveas(gcf, sprintf('../report0/boundary/C2d(%d).png',k))
end