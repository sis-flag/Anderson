clear;

%% example 1
rng(0);

% peremeters
K = 8000;
V = rand(1, 20);
V = K * V;

% plot potential
figure();
cf = bar((1:20)-0.5, V);
cf.BarWidth = 1;
saveas(gcf, '../report0/parameter/V1.png')

[U, lam] = eigR1d(V, 0, 6);

for x0 = [0.3, 0.5, 0.7]
    for alpha = [1000, 3000, 10000]
        W = solveN1d(V, 0, x0, 1/alpha);

        w = getval1d(W);
        x = linspace(0, 1, length(w));

        figure();
        hold on
        plot(x, w, 'k')
        for k = 1:6
            uk = my_nmlz(getval1d(U(:,k)));
            plot(x, uk / (lam(k) + alpha) )
        end
        title(sprintf('\\alpha = %g', alpha))
        saveas(gcf, sprintf('../report0/parameter/W1X%dA%d.png', x0*10, alpha))
    end
end

figure();
hold on
for x0 = [0.3, 0.5, 0.7]
    for alpha = [1000, 3000, 10000]
        W = solveN1d(V, 0, x0, 1/alpha);

        w = getval1d(W);
        x = linspace(0, 1, length(w));

        plot(x, w)
        saveas(gcf, '../report0/parameter/ALL1.png')
    end
end

figure();
hold on
alpha = 300;
for x0 = 0.1: 0.1: 0.9
    W = solveN1d(V, 0, x0, 1/alpha);

    w = getval1d(W);
    x = linspace(0, 1, length(w));

    plot(x, w)
end
W = solveN1d(V, 0);
w = getval1d(W);
x = linspace(0, 1, length(w));
plot(x, w, 'k')
saveas(gcf, '../report0/parameter/INF1.png')

%% example2
rng(0);

% peremeters
K = 8000;
V = rand(1, 20);
p = 0.3;
V(V<p) = 0; V(V>=p) = 1;
V = K * V;

% plot potential
figure();
cf = bar((1:20)-0.5, V);
cf.BarWidth = 1;
saveas(gcf, '../report0/parameter/V2.png')

[U, lam] = eigR1d(V, 0, 6);

for x0 = [0.3, 0.5, 0.7]
    for alpha = [1000, 3000, 10000]
        W = solveN1d(V, 0, x0, 1/alpha);

        w = getval1d(W);
        x = linspace(0, 1, length(w));

        figure();
        hold on
        plot(x, w, 'k')
        for k = 1:6
            uk = my_nmlz(getval1d(U(:,k)));
            plot(x, uk / (lam(k) + alpha) )
        end
        title(sprintf('\\alpha = %g', alpha))
        saveas(gcf, sprintf('../report0/parameter/W2X%dA%d.png', x0*10, alpha))
    end
end

figure();
hold on
for x0 = [0.3, 0.5, 0.7]
    for alpha = [1000, 3000, 10000]
        W = solveN1d(V, 0, x0, 1/alpha);

        w = getval1d(W);
        x = linspace(0, 1, length(w));

        plot(x, w)
        saveas(gcf, '../report0/parameter/ALL2.png')
    end
end

figure();
hold on
alpha = 300;
for x0 = 0.1: 0.1: 0.9
    W = solveN1d(V, 0, x0, 1/alpha);

    w = getval1d(W);
    x = linspace(0, 1, length(w));

    plot(x, w)
end
W = solveN1d(V, 0);
w = getval1d(W);
x = linspace(0, 1, length(w));
plot(x, w, 'k')
saveas(gcf, '../report0/parameter/INF2.png')

%% example 2d alpha
rng(0);

K = 8000;
V = rand(20);
V = K * V;

W = solveN2d(V, 0);

w = getval2d(W);
x = linspace(0, 1, size(w,1));
[x2, x1] = meshgrid(x, x); % caution!

figure()
mesh(x1, x2, w);
saveas(gcf, '../report0/parameter/2DA1.png')

W = solveN2d(V, 0, 0.5, 0.5, 1/700);
w = getval2d(W);

figure()
mesh(x1, x2, w);
saveas(gcf, '../report0/parameter/2DA2.png')
