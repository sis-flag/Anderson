%% 1-d case
clear;
rng(0);

% peremeters
K = 3000;
V =[ 1, 5, 16, 12, 4, 18, 13, 7, 14, 5,...
    11, 17, 4, 19, 5, 13, 8, 10, 6, 2] / 20;
V = rand(1,20);
V = V*K;
bar(1/40:1/20:1, V)
h = 0;

[UR, lamR] = eigR1d(V, h, 10, 50);

[UD, lamD] = eigD1d(V, 10, 50);

uRk = getval1d(UR(:,1), 50);
x = linspace(0,1,length(uRk));

%% many figures
for k = [1,2,3,5,7,10,20]
    uRk = getval1d(UR(:,k), 50);
    uRk = my_nmlz(uRk);
    uDk = getval1d(UD(:,k), 50);
    uDk = my_nmlz(uDk);
    figure();
    hold on
    plot(x, uRk);
    plot(x, uDk);
    title(sprintf('eigen mode %d', k))
    legend('Neumann boundary', 'Dirichlet boundary')
end

%% two figures
figure();
hold on
for k = [1,2,3,5,7,10,20]
    uRk = getval1d(UR(:,k), 50);
    uRk = my_nmlz(uRk);
    plot(x, uRk);
end

figure();
hold on
for k = [1,2,3,5,7,10,20]
    uDk = getval1d(UD(:,k), 50);
    uDk = my_nmlz(uDk);
    plot(x, uDk);
end

%% 2-d case
clear;
rng(0);

% peremeters
K = 3000;
V = rand(20)*K;
h = 0;

tic
[UR, lamR] = eigR2d(V, h, 10, 50);
toc
[UD, lamD] = eigD2d(V, 10, 50);
toc

u1 = getval2d(UR(:,1), 20);
x = linspace(0,1,length(u1));
[x2, x1] = meshgrid(x, x); % caution!

%% many figures
for k = [2,9]
    uRk = getval2d(UR(:,k), 20);
    uRk = my_nmlz(uRk);
    uDk = getval2d(UD(:,k), 20);
    uDk = my_nmlz(uDk);
    figure();
    subplot(1,2,1);
    mesh(x1, x2, uRk);
    set(gca, 'ZLim', [-1, 1]);
    set(gca, 'CLim', [-1, 1]);
    title('Neumann boundary')
    subplot(1,2,2);
    mesh(x1, x2, uDk);
    set(gca, 'ZLim', [-1, 1]);
    set(gca, 'CLim', [-1, 1]);
    title('Dirichlet boundary')
end
