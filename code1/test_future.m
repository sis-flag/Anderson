clear
rng(0)
tic

N_samp = 1000;
N1d = 50;
N2d = 15;

thes = 0.5;

all_h = [3e-1, 1e0, 3e0, 1e1, 3e1, 1e2, 3e2];

all_genV = cell(12,1);
all_genV{1} = @(sz) 1e4 * binornd(1, 1/2, sz);
all_genV{2} = @(sz) normrnd(5e3, 5e3, sz);
all_genV{3} = @(sz) gamrnd(1, 5e3, sz);
all_genV{4} = @(sz) (1-sqrt(3))*5e3 + 1e4*sqrt(3)*rand(sz);
all_genV{5} = @(sz) 2e4/3 * binornd(1, 3/4, sz);
all_genV{6} = @(sz) normrnd(5e3, 5e3/sqrt(3), sz);
all_genV{7} = @(sz) gamrnd(3, 5e3/3, sz);
all_genV{8} = @(sz) 1e4 * rand(sz);
all_genV{9} = @(sz) 5e4/9 * binornd(1, 0.9, sz);
all_genV{10} = @(sz) normrnd(5e3, 5e3/3, sz);
all_genV{11} = @(sz) gamrnd(9, 5e3/9, sz);
all_genV{12} = @(sz) (1-1/sqrt(3))*5e3 + 1e4/sqrt(3)*rand(sz);

%% run!

Pb_h = zeros(N_samp, length(all_h), length(all_genV), 'logical');
Pe_h = zeros(N_samp, length(all_h), length(all_genV), 'logical');
Pc_h = zeros(N_samp, length(all_h), length(all_genV), 'logical');

for jv = 9:12
    genV = all_genV{jv};
    for jh = 1:length(all_h)
        h = all_h(jh);
        for js = 1:N_samp
            V = genV([N1d, 1]);
            U = eig1d(V, h, 1);
            u = abs(getval1d(U));
            tPb = max(u(1), u(end)) / max(u) > thes;
            Pb_h(js, jh, jv) = tPb;
            
            V = genV([N2d, N2d]);
            U = eig2d(V, h, 1);
            u = abs(getval2d(U));
            tPe = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(u(:)) > thes;
            tPc = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(u(:)) > thes;
            Pe_h(js, jh, jv) = tPe;
            Pc_h(js, jh, jv) = tPc;
        end
        toc
    end
    fprintf('%d / %d V finished\n', jv, length(all_genV));
end

Pb_hp = reshape(mean(Pb_h, 1), [length(all_h), length(all_genV)]);
Pe_hp = reshape(mean(Pe_h, 1), [length(all_h), length(all_genV)]);
Pc_hp = reshape(mean(Pc_h, 1), [length(all_h), length(all_genV)]);

save('future.mat')

%% plot
load('future.mat')

for jv = 1:length(all_genV)
    figure
    hold on
    plot(all_h, Pb_hp(:,jv), '.-', 'MarkerSize', 15, 'LineWidth', 1)
    plot(all_h, Pe_hp(:,jv), '.-', 'MarkerSize', 15, 'LineWidth', 1)
    plot(all_h, Pc_hp(:,jv), '.-', 'MarkerSize', 15, 'LineWidth', 1)
    legend('P_b', 'P_e', 'P_c')
    ylim([0, 1])
    xlabel('h')
    xticks([1e-1, 1e0, 1e1, 1e2, 1e3])
    xlim([1e-1, 1e3])
    ylabel('probability')
    set(gcf, 'Position', [300 300 400 300])
    set(gca, 'FontSize', 15)
    set(gca, 'XScale', 'log')
end