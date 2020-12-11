%%
clear
% rng(0)

K = 1e6;
p = 0.5;
h = inf;

N = 50;
N_samp = 1000;
N_eig = 6;

%% simulation
Ds = zeros(N_eig, N_samp);
for js = 1:N_samp
    V = rand(N, 1);
    if p >= 0
        V(V < p) = 0; V(V >= p) = 1;
    end

    if isinf(h)
        [U, lam] = eigD1d(K*V, N_eig);
        W = solveD1d(K*V);
    else
        [U, lam] = eigR1d(K*V, h, N_eig);
        W = solveR1d(K*V, h);
    end
    
    w = abs(getval1d(W, 20));
    wn = nlocmax(w, 0);
    wn = wn(1:N_eig);
    
    un = zeros(N_eig, 1);
    for je = 1:N_eig
        u = abs(getval1d(U(:,je), 20));
        [~, unje] = max(u);
        un(je) = unje;
    end
    Ds(:,js) = wn - un;
    
    if abs(Ds(:,js)) / (20*N) > 0.01
        figure
        hold on
        x = linspace(0, 1, length(w));
        plot(x, w, 'k-');
        for je = 1:N_eig
            uk = getval1d(U(:,je), 20);
            uk = abs(my_nmlz(uk)) / lam(je);
            plot(x, uk)
        end
        xlim([0,1])
        error('666')
    end
end
Ds = Ds / (20*N);

%% plot
figure
boxplot(Ds')

figure
bar(mean(Ds.^2, 2));