clear
rng(0)

all_p = [0.2, 0.5, 0.8];
all_K = [1e2, 1e3, 1e4, 1e5];
all_h = [0.1, 1, 10, 100, 1000];
N_eig = 5;
N_samp = 1000;

tic
Data_e = zeros(length(all_p), length(all_K), length(all_h), N_eig, N_samp);
Data_c = zeros(length(all_p), length(all_K), length(all_h), N_eig, N_samp);
for jp = 1:length(all_p)
    for jk = 1:length(all_K)
        for jh = 1:length(all_h)
            
            K = all_K(jk);
            p = all_p(jp);
            h = all_h(jh);
            
            % test many times and count
            for js = 1:N_samp
                
                V = rand(10);
                if p > 0
                    V(V<p) = 0; V(V>=p) = 1;
                end
                
                U = eigR2d(K*V, h, N_eig);
                
                for je = 1:N_eig
                    u = abs(getval2d(U(:,je)));
                    P_e = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(max(u));
                    P_c = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(max(u));
                    
                    Data_e(jp, jk, jh, je, js) = P_e;
                    Data_c(jp, jk, jh, je, js) = P_c;
                end
            end
            
            fprintf('K:%g, h:%g, p:%g finished\n', K, h, p);
            toc
        end
    end
end
time = toc
save('Data_2d.mat', 'all_p', 'all_K', 'all_h', 'N_eig', 'N_samp', 'Data_e', 'Data_c', 'time')
