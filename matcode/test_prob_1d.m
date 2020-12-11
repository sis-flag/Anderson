clear
rng(0)

all_p = [0.2, 0.5, 0.8];
all_K = [1e2, 1e3, 1e4, 1e5];
all_h = [0.01, 0.1, 1, 10, 100];
N_eig = 5;
N_samp = 1000;

tic
Data = zeros(length(all_p), length(all_K), length(all_h), N_eig, N_samp);
for jp = 1:length(all_p)
    for jk = 1:length(all_K)
        for jh = 1:length(all_h)
            
            K = all_K(jk);
            p = all_p(jp);
            h = all_h(jh);
            
            % test many times and count
            for js = 1:N_samp
                
                V = rand(20, 1);
                if p > 0
                    V(V<p) = 0; V(V>=p) = 1;
                end
                
                U = eigR1d(K*V, h, N_eig);
                
                for je = 1:N_eig
                    u = abs(getval1d(U(:,je), 20));
                    P_b = max(u(1), u(end)) / max(u);
                    
                    Data(jp, jk, jh, je, js) = P_b;
                end
            end
            
            fprintf('K:%g, h:%g, p:%g finished\n', K, h, p);
            toc
        end
    end
end
time = toc
save('Data_1d.mat', 'all_p', 'all_K', 'all_h', 'N_eig', 'N_samp', 'Data', 'time')
