%%
tic
clear
rng(0)

N_eig = 5;
N_samp = 3;

%% 1d-h
p = 0.5;
K = 1e3;
all_h = [1e-1, 3e-1, 1e0, 3e0, 1e1, 3e1, 1e2, 3e2, 1e3];

Data_hb = zeros(length(all_h), N_eig, N_samp);
for jh = 1:length(all_h)
    h = all_h(jh);
    for js = 1:N_samp
        V = rand(20, 1);
        V(V<p) = 0; V(V>=p) = 1;
        U = eigR1d(K*V, h, N_eig);
        for je = 1:N_eig
            u = abs(getval1d(U(:,je), 20));
            P_b = max(u(1), u(end)) / max(u);
            Data_hb(jh, je, js) = P_b;
        end
    end
    
    fprintf('h-1d %.2f finished\n', jh/length(all_h));
    toc
end

%% 1d-p
all_p = 0.1:0.1:0.9;
K = 1e3;
h = 0.1;

Data_pb1 = zeros(length(all_p), N_eig, N_samp);
for jp = 1:length(all_p)
    p = all_p(jp);
    for js = 1:N_samp
        V = rand(20, 1);
        V(V<p) = 0; V(V>=p) = 1;
        U = eigR1d(K*V, h, N_eig);
        for je = 1:N_eig
            u = abs(getval1d(U(:,je), 20));
            P_b = max(u(1), u(end)) / max(u);
            Data_pb1(jp, je, js) = P_b;
        end
    end
    
    fprintf('p-1d %.2f finished\n', jp/length(all_p));
    toc
end

all_p = 0.1:0.1:0.9;
K = 1e3;
h = 100;

Data_pb2 = zeros(length(all_p), N_eig, N_samp);
for jp = 1:length(all_p)
    p = all_p(jp);
    for js = 1:N_samp
        V = rand(20, 1);
        V(V<p) = 0; V(V>=p) = 1;
        U = eigR1d(K*V, h, N_eig);
        for je = 1:N_eig
            u = abs(getval1d(U(:,je), 20));
            P_b = max(u(1), u(end)) / max(u);
            Data_pb2(jp, je, js) = P_b;
        end
    end
    
    fprintf('p-1d %.2f finished\n', jp/length(all_p));
    toc
end

%% 1d-K
p = 0.5;
all_K = [1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5];
h = 0.1;

Data_kb1 = zeros(length(all_K), N_eig, N_samp);
for jk = 1:length(all_K)
    K = all_K(jk);
    for js = 1:N_samp
        V = rand(20, 1);
        V(V<p) = 0; V(V>=p) = 1;
        U = eigR1d(K*V, h, N_eig);
        for je = 1:N_eig
            u = abs(getval1d(U(:,je), 20));
            P_b = max(u(1), u(end)) / max(u);
            Data_kb1(jk, je, js) = P_b;
        end
    end
    
    fprintf('K-1d %.2f finished\n', jk/length(all_K));
    toc
end

p = 0.5;
all_K = [1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5];
h = 100;

Data_kb2 = zeros(length(all_K), N_eig, N_samp);
for jk = 1:length(all_K)
    K = all_K(jk);
    for js = 1:N_samp
        V = rand(20, 1);
        V(V<p) = 0; V(V>=p) = 1;
        U = eigR1d(K*V, h, N_eig);
        for je = 1:N_eig
            u = abs(getval1d(U(:,je), 20));
            P_b = max(u(1), u(end)) / max(u);
            Data_kb2(jk, je, js) = P_b;
        end
    end
    
    fprintf('K-1d %.2f finished\n', jk/length(all_K));
    toc
end

%% save
save('Data_1d_s.mat', 'all_p', 'all_K', 'all_h', 'N_eig', 'N_samp',...
    'Data_hb', 'Data_pb1', 'Data_pb2', 'Data_kb1', 'Data_kb2')

%%
clear
rng(0)

N_eig = 5;
N_samp = 3;

%% 2d-h
p = 0.5;
K = 1e3;
all_h = [1e-1, 3e-1, 1e0, 3e0, 1e1, 3e1, 1e2, 3e2, 1e3];

Data_he = zeros(length(all_h), N_eig, N_samp);
Data_hc = zeros(length(all_h), N_eig, N_samp);
for jh = 1:length(all_h)
    h = all_h(jh);
    for js = 1:N_samp
        V = rand(10);
        V(V<p) = 0; V(V>=p) = 1;
        U = eigR2d(K*V, h, N_eig);
        for je = 1:N_eig
            u = abs(getval2d(U(:,je)));
            P_e = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(max(u));
            P_c = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(max(u));
            
            Data_he(jh, je, js) = P_e;
            Data_hc(jh, je, js) = P_c;
        end
    end
    
    fprintf('h-2d %.2f finished\n', jh/length(all_h));
    toc
end

%% 2d-p
all_p = 0.1:0.2:0.9;
K = 1e4;
h = 100;

Data_pe1 = zeros(length(all_p), N_eig, N_samp);
Data_pc1 = zeros(length(all_p), N_eig, N_samp);
for jp = 1:length(all_p)
    p = all_p(jp);
    for js = 1:N_samp
        V = rand(10);
        V(V<p) = 0; V(V>=p) = 1;
        U = eigR2d(K*V, h, N_eig);
        for je = 1:N_eig
            u = abs(getval2d(U(:,je)));
            P_e = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(max(u));
            P_c = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(max(u));
            
            Data_pe1(jp, je, js) = P_e;
            Data_pc1(jp, je, js) = P_c;
        end
    end
    
    fprintf('p-2d %.2f finished\n', jp/length(all_p));
    toc
end

all_p = 0.1:0.2:0.9;
K = 1e4;
h = 1;

Data_pe2 = zeros(length(all_p), N_eig, N_samp);
Data_pc2 = zeros(length(all_p), N_eig, N_samp);
for jp = 1:length(all_p)
    p = all_p(jp);
    for js = 1:N_samp
        V = rand(10);
        V(V<p) = 0; V(V>=p) = 1;
        U = eigR2d(K*V, h, N_eig);
        for je = 1:N_eig
            u = abs(getval2d(U(:,je)));
            P_e = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(max(u));
            P_c = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(max(u));
            
            Data_pe2(jp, je, js) = P_e;
            Data_pc2(jp, je, js) = P_c;
        end
    end
    
    fprintf('p-2d %.2f finished\n', jp/length(all_p));
    toc
end

%% 2d-K
p = 0.5;
all_K = [1e1, 1e2, 1e3, 1e4, 1e5];
h = 1;

Data_ke1 = zeros(length(all_K), N_eig, N_samp);
Data_kc1 = zeros(length(all_K), N_eig, N_samp);
for jk = 1:length(all_K)
    K = all_K(jk);
    for js = 1:N_samp
        V = rand(10);
        V(V<p) = 0; V(V>=p) = 1;
        U = eigR2d(K*V, h, N_eig);
        for je = 1:N_eig
            u = abs(getval2d(U(:,je)));
            P_e = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(max(u));
            P_c = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(max(u));
            
            Data_ke1(jk, je, js) = P_e;
            Data_kc1(jk, je, js) = P_c;
        end
    end
    
    fprintf('k-2d %.2f finished\n', jk/length(all_K));
    toc
end


p = 0.5;
all_K = [1e1, 1e2, 1e3, 1e4, 1e5];
h = 100;

Data_ke2 = zeros(length(all_K), N_eig, N_samp);
Data_kc2 = zeros(length(all_K), N_eig, N_samp);
for jk = 1:length(all_K)
    K = all_K(jk);
    for js = 1:N_samp
        V = rand(10);
        V(V<p) = 0; V(V>=p) = 1;
        U = eigR2d(K*V, h, N_eig);
        for je = 1:N_eig
            u = abs(getval2d(U(:,je)));
            P_e = max([u(1,:), u(end,:), u(:,1)', u(:,end)']) / max(max(u));
            P_c = max([u(1,1), u(1,end), u(end,1), u(end,end)]) / max(max(u));
            
            Data_ke2(jk, je, js) = P_e;
            Data_kc2(jk, je, js) = P_c;
        end
    end
    
    fprintf('K-2d %.2f finished\n', jk/length(all_K));
    toc
end

%% save
save('Data_2d_s.mat', 'all_p', 'all_K', 'all_h', 'N_eig', 'N_samp',...
    'Data_hc', 'Data_pc1', 'Data_pc2', 'Data_kc1', 'Data_kc2',...
    'Data_he', 'Data_pe1', 'Data_pe2', 'Data_ke1', 'Data_ke2')
