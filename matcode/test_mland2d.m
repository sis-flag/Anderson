%%
clear
% rng(0)

K = 1e6;
p = 0.5;
h = 10;

N = 10;
N_samp = 200;
N_eig = 6;

%% simulation
Ds = zeros(N_eig, N_samp);
for js = 1:N_samp
    V = rand(N, N);
    if p >= 0
        V(V < p) = 0; V(V >= p) = 1;
    end

    if isinf(h)
        [U, lam] = eigD2d(K*V, N_eig);
        W = solveD2d(K*V);
    else
        [U, lam] = eigR2d(K*V, h, N_eig);
        W = solveR2d(K*V, h);
    end
    
    w = abs(getval2d(W, 20));
    [wn1, wn2] = nlocmax2d(w, 0);
    wn = [wn1(1:N_eig), wn2(1:N_eig)];
    
    un = zeros(N_eig, 2);
    for je = 1:N_eig
        u = abs(getval2d(U(:,je), 20));
        [unje1, unje2] = find(u == max(max(u)));
        un(je,1) = unje1;
        un(je,2) = unje2;
    end
    Ds(:,js) = sqrt((wn(:,1)-un(:,1)).^2+(wn(:,2)-un(:,2)).^2);
    
%     if Ds(1,js) / (20*N) > 0.01
%         figure
%         x = linspace(0, 1, size(w,1));
%         y = linspace(0, 1, size(w,2));
%         subplot(1,2,1)
%         cf = surf(x,y,w);
%         cf.LineStyle = 'None';
%         subplot(1,2,2)
%         u1 = abs(getval2d(U(:,1)));
%         cf = surf(x,y,u1);
%         cf.LineStyle = 'None';
%         error('666')
%     end
end
Ds = Ds / (20*N);

%% plot mse
figure
bar(mean(Ds.^2, 2));