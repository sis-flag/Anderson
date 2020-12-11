%%
tic
clear
rng(0)

N = 50;
N_samp = 1000;

%% simulation
K = 1e4;
h = 0.01;
all_p = (0.05:0.05:0.95)';

Pb = zeros(length(all_p), N_samp);
Bb1 = zeros(length(all_p), N_samp, 'logical');
Bb2 = zeros(length(all_p), N_samp, 'logical');
Bb3 = zeros(length(all_p), N_samp, 'logical');
Bb4 = zeros(length(all_p), N_samp, 'logical');
for jp = 1:length(all_p)
    p = all_p(jp);
    
    for js = 1:N_samp
        V = rand(N, 1);
        V(V<p) = 0; V(V>=p) = 1;
        
        intv = [0; find(diff(V)~=0); N];
        duan = diff(intv);
        
        if V(1) == 0 && V(N) == 0
            o_o = max(duan(3:2:end-2));
            u_u = 2 * max(duan(1), duan(end));
            if isempty(o_o) || o_o < u_u
                Bb1(jp, js) = 1;
            end
        elseif V(1) == 0 && V(N) == 1
            o_o = max(duan(3:2:end-1));
            if isempty(o_o) || o_o < 2*duan(1)
                Bb2(jp, js) = 1;
            end
        elseif V(1) == 1 && V(N) == 0
            o_o = max(duan(end-2:-2:1));
            if isempty(o_o) || o_o < 2*duan(end)
                Bb3(jp, js) = 1;
            end
        elseif all(V==1)
            Bb4(jp, js) = 1;
        end
        
        U = eigR1d(K*V, h, 1);
        u = abs(getval1d(U, 20));
        Pb(jp, js) = max(u(1), u(end)) / max(u);
        
        if (Bb1(jp,js)||Bb2(jp,js)||Bb3(jp,js)) && (Pb(jp,js)<0.5)
            error('666')
        end
    end
    
    fprintf('p = %g finished\n', all_p(jp));
    toc
end

%% save
save('Data_1d_extra50.mat','all_p','K','h','N','Pb','Bb3','Bb2','Bb1')

%% predict
pB1 = zeros(length(all_p), 1);
pB2 = zeros(length(all_p), 1);
for nn = 1:length(all_p)
    p = all_p(nn);
    [tp1, tp2] = pred_prob(p, N);
    pB1(nn) = p*p*tp1;
    pB2(nn) = p*(1-p)*tp2;
end

%% plot
mPb = mean(Pb, 2);
mB1 = mean(Bb1, 2);
mB2 = mean(Bb2, 2);
mB3 = mean(Bb3, 2);
mB4 = mean(Bb4, 2);

figure
hold on
plot(all_p, mPb, 'r+--')
plot(all_p, mB1+mB2+mB3+mB4, 'mx--')
plot(all_p, pB1 + 2*pB2 + (1-all_p).^N, 'b-')
% xlabel('p'); ylabel('Pb'); ylim([0, 1]);